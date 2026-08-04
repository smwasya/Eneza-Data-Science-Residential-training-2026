# Beta diversity analysis by disease status
#
# Bray-Curtis dissimilarity, PCoA, 95% confidence ellipses, PERMANOVA,
# and a homogeneity-of-multivariate-dispersion test.
#
# Required packages: phyloseq, vegan, ggplot2

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
})

# ---- Configuration ----------------------------------------------------------
phyloseq_file <- NA_character_  # Leave as NA to use loaded `ps` or a file picker
metadata_file <- NA_character_  # Example: "sample_metadata.csv"
sample_id_column <- "SampleID"  # Used only when metadata_file is supplied
group_column <- "DISEASE"
output_dir <- "beta_diversity_results"
permutation_count <- 999
minimum_prevalence_samples <- 3
set.seed(128)

# ---- Validation and metadata ------------------------------------------------
if (exists("ps", inherits = TRUE) && inherits(get("ps", inherits = TRUE), "phyloseq")) {
  ps <- get("ps", inherits = TRUE)
  message("Using the phyloseq object `ps` already loaded in R.")
} else {
  if (is.na(phyloseq_file) || !nzchar(phyloseq_file)) {
    message("Select the phyloseq .rds file in the file-selection window.")
    phyloseq_file <- file.choose()
  }
  if (!file.exists(phyloseq_file)) {
    stop("Phyloseq file not found: ", normalizePath(phyloseq_file, mustWork = FALSE))
  }
  ps <- readRDS(phyloseq_file)
}

if (!inherits(ps, "phyloseq")) {
  stop("The RDS file does not contain a phyloseq object.")
}

if (!is.na(metadata_file)) {
  if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)
  metadata_external <- read.csv(metadata_file, check.names = FALSE)
  required_columns <- c(sample_id_column, group_column)
  missing_columns <- setdiff(required_columns, names(metadata_external))
  if (length(missing_columns)) {
    stop("Metadata is missing: ", paste(missing_columns, collapse = ", "))
  }
  if (anyDuplicated(as.character(metadata_external[[sample_id_column]]))) {
    stop("Sample IDs in the metadata file must be unique.")
  }

  positions <- match(sample_names(ps), as.character(metadata_external[[sample_id_column]]))
  if (anyNA(positions)) {
    missing_ids <- sample_names(ps)[is.na(positions)]
    stop(
      "Metadata could not be matched to ", length(missing_ids),
      " phyloseq sample(s). Examples: ", paste(head(missing_ids, 5), collapse = ", ")
    )
  }
  sample_data(ps)[[group_column]] <- metadata_external[[group_column]][positions]
}

if (!group_column %in% sample_variables(ps)) {
  stop(
    "'", group_column, "' is absent from sample_data(ps). ",
    "Add it to the phyloseq object or configure metadata_file."
  )
}

group_values <- as.character(data.frame(sample_data(ps))[[group_column]])
keep <- !is.na(group_values) & trimws(group_values) != "" & sample_sums(ps) > 0
if (sum(!keep)) message("Removing ", sum(!keep), " sample(s) with missing group or zero reads.")

ps_beta <- prune_samples(keep, ps)
ps_beta <- prune_taxa(taxa_sums(ps_beta) > 0, ps_beta)

# ---- Filter low-prevalence and unwanted taxa --------------------------------
# Retain taxa observed in at least `minimum_prevalence_samples` samples.
# The orientation check avoids assuming that taxa are stored as columns.
taxa_before_prevalence <- ntaxa(ps_beta)
otu_matrix <- as(otu_table(ps_beta), "matrix")
taxon_prevalence <- if (taxa_are_rows(ps_beta)) {
  rowSums(otu_matrix > 0)
} else {
  colSums(otu_matrix > 0)
}

ps_beta <- prune_taxa(
  taxon_prevalence >= minimum_prevalence_samples,
  ps_beta
)

cat("Taxa before prevalence filtering:", taxa_before_prevalence, "\n")
cat("Taxa after prevalence filtering :", ntaxa(ps_beta), "\n")

# Remove mitochondrial and chloroplast sequences. Unclassified taxa are kept.
available_ranks <- rank_names(ps_beta)
if (!all(c("Family", "Order") %in% available_ranks)) {
  stop(
    "Removing mitochondria/chloroplasts requires Family and Order ranks. ",
    "Available ranks: ", paste(available_ranks, collapse = ", ")
  )
}

taxonomy <- as(tax_table(ps_beta), "matrix")
keep_non_organelle <- (
  is.na(taxonomy[, "Family"]) |
    taxonomy[, "Family"] != "Mitochondria"
) & (
  is.na(taxonomy[, "Order"]) |
    taxonomy[, "Order"] != "Chloroplast"
)

taxa_before_organelle_filter <- ntaxa(ps_beta)
ps_beta <- prune_taxa(keep_non_organelle, ps_beta)
ps_beta <- prune_taxa(taxa_sums(ps_beta) > 0, ps_beta)

cat("Taxa before organelle filtering:", taxa_before_organelle_filter, "\n")
cat("Taxa after organelle filtering :", ntaxa(ps_beta), "\n")

if (ntaxa(ps_beta) == 0) {
  stop("No taxa remain after prevalence and organelle filtering.")
}

# Filtering can rarely leave a sample with no reads; remove such samples before
# converting counts to proportions.
zero_read_samples <- sample_sums(ps_beta) == 0
if (any(zero_read_samples)) {
  message("Removing ", sum(zero_read_samples), " sample(s) with zero reads after filtering.")
  ps_beta <- prune_samples(!zero_read_samples, ps_beta)
}

# Relative abundance is used here to make library sizes comparable.
ps_relative <- transform_sample_counts(ps_beta, function(x) x / sum(x))

# ---- Bray-Curtis distance and PCoA ------------------------------------------
bray_distance <- phyloseq::distance(ps_relative, method = "bray")
pcoa <- ordinate(ps_relative, method = "PCoA", distance = bray_distance)

metadata <- data.frame(sample_data(ps_relative), check.names = FALSE)
metadata <- metadata[labels(bray_distance), , drop = FALSE]
metadata[[group_column]] <- droplevels(factor(metadata[[group_column]]))

stopifnot(
  !anyNA(metadata[[group_column]]),
  nlevels(metadata[[group_column]]) >= 2,
  identical(rownames(metadata), labels(bray_distance))
)

pcoa_df <- data.frame(
  SampleID = rownames(pcoa$vectors),
  Axis1 = pcoa$vectors[, 1],
  Axis2 = pcoa$vectors[, 2],
  stringsAsFactors = FALSE
)
pcoa_df[[group_column]] <- metadata[
  match(pcoa_df$SampleID, rownames(metadata)),
  group_column
]
pcoa_df[[group_column]] <- droplevels(factor(pcoa_df[[group_column]]))

relative_eigenvalues <- pcoa$values$Relative_eig
axis1_percent <- round(100 * relative_eigenvalues[1], 1)
axis2_percent <- round(100 * relative_eigenvalues[2], 1)

# ---- PERMANOVA ---------------------------------------------------------------
set.seed(128)
permanova <- adonis2(
  reformulate(group_column, response = "bray_distance"),
  data = metadata,
  permutations = permutation_count,
  by = "margin"
)
permanova_table <- data.frame(term = rownames(permanova), permanova, row.names = NULL)

# ---- Homogeneity of multivariate dispersion ---------------------------------
dispersion <- betadisper(bray_distance, group = metadata[[group_column]])
dispersion_anova <- anova(dispersion)
set.seed(128)
dispersion_permutation <- permutest(dispersion, permutations = permutation_count)

dispersion_anova_table <- data.frame(
  term = rownames(dispersion_anova), dispersion_anova, row.names = NULL
)
dispersion_permutation_table <- data.frame(
  term = rownames(dispersion_permutation$tab),
  dispersion_permutation$tab,
  row.names = NULL
)

# Record variance explained by every PCoA axis.
pcoa_variance_table <- data.frame(
  Axis = seq_along(relative_eigenvalues),
  Relative_eigenvalue = relative_eigenvalues,
  Variance_explained_percent = 100 * relative_eigenvalues
)

# Distances to group centroids are the values tested by PERMDISP.
dispersion_ids <- names(dispersion$distances)
if (is.null(dispersion_ids)) dispersion_ids <- labels(bray_distance)
dispersion_values <- data.frame(
  SampleID = dispersion_ids,
  Group = metadata[
    match(dispersion_ids, rownames(metadata)),
    group_column
  ],
  Distance_to_centroid = as.numeric(dispersion$distances),
  stringsAsFactors = FALSE
)
dispersion_group_summary <- do.call(
  rbind,
  lapply(split(dispersion_values$Distance_to_centroid, dispersion_values$Group), function(x) {
    data.frame(
      n = length(x),
      mean = mean(x),
      SD = sd(x),
      median = median(x),
      IQR = IQR(x)
    )
  })
)
dispersion_group_summary[[group_column]] <- rownames(dispersion_group_summary)
rownames(dispersion_group_summary) <- NULL
dispersion_group_summary <- dispersion_group_summary[
  ,
  c(group_column, "n", "mean", "SD", "median", "IQR"),
  drop = FALSE
]

# ---- Publication-quality PCoA plot ------------------------------------------
group_sizes <- table(pcoa_df[[group_column]])
can_draw_ellipses <- all(group_sizes >= 3)
if (!can_draw_ellipses) {
  warning("Ellipses omitted because at least one group has fewer than three samples.")
}

p_beta <- ggplot(
  pcoa_df,
  aes(x = Axis1, y = Axis2,
      colour = .data[[group_column]], fill = .data[[group_column]])
)

if (can_draw_ellipses) {
  p_beta <- p_beta +
    stat_ellipse(
      aes(group = .data[[group_column]]),
      geom = "polygon", type = "t", level = 0.95,
      alpha = 0.15, colour = NA, show.legend = FALSE
    ) +
    stat_ellipse(
      aes(group = .data[[group_column]]),
      geom = "path", type = "t", level = 0.95,
      linewidth = 0.9, linetype = 2, show.legend = FALSE
    )
}

p_beta <- p_beta +
  geom_point(shape = 21, size = 3.2, stroke = 0.8, alpha = 0.90) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = paste0("PCoA 1 (", axis1_percent, "%)"),
    y = paste0("PCoA 2 (", axis2_percent, "%)"),
    colour = "Disease status",
    fill = "Disease status",
    title = "Beta diversity by disease status",
    subtitle = "Bray-Curtis dissimilarity; 95% confidence ellipses"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, colour = "black"),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    plot.margin = margin(10, 15, 10, 10)
  )

# ---- Save outputs ------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(pcoa_df, file.path(output_dir, "bray_curtis_pcoa_coordinates.csv"), row.names = FALSE)
write.csv(permanova_table, file.path(output_dir, "permanova_results.csv"), row.names = FALSE)
write.csv(
  pcoa_variance_table,
  file.path(output_dir, "pcoa_variance_explained.csv"),
  row.names = FALSE
)
write.csv(
  dispersion_anova_table,
  file.path(output_dir, "dispersion_anova_results.csv"),
  row.names = FALSE
)
write.csv(
  dispersion_permutation_table,
  file.path(output_dir, "dispersion_permutation_results.csv"),
  row.names = FALSE
)
write.csv(
  dispersion_values,
  file.path(output_dir, "dispersion_distance_to_centroid_by_sample.csv"),
  row.names = FALSE
)
write.csv(
  dispersion_group_summary,
  file.path(output_dir, "dispersion_group_summary.csv"),
  row.names = FALSE
)
beta_summary_text <- capture.output({
  cat("BETA-DIVERSITY STATISTICAL SUMMARY\n")
  cat("Generated:", format(Sys.time()), "\n")
  cat("Distance: Bray-Curtis\n")
  cat("Permutations:", permutation_count, "\n\n")
  cat("Samples per group\n")
  print(table(metadata[[group_column]]))
  cat("\nPCoA variance explained\n")
  print(head(pcoa_variance_table, 10), row.names = FALSE)
  cat("\nPERMANOVA\n")
  print(permanova)
  cat("\nHomogeneity of dispersion: parametric ANOVA\n")
  print(dispersion_anova)
  cat("\nHomogeneity of dispersion: permutation test\n")
  print(dispersion_permutation)
  cat("\nDistance-to-centroid summary\n")
  print(dispersion_group_summary, row.names = FALSE)
})
writeLines(
  beta_summary_text,
  file.path(output_dir, "beta_statistical_summary.txt")
)
saveRDS(bray_distance, file.path(output_dir, "bray_curtis_distance.rds"))
ggsave(
  file.path(output_dir, "beta_diversity_bray_pcoa.png"),
  p_beta, width = 7, height = 5.5, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(output_dir, "beta_diversity_bray_pcoa.pdf"),
  p_beta, width = 7, height = 5.5, units = "in"
)

print(table(metadata[[group_column]]))
print(permanova)
print(dispersion_anova)
print(dispersion_permutation)
print(p_beta)
message("Beta-diversity outputs saved in: ", normalizePath(output_dir))
