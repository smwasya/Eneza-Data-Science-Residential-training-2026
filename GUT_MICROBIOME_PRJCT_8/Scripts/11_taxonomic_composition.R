# Descriptive taxonomic composition at phylum and genus levels
#
# Creates individual-sample and group-mean relative-abundance plots.
# Required packages: phyloseq, ggplot2

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
})

# ---- Configuration ----------------------------------------------------------
phyloseq_file <- NA_character_  # Use loaded `ps`; otherwise open file picker
metadata_file <- NA_character_  # Example: "sample_metadata.csv"
sample_id_column <- "SampleID"
group_column <- "DISEASE"
taxonomic_ranks <- c("Phylum", "Genus")
top_taxa_to_show <- c(Phylum = 8, Genus = 15)
top_genus_grouped_bars <- 15
minimum_prevalence_samples <- 3
output_dir <- "taxonomic_composition_results"

# ---- Load phyloseq object ----------------------------------------------------
if (exists("ps", inherits = TRUE) && inherits(get("ps", inherits = TRUE), "phyloseq")) {
  ps <- get("ps", inherits = TRUE)
  message("Using the phyloseq object `ps` already loaded in R.")
} else {
  if (is.na(phyloseq_file) || !nzchar(phyloseq_file)) {
    message("Select the phyloseq .rds file.")
    phyloseq_file <- file.choose()
  }
  if (!file.exists(phyloseq_file)) stop("Phyloseq file not found: ", phyloseq_file)
  ps <- readRDS(phyloseq_file)
}
if (!inherits(ps, "phyloseq")) stop("Input is not a phyloseq object.")

# Optionally attach external metadata using exact sample-ID matching.
if (!is.na(metadata_file)) {
  metadata_external <- read.csv(metadata_file, check.names = FALSE)
  required <- c(sample_id_column, group_column)
  if (!all(required %in% names(metadata_external))) {
    stop("Metadata must contain: ", paste(required, collapse = ", "))
  }
  positions <- match(sample_names(ps), as.character(metadata_external[[sample_id_column]]))
  if (anyNA(positions)) stop("Some phyloseq sample IDs do not match the metadata.")
  ps_metadata <- data.frame(sample_data(ps), check.names = FALSE)
  ps_metadata[[group_column]] <- metadata_external[[group_column]][positions]
  rownames(ps_metadata) <- sample_names(ps)
  sample_data(ps) <- sample_data(ps_metadata)
}

if (!group_column %in% sample_variables(ps)) {
  stop("'", group_column, "' is absent from sample_data(ps).")
}
missing_ranks <- setdiff(taxonomic_ranks, rank_names(ps))
if (length(missing_ranks)) {
  stop("Missing taxonomy ranks: ", paste(missing_ranks, collapse = ", "))
}

# ---- Sample and taxon filtering ---------------------------------------------
groups <- as.character(data.frame(sample_data(ps))[[group_column]])
keep_samples <- !is.na(groups) & trimws(groups) != "" & sample_sums(ps) > 0
ps_filtered <- prune_samples(keep_samples, ps)
ps_filtered <- prune_taxa(taxa_sums(ps_filtered) > 0, ps_filtered)

otu_matrix <- as(otu_table(ps_filtered), "matrix")
prevalence <- if (taxa_are_rows(ps_filtered)) {
  rowSums(otu_matrix > 0)
} else {
  colSums(otu_matrix > 0)
}
cat("Taxa before prevalence filtering:", ntaxa(ps_filtered), "\n")
ps_filtered <- prune_taxa(prevalence >= minimum_prevalence_samples, ps_filtered)
cat("Taxa after prevalence filtering :", ntaxa(ps_filtered), "\n")

# Remove mitochondria and chloroplasts when these ranks are available.
taxonomy <- as(tax_table(ps_filtered), "matrix")
keep_taxa <- rep(TRUE, ntaxa(ps_filtered))
if ("Family" %in% colnames(taxonomy)) {
  keep_taxa <- keep_taxa & (
    is.na(taxonomy[, "Family"]) | taxonomy[, "Family"] != "Mitochondria"
  )
}
if ("Order" %in% colnames(taxonomy)) {
  keep_taxa <- keep_taxa & (
    is.na(taxonomy[, "Order"]) | taxonomy[, "Order"] != "Chloroplast"
  )
}
ps_filtered <- prune_taxa(keep_taxa, ps_filtered)
ps_filtered <- prune_taxa(taxa_sums(ps_filtered) > 0, ps_filtered)

zero_samples <- sample_sums(ps_filtered) == 0
if (any(zero_samples)) ps_filtered <- prune_samples(!zero_samples, ps_filtered)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Composition function ----------------------------------------------------
make_composition_plots <- function(ps_input, rank_name, top_n) {
  ps_rank <- tax_glom(ps_input, taxrank = rank_name, NArm = FALSE)
  ps_relative <- transform_sample_counts(ps_rank, function(x) 100 * x / sum(x))
  composition <- psmelt(ps_relative)

  composition$Taxon <- as.character(composition[[rank_name]])
  composition$Taxon[
    is.na(composition$Taxon) | trimws(composition$Taxon) == ""
  ] <- paste("Unclassified", rank_name)
  composition$Group <- factor(composition[[group_column]])

  mean_by_taxon <- aggregate(Abundance ~ Taxon, composition, mean)
  top_taxa <- head(
    mean_by_taxon$Taxon[order(mean_by_taxon$Abundance, decreasing = TRUE)],
    top_n
  )
  composition$Taxon_display <- ifelse(
    composition$Taxon %in% top_taxa,
    composition$Taxon,
    "Other"
  )

  # Sum agglomerated rows that were combined into Other/unclassified labels.
  sample_table <- aggregate(
    Abundance ~ Sample + Group + Taxon_display,
    composition,
    sum
  )
  group_table <- aggregate(
    Abundance ~ Group + Taxon_display,
    sample_table,
    mean
  )

  taxon_order <- c(top_taxa, "Other")
  sample_table$Taxon_display <- factor(sample_table$Taxon_display, levels = taxon_order)
  group_table$Taxon_display <- factor(group_table$Taxon_display, levels = taxon_order)
  sample_table$Sample <- factor(
    sample_table$Sample,
    levels = unique(sample_table$Sample[order(sample_table$Group)])
  )

  colours <- setNames(
    hcl.colors(length(taxon_order), palette = "Dynamic"),
    taxon_order
  )
  colours["Other"] <- "grey80"

  p_samples <- ggplot(
    sample_table,
    aes(x = Sample, y = Abundance, fill = Taxon_display)
  ) +
    geom_col(width = 1, colour = NA) +
    facet_grid(~Group, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = colours, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      x = "Individual samples",
      y = "Relative abundance (%)",
      fill = rank_name,
      title = paste(rank_name, "composition by disease status"),
      subtitle = paste0("Top ", top_n, " taxa; remaining taxa grouped as Other")
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )

  p_groups <- ggplot(
    group_table,
    aes(x = Group, y = Abundance, fill = Taxon_display)
  ) +
    geom_col(width = 0.70, colour = "white", linewidth = 0.15) +
    scale_fill_manual(values = colours, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      x = "Disease status",
      y = "Mean relative abundance (%)",
      fill = rank_name,
      title = paste("Mean", tolower(rank_name), "composition"),
      subtitle = paste0("Top ", top_n, " taxa; descriptive summary")
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      legend.position = "right"
    )

  prefix <- tolower(rank_name)
  write.csv(
    sample_table,
    file.path(output_dir, paste0(prefix, "_sample_relative_abundance.csv")),
    row.names = FALSE
  )
  write.csv(
    group_table,
    file.path(output_dir, paste0(prefix, "_group_mean_relative_abundance.csv")),
    row.names = FALSE
  )
  ggsave(
    file.path(output_dir, paste0(prefix, "_composition_by_sample.png")),
    p_samples, width = 11, height = 6, units = "in", dpi = 600, bg = "white"
  )
  ggsave(
    file.path(output_dir, paste0(prefix, "_composition_by_sample.pdf")),
    p_samples, width = 11, height = 6, units = "in"
  )
  ggsave(
    file.path(output_dir, paste0(prefix, "_mean_composition_by_group.png")),
    p_groups, width = 8, height = 6, units = "in", dpi = 600, bg = "white"
  )
  ggsave(
    file.path(output_dir, paste0(prefix, "_mean_composition_by_group.pdf")),
    p_groups, width = 8, height = 6, units = "in"
  )

  print(p_samples)
  print(p_groups)
  invisible(list(sample = sample_table, group = group_table))
}

composition_results <- lapply(taxonomic_ranks, function(rank_name) {
  make_composition_plots(
    ps_filtered,
    rank_name,
    unname(top_taxa_to_show[rank_name])
  )
})
names(composition_results) <- taxonomic_ranks

# ---- Top-20 genus grouped bar chart -----------------------------------------
# Rank genera by their mean relative abundance across all retained samples.
# Then display one adjacent bar per disease group for each genus.
ps_genus <- tax_glom(ps_filtered, taxrank = "Genus", NArm = FALSE)
ps_genus_relative <- transform_sample_counts(
  ps_genus,
  function(x) 100 * x / sum(x)
)
genus_long <- psmelt(ps_genus_relative)
genus_long$Genus_label <- as.character(genus_long$Genus)
genus_long$Genus_label[
  is.na(genus_long$Genus_label) | trimws(genus_long$Genus_label) == ""
] <- "Unclassified genus"
genus_long$Group <- factor(genus_long[[group_column]])

# Sum rows sharing the same displayed genus within each sample.
genus_by_sample <- aggregate(
  Abundance ~ Sample + Group + Genus_label,
  genus_long,
  sum
)

overall_genus_mean <- aggregate(
  Abundance ~ Genus_label,
  genus_by_sample,
  mean
)
overall_genus_mean <- overall_genus_mean[
  order(overall_genus_mean$Abundance, decreasing = TRUE),
  ,
  drop = FALSE
]
top_20_genera <- head(
  overall_genus_mean$Genus_label,
  top_genus_grouped_bars
)

top_genus_samples <- genus_by_sample[
  genus_by_sample$Genus_label %in% top_20_genera,
  ,
  drop = FALSE
]

# Mean, standard deviation, standard error, and sample count by group and genus.
genus_mean <- aggregate(
  Abundance ~ Group + Genus_label,
  top_genus_samples,
  mean
)
genus_sd <- aggregate(
  Abundance ~ Group + Genus_label,
  top_genus_samples,
  sd
)
genus_n <- aggregate(
  Abundance ~ Group + Genus_label,
  top_genus_samples,
  length
)
names(genus_mean)[names(genus_mean) == "Abundance"] <- "Mean_abundance"
names(genus_sd)[names(genus_sd) == "Abundance"] <- "SD"
names(genus_n)[names(genus_n) == "Abundance"] <- "n"

genus_group_summary <- merge(
  merge(genus_mean, genus_sd, by = c("Group", "Genus_label"), all = TRUE),
  genus_n,
  by = c("Group", "Genus_label"),
  all = TRUE
)
genus_group_summary$SE <- genus_group_summary$SD / sqrt(genus_group_summary$n)
genus_group_summary$Genus_label <- factor(
  genus_group_summary$Genus_label,
  levels = rev(top_20_genera)
)

# Median and IQR complement mean and SE because abundance distributions are
# commonly skewed and zero-inflated.
genus_median_iqr <- do.call(
  rbind,
  lapply(
    split(
      top_genus_samples$Abundance,
      interaction(
        top_genus_samples$Group,
        top_genus_samples$Genus_label,
        drop = TRUE,
        sep = "___"
      )
    ),
    function(x) {
      q <- quantile(x, c(0.25, 0.50, 0.75), names = FALSE, na.rm = TRUE)
      data.frame(median = q[2], Q1 = q[1], Q3 = q[3], IQR = q[3] - q[1])
    }
  )
)
median_keys <- strsplit(rownames(genus_median_iqr), "___", fixed = TRUE)
genus_median_iqr$Group <- vapply(median_keys, `[`, character(1), 1)
genus_median_iqr$Genus_label <- vapply(median_keys, `[`, character(1), 2)
rownames(genus_median_iqr) <- NULL
genus_median_iqr <- genus_median_iqr[
  ,
  c("Group", "Genus_label", "median", "Q1", "Q3", "IQR"),
  drop = FALSE
]

group_levels <- levels(genus_group_summary$Group)
if (length(group_levels) == 2) {
  group_colours <- setNames(c("#2166AC", "#B2182B"), group_levels)
} else {
  group_colours <- setNames(
    hcl.colors(length(group_levels), palette = "Dark 3"),
    group_levels
  )
}

p_top15_genus <- ggplot(
  genus_group_summary,
  aes(
    x = Genus_label,
    y = Mean_abundance,
    fill = Group
  )
) +
  geom_col(
    position = position_dodge(width = 0.66),
    width = 0.50,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(0, Mean_abundance - SE),
      ymax = Mean_abundance + SE
    ),
    position = position_dodge(width = 0.66),
    width = 0.22,
    linewidth = 0.45
  ) +
  coord_flip() +
  scale_fill_manual(values = group_colours) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = "Mean relative abundance (%)",
    fill = "Disease status",
    title = "Top 15 most abundant genera by disease status",
    subtitle = "Bars show group means; error bars show ±1 standard error"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(
      colour = "black",
      face = "italic",
      size = 10
    ),
    axis.text.x = element_text(colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.position = "top",
    plot.margin = margin(10, 20, 10, 10)
  )

write.csv(
  genus_group_summary,
  file.path(output_dir, "top15_genus_group_mean_abundance.csv"),
  row.names = FALSE
)
write.csv(
  top_genus_samples,
  file.path(output_dir, "top15_genus_sample_abundance.csv"),
  row.names = FALSE
)
write.csv(
  genus_median_iqr,
  file.path(output_dir, "top15_genus_group_median_IQR.csv"),
  row.names = FALSE
)
composition_summary_text <- capture.output({
  cat("TAXONOMIC-COMPOSITION DESCRIPTIVE SUMMARY\n")
  cat("Generated:", format(Sys.time()), "\n")
  cat("Abundance unit: percent relative abundance\n")
  cat("Prevalence threshold:", minimum_prevalence_samples, "samples\n")
  cat("Top genera displayed:", top_genus_grouped_bars, "\n\n")
  cat("Group mean, SD, SE, and sample size by genus\n")
  print(genus_group_summary, row.names = FALSE)
  cat("\nGroup median and IQR by genus\n")
  print(genus_median_iqr, row.names = FALSE)
  cat(
    "\nNote: These are descriptive statistics. ",
    "ANCOM-BC2 provides the inferential taxon-level comparisons.\n",
    sep = ""
  )
})
writeLines(
  composition_summary_text,
  file.path(output_dir, "taxonomic_composition_statistical_summary.txt")
)
ggsave(
  file.path(output_dir, "top15_genus_grouped_bar_chart.png"),
  p_top15_genus,
  width = 9,
  height = 9.5,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(output_dir, "top15_genus_grouped_bar_chart.pdf"),
  p_top15_genus,
  width = 9,
  height = 9.5,
  units = "in"
)
print(p_top15_genus)

message("Taxonomic-composition outputs saved in: ", normalizePath(output_dir))
