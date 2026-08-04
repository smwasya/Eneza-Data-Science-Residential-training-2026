# Differential abundance analysis with ANCOM-BC2
#
# This script compares microbial abundance between DISEASE groups using raw
# counts, bias correction, prevalence/library-size filtering, BH-adjusted
# p-values, and structural-zero detection.
#
# Required packages: phyloseq, ANCOMBC, ggplot2

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(ggplot2)
})

# ---- Configuration ----------------------------------------------------------
phyloseq_file <- NA_character_  # Leave as NA to use loaded `ps` or a file picker
metadata_file <- NA_character_  # Example: "sample_metadata.csv"
sample_id_column <- "SampleID"  # Used only when metadata_file is supplied
group_column <- "DISEASE"
taxonomic_rank <- "Genus"       # Change to "Species", "Family", etc. if needed
reference_group <- "Healthy"
prevalence_cutoff <- 0.10
library_size_cutoff <- 1000
fdr_cutoff <- 0.05
output_dir <- "differential_abundance_results"
set.seed(128)

# ---- Load and validate -------------------------------------------------------
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

if (!inherits(ps, "phyloseq")) stop("The RDS file is not a phyloseq object.")

# Optionally attach external metadata by exact sample-ID matching.
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
  ps_metadata <- data.frame(sample_data(ps), check.names = FALSE)
  ps_metadata[[group_column]] <- metadata_external[[group_column]][positions]
  rownames(ps_metadata) <- sample_names(ps)
  sample_data(ps) <- sample_data(ps_metadata)
}

if (!group_column %in% sample_variables(ps)) {
  stop(
    "'", group_column, "' is absent from sample_data(ps). ",
    "Add it to the phyloseq object or configure metadata_file."
  )
}
if (!taxonomic_rank %in% rank_names(ps)) {
  stop(
    "Taxonomic rank '", taxonomic_rank, "' is unavailable. Available ranks: ",
    paste(rank_names(ps), collapse = ", ")
  )
}
rank_output_dir <- file.path(output_dir, tolower(taxonomic_rank))

# ANCOM-BC2 requires untransformed count data.
count_matrix <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) count_matrix <- t(count_matrix)
if (any(count_matrix < 0, na.rm = TRUE)) stop("Counts cannot be negative.")
if (any(abs(count_matrix - round(count_matrix)) > .Machine$double.eps^0.5, na.rm = TRUE)) {
  stop("The phyloseq object contains non-integer values. Use the raw-count object.")
}

group_values <- as.character(data.frame(sample_data(ps))[[group_column]])
keep_samples <- !is.na(group_values) &
  trimws(group_values) != "" &
  sample_sums(ps) >= library_size_cutoff
if (sum(!keep_samples)) {
  message(
    "Removing ", sum(!keep_samples),
    " sample(s) with missing group or library size below ", library_size_cutoff, "."
  )
}

ps_da <- prune_samples(keep_samples, ps)
ps_da <- prune_taxa(taxa_sums(ps_da) > 0, ps_da)

metadata <- data.frame(sample_data(ps_da), check.names = FALSE)
metadata[[group_column]] <- factor(metadata[[group_column]])
if (!reference_group %in% levels(metadata[[group_column]])) {
  stop(
    "Reference group '", reference_group, "' was not found. Available groups: ",
    paste(levels(metadata[[group_column]]), collapse = ", ")
  )
}
metadata[[group_column]] <- relevel(metadata[[group_column]], ref = reference_group)
sample_data(ps_da) <- sample_data(metadata)
if (nlevels(metadata[[group_column]]) < 2) stop("At least two groups are required.")

# Remove taxa without an annotation at the requested rank, then agglomerate.
rank_values <- as.character(tax_table(ps_da)[, taxonomic_rank])
keep_taxa <- !is.na(rank_values) & trimws(rank_values) != ""
if (!any(keep_taxa)) stop("No taxa have a valid ", taxonomic_rank, " annotation.")
ps_da <- prune_taxa(keep_taxa, ps_da)
ps_aggregated <- tax_glom(ps_da, taxrank = taxonomic_rank, NArm = TRUE)

# Use readable, unique taxon identifiers in result tables.
taxon_labels <- as.character(tax_table(ps_aggregated)[, taxonomic_rank])
taxon_labels <- make.unique(taxon_labels, sep = "_")
taxa_names(ps_aggregated) <- taxon_labels

message("Samples per group:")
print(table(data.frame(sample_data(ps_aggregated))[[group_column]]))
message("Taxa entering ANCOM-BC2 before its prevalence filter: ", ntaxa(ps_aggregated))

# ---- ANCOM-BC2 ---------------------------------------------------------------
ancombc_fit <- ancombc2(
  data = ps_aggregated,
  fix_formula = group_column,
  rand_formula = NULL,
  p_adj_method = "BH",
  pseudo_sens = TRUE,
  prv_cut = prevalence_cutoff,
  lib_cut = 0,
  s0_perc = 0.05,
  group = group_column,
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = fdr_cutoff,
  n_cl = 1,
  verbose = TRUE,
  global = TRUE,
  pairwise = TRUE,
  dunnet = FALSE,
  trend = FALSE
)

result_wide <- ancombc_fit$res

# ANCOM-BC2 normally stores biological labels in a column named `taxon`.
# Its data-frame row names may only be 1, 2, 3, ... and must not be used as
# phylum/genus labels when the `taxon` column is available.
if ("taxon" %in% names(result_wide)) {
  result_wide$Taxon <- as.character(result_wide$taxon)
} else if (
  !is.null(rownames(result_wide)) &&
    all(rownames(result_wide) %in% taxa_names(ps_aggregated))
) {
  result_wide$Taxon <- rownames(result_wide)
} else if (nrow(result_wide) == ntaxa(ps_aggregated)) {
  result_wide$Taxon <- taxa_names(ps_aggregated)
} else {
  stop(
    "Taxon names could not be mapped from the ANCOM-BC2 results. ",
    "Inspect names(ancombc_fit$res) and head(ancombc_fit$res)."
  )
}

# Convert coefficient-specific columns to a convenient long table.
lfc_columns <- grep("^lfc_", names(result_wide), value = TRUE)
lfc_columns <- setdiff(lfc_columns, "lfc_(Intercept)")
if (!length(lfc_columns)) {
  stop("ANCOM-BC2 returned no non-intercept group coefficient.")
}

extract_contrast <- function(lfc_column) {
  contrast <- sub("^lfc_", "", lfc_column)
  get_column <- function(prefix, default = NA) {
    column <- paste0(prefix, contrast)
    if (column %in% names(result_wide)) result_wide[[column]] else rep(default, nrow(result_wide))
  }
  data.frame(
    Taxon = result_wide$Taxon,
    Contrast = contrast,
    log_fold_change = result_wide[[lfc_column]],
    standard_error = get_column("se_"),
    W = get_column("W_"),
    p_value = get_column("p_"),
    q_value = get_column("q_"),
    significant = as.logical(get_column("diff_", FALSE)),
    passed_sensitivity = as.logical(get_column("passed_ss_", FALSE)),
    stringsAsFactors = FALSE
  )
}

result_long <- do.call(rbind, lapply(lfc_columns, extract_contrast))
result_long$significant[is.na(result_long$significant)] <- FALSE
result_long$significant <- result_long$significant &
  !is.na(result_long$q_value) &
  result_long$q_value < fdr_cutoff
result_long$direction <- ifelse(
  result_long$log_fold_change > 0,
  "Higher in comparison group",
  "Higher in reference group"
)

significant_results <- result_long[result_long$significant, , drop = FALSE]
significant_results <- significant_results[
  order(significant_results$q_value, -abs(significant_results$log_fold_change)),
  ,
  drop = FALSE
]

# ---- Plots -------------------------------------------------------------------
plot_df <- result_long[
  is.finite(result_long$log_fold_change) & !is.na(result_long$q_value),
  ,
  drop = FALSE
]
plot_df$minus_log10_q <- -log10(pmax(plot_df$q_value, .Machine$double.xmin))

p_volcano <- ggplot(
  plot_df,
  aes(x = log_fold_change, y = minus_log10_q, colour = significant)
) +
  geom_hline(
    yintercept = -log10(fdr_cutoff),
    linetype = 2, linewidth = 0.6, colour = "grey40"
  ) +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.5, colour = "grey50") +
  geom_point(size = 2.4, alpha = 0.80) +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#B2182B")) +
  facet_wrap(~Contrast, scales = "free") +
  labs(
    x = "Bias-corrected log fold change",
    y = expression(-log[10]("BH-adjusted p-value")),
    colour = paste0("FDR < ", fdr_cutoff),
    title = paste("Differential abundance at", tolower(taxonomic_rank), "level"),
    subtitle = paste("ANCOM-BC2; reference group:", reference_group)
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    legend.title = element_text(face = "bold")
  )

# Label significant taxa in the volcano plot, as commonly done for genes in
# RNA-seq volcano plots. ggrepel is preferred because it reduces label overlap.
volcano_labels <- plot_df[plot_df$significant, , drop = FALSE]
volcano_labels <- volcano_labels[
  order(volcano_labels$q_value, -abs(volcano_labels$log_fold_change)),
  ,
  drop = FALSE
]

# Label five taxa while ensuring that both sides of the volcano are represented
# whenever significant positive and negative effects are available. Start with
# the best-supported taxon from each direction, then fill the remaining places
# by FDR rank.
volcano_labels$label_key <- paste(
  volcano_labels$Taxon,
  volcano_labels$Contrast,
  sep = "___"
)
best_negative <- head(
  volcano_labels[volcano_labels$log_fold_change < 0, , drop = FALSE],
  1
)
best_positive <- head(
  volcano_labels[volcano_labels$log_fold_change > 0, , drop = FALSE],
  1
)
required_labels <- rbind(best_negative, best_positive)
remaining_labels <- volcano_labels[
  !volcano_labels$label_key %in% required_labels$label_key,
  ,
  drop = FALSE
]
volcano_labels <- head(
  rbind(required_labels, remaining_labels),
  5
)
if (nrow(volcano_labels)) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_volcano <- p_volcano +
      ggrepel::geom_label_repel(
        data = volcano_labels,
        aes(label = Taxon),
        colour = "black",
        fill = "white",
        size = 3.1,
        fontface = "italic",
        box.padding = 0.75,
        point.padding = 0.45,
        label.padding = 0.20,
        label.r = grid::unit(0.12, "lines"),
        force = 8,
        force_pull = 0.5,
        nudge_x = ifelse(
          volcano_labels$log_fold_change > 0,
          0.25,
          -0.25
        ),
        min.segment.length = 0,
        segment.colour = "grey35",
        segment.linewidth = 0.4,
        max.overlaps = Inf,
        max.iter = 20000,
        max.time = 5,
        seed = 128,
        show.legend = FALSE
      )
  } else {
    warning(
      "Install ggrepel for non-overlapping volcano labels: ",
      "install.packages('ggrepel'). Using basic labels for now."
    )
    p_volcano <- p_volcano +
      geom_text(
        data = volcano_labels,
        aes(label = Taxon),
        colour = "black",
        size = 3,
        check_overlap = TRUE,
        nudge_x = ifelse(
          volcano_labels$log_fold_change > 0,
          0.20,
          -0.20
        ),
        show.legend = FALSE
      )
  }
}

# Extra vertical and horizontal room prevents labels from being clipped at the
# plot boundaries.
p_volcano <- p_volcano +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.22))) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(15, 35, 15, 20))

p_effects <- NULL
ranked_results <- result_long[
  !is.na(result_long$q_value) &
    is.finite(result_long$log_fold_change) &
    is.finite(result_long$standard_error),
  ,
  drop = FALSE
]
ranked_results <- ranked_results[
  order(!ranked_results$significant, ranked_results$q_value),
  ,
  drop = FALSE
]
if (nrow(ranked_results)) {
  # Display significant taxa first. If fewer than 20 are significant, include
  # the best-supported non-significant taxa so the taxon-level axis is retained.
  top_results <- head(ranked_results, 20)
  top_results$label <- factor(
    top_results$Taxon,
    levels = rev(unique(top_results$Taxon))
  )
  top_results$lower <- top_results$log_fold_change - 1.96 * top_results$standard_error
  top_results$upper <- top_results$log_fold_change + 1.96 * top_results$standard_error

  p_effects <- ggplot(
    top_results,
    aes(x = log_fold_change, y = label, colour = direction, shape = significant)
  ) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbarh(
      aes(xmin = lower, xmax = upper),
      height = 0.20, linewidth = 0.7
    ) +
    geom_point(size = 2.8) +
    facet_wrap(~Contrast, scales = "free_y") +
    scale_colour_manual(
      values = c(
        "Higher in comparison group" = "#B2182B",
        "Higher in reference group" = "#2166AC"
      )
    ) +
    scale_shape_manual(
      values = c(`FALSE` = 1, `TRUE` = 16),
      labels = c(`FALSE` = "No", `TRUE` = "Yes")
    ) +
    labs(
      x = "Bias-corrected log fold change (95% CI)",
      y = taxonomic_rank,
      colour = "Direction",
      shape = paste0("FDR < ", fdr_cutoff),
      title = paste("Differential abundance at", tolower(taxonomic_rank), "level"),
      subtitle = paste0(
        "ANCOM-BC2; taxa are shown on the y-axis; reference: ", reference_group
      )
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      legend.position = "bottom"
    )
}

# ---- Save outputs ------------------------------------------------------------
dir.create(rank_output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(result_wide, file.path(rank_output_dir, "ancombc2_results_wide.csv"), row.names = FALSE)
write.csv(result_long, file.path(rank_output_dir, "ancombc2_results_long.csv"), row.names = FALSE)
write.csv(
  significant_results,
  file.path(rank_output_dir, "ancombc2_significant_taxa.csv"),
  row.names = FALSE
)

if (!is.null(ancombc_fit$res_global)) {
  global_results <- ancombc_fit$res_global
  global_results$Taxon <- rownames(global_results)
  write.csv(
    global_results,
    file.path(rank_output_dir, "ancombc2_global_test.csv"),
    row.names = FALSE
  )
}
if (!is.null(ancombc_fit$zero_ind)) {
  write.csv(
    ancombc_fit$zero_ind,
    file.path(rank_output_dir, "ancombc2_structural_zeros.csv"),
    row.names = TRUE
  )
}

da_summary_text <- capture.output({
  cat("DIFFERENTIAL-ABUNDANCE STATISTICAL SUMMARY\n")
  cat("Generated:", format(Sys.time()), "\n")
  cat("Method: ANCOM-BC2\n")
  cat("Taxonomic rank:", taxonomic_rank, "\n")
  cat("Reference group:", reference_group, "\n")
  cat("Prevalence cutoff:", prevalence_cutoff, "\n")
  cat("FDR cutoff:", fdr_cutoff, "\n\n")
  cat("Samples per group\n")
  print(table(data.frame(sample_data(ps_aggregated))[[group_column]]))
  cat("\nNumber of taxa tested:", length(unique(result_long$Taxon)), "\n")
  cat("Significant taxon-contrast results:", nrow(significant_results), "\n\n")
  if (nrow(significant_results)) {
    cat("Significant results ranked by BH-adjusted p-value\n")
    print(
      significant_results[
        ,
        c(
          "Taxon", "Contrast", "log_fold_change", "standard_error",
          "W", "p_value", "q_value", "direction"
        ),
        drop = FALSE
      ],
      row.names = FALSE
    )
  } else {
    cat("No taxa met the specified FDR threshold.\n")
  }
})
writeLines(
  da_summary_text,
  file.path(rank_output_dir, "ancombc2_statistical_summary.txt")
)

ggsave(
  file.path(rank_output_dir, "ancombc2_volcano.png"),
  p_volcano, width = 7, height = 5.5, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(rank_output_dir, "ancombc2_volcano.pdf"),
  p_volcano, width = 7, height = 5.5, units = "in"
)

if (!is.null(p_effects)) {
  ggsave(
    file.path(rank_output_dir, paste0("ancombc2_", tolower(taxonomic_rank), "_taxa_effects.png")),
    p_effects, width = 8, height = 6, units = "in", dpi = 600, bg = "white"
  )
  ggsave(
    file.path(rank_output_dir, paste0("ancombc2_", tolower(taxonomic_rank), "_taxa_effects.pdf")),
    p_effects, width = 8, height = 6, units = "in"
  )
}

print(p_volcano)
if (!is.null(p_effects)) print(p_effects)
message("Significant taxon-contrast results: ", nrow(significant_results))
message("Differential-abundance outputs saved in: ", normalizePath(rank_output_dir))
