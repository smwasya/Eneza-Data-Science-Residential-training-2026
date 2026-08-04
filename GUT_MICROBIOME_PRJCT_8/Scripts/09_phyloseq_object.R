#!/usr/bin/env Rscript

############################################################
# Create phyloseq object
############################################################

library(phyloseq)

#-----------------------------------------------------------
# Paths
#-----------------------------------------------------------

results_dir <- "/home/sam/ENEZA/dada2/results_2"

metadata_file <- "/home/sam/ENEZA/dada2/results_2/metadata_all.tsv"

#-----------------------------------------------------------
# Load data
#-----------------------------------------------------------

seqtab <- readRDS(
    file.path(results_dir, "131_seqtab_nochim_merged.rds")
)

taxa <- readRDS(
    file.path(results_dir, "taxonomy_merged.rds")
)

metadata <- read.table(
  metadata_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

#-----------------------------------------------------------
# Match sample order
#-----------------------------------------------------------

metadata <- metadata[match(rownames(seqtab), rownames(metadata)), , drop = FALSE]

all(rownames(seqtab) == rownames(metadata))

#######################

metadata <- metadata[rownames(seqtab), ]

stopifnot(all(rownames(seqtab) == rownames(metadata)))

#-----------------------------------------------------------
# Build phyloseq object
#-----------------------------------------------------------
library(phyloseq)

OTU <- otu_table(seqtab, taxa_are_rows = FALSE)

TAX <- tax_table(as.matrix(taxa))

SAM <- sample_data(metadata)

ps <- phyloseq(
  OTU,
  TAX,
  SAM
)

ps

#------------------------------------
nsamples(ps)
ntaxa(ps)
sample_variables(ps)
rank_names(ps)
#-----------------------------------------------------------
# Save
#-----------------------------------------------------------

saveRDS(
    ps,
    file.path(results_dir, "phyloseq.rds")
)

cat("=====================================\n")
cat("Phyloseq object created successfully\n")
cat("Samples :", nsamples(ps), "\n")
cat("ASVs    :", ntaxa(ps), "\n")
cat("=====================================\n")



################################



seqtab <- readRDS("results/seqtab_nochim_merged.rds")
taxa <- readRDS("results/taxonomy.rds")

# Original sequences
seqs <- colnames(seqtab)

# Create ASV IDs
asv_ids <- paste0("ASV", seq_along(seqs))

# Rename abundance table columns
colnames(seqtab) <- asv_ids

# Rename taxonomy rows
rownames(taxa) <- asv_ids

###############

write.csv(
  seqtab,
  "tables/asv_abundance_table.csv"
)
####################3

write.csv(
  taxa,
  "tables/taxonomy_table.csv"
)
#######################################


con <- file("tables/rep_seqs.fasta", "w")

for (i in seq_along(seqs)) {
  writeLines(paste0(">", asv_ids[i]), con)
  writeLines(seqs[i], con)
}

close(con)


