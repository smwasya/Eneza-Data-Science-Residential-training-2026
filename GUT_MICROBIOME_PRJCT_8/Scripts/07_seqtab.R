#!/usr/bin/env Rscript

############################################################
# DADA2 Step 5: Create ASV table and remove chimeras
############################################################

library(dada2)

threads <- 4

#-----------------------------------------------------------
# Paths
#-----------------------------------------------------------

results_dir <- "/home/sam/ENEZA/dada2/results_2"

#-----------------------------------------------------------
# Load merged reads
#-----------------------------------------------------------

mergers <- readRDS(
    file.path(results_dir, "mergers.rds")
)

#-----------------------------------------------------------
# Construct sequence table
#-----------------------------------------------------------

cat("Creating sequence table...\n")

seqtab <- makeSequenceTable(mergers)

cat("Sequence table dimensions:\n")
print(dim(seqtab))

saveRDS(
    seqtab,
    file.path(results_dir, "seqtab.rds")
)

#-----------------------------------------------------------
# Remove chimeras
#-----------------------------------------------------------

cat("Removing chimeras...\n")

seqtab.nochim <- removeBimeraDenovo(
    seqtab,
    method = "consensus",
    multithread = threads,
    verbose = TRUE
)

cat("Non-chimeric sequence table dimensions:\n")
print(dim(seqtab.nochim))

saveRDS(
    seqtab.nochim,
    file.path(results_dir, "seqtab_nochim.rds")
)

#-----------------------------------------------------------
# Summary
#-----------------------------------------------------------

cat("\n=====================================\n")
cat("Total ASVs before chimera removal :", ncol(seqtab), "\n")
cat("Total ASVs after chimera removal  :", ncol(seqtab.nochim), "\n")
cat("=====================================\n")
