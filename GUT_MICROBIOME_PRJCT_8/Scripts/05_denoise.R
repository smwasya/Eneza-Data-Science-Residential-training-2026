#!/usr/bin/env Rscript

############################################################
# DADA2 Step 3: Dereplication + Denoising
# Test run on first 20 samples
############################################################

library(dada2)

threads <- 4

#-----------------------------------------------------------
# Paths
#-----------------------------------------------------------

filt_path <- "/home/sam/ENEZA/results/fastp_results/disease_25_2"
results_dir <- "/home/sam/ENEZA/dada2/results_2"

#-----------------------------------------------------------
# Read filtered FASTQs
#-----------------------------------------------------------

filtFs <- sort(list.files(
  filt_path,
  pattern = "_F_filt.fastq.gz$",
  full.names = TRUE
))

filtRs <- sort(list.files(
  filt_path,
  pattern = "_R_filt.fastq.gz$",
  full.names = TRUE
))

sample.names <- sub("_F_filt.fastq.gz$", "", basename(filtFs))


#-----------------------------------------------------------
# Process all samples
#-----------------------------------------------------------

stopifnot(length(filtFs) == length(filtRs))

n.samples <- length(filtFs)

cat("Processing", n.samples, "samples\n")
#-----------------------------------------------------------
# Load error models
#-----------------------------------------------------------

errF <- readRDS(file.path(results_dir, "errF.rds"))
errR <- readRDS(file.path(results_dir, "errR.rds"))

#-----------------------------------------------------------
# Process one sample at a time
#-----------------------------------------------------------

dadaFs <- vector("list", n.samples)
dadaRs <- vector("list", n.samples)

for(i in seq_len(n.samples)){

    cat("\n=====================================\n")
    cat("Sample", i, "of", n.samples, "\n")
    cat(sample.names[i], "\n")
    cat("=====================================\n")

    # Dereplicate
    derepF <- derepFastq(filtFs[i])
    derepR <- derepFastq(filtRs[i])

    names(derepF) <- sample.names[i]
    names(derepR) <- sample.names[i]

    # Denoise
    dadaFs[[i]] <- dada(
        derepF,
        err = errF,
        multithread = threads
    )

    dadaRs[[i]] <- dada(
        derepR,
        err = errR,
        multithread = threads
    )

    # Free memory
    rm(derepF, derepR)
    gc()
}

names(dadaFs) <- sample.names
names(dadaRs) <- sample.names

#-----------------------------------------------------------
# Save results
#-----------------------------------------------------------

saveRDS(
  dadaFs,
  file.path(results_dir, "dadaFs.rds")
)

saveRDS(
  dadaRs,
  file.path(results_dir, "dadaRs.rds")
)

saveRDS(
  sample.names,
  file.path(results_dir, "sample_names.rds")
)

cat("\n=====================================\n")
cat("Finished processing first", n.samples, "samples.\n")
cat("Results saved in:", results_dir, "\n")
cat("=====================================\n")
