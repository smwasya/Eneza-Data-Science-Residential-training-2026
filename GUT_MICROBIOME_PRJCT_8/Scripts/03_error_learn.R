#!/usr/bin/env Rscript

library(dada2)

threads <- 4

path <- "/home/sam/ENEZA/results/fastp_results/trimmed_fq"
filt_path <- file.path(path, "filtered")

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

cat("Forward files:", length(filtFs), "\n")
cat("Reverse files:", length(filtRs), "\n")

dir.create("results", showWarnings = FALSE)

cat("Learning forward errors...\n")
errF <- learnErrors(filtFs, multithread = threads)

saveRDS(errF, "results/errF.rds")

cat("Forward errors saved.\n")

cat("Learning reverse errors...\n")
errR <- learnErrors(filtRs, multithread = threads)

saveRDS(errR, "results/errR.rds")

cat("Reverse errors saved.\n")

pdf("results/error_rates.pdf")

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

dev.off()

cat("Finished successfully!\n")
