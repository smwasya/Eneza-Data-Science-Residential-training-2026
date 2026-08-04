#!/usr/bin/env Rscript

############################################################
# DADA2 Step 4: Merge paired-end reads
############################################################

library(dada2)

#-----------------------------------------------------------
# Paths
#-----------------------------------------------------------

results_dir <- "/home/sam/ENEZA/dada2/results_2"

#-----------------------------------------------------------
# Load objects
#-----------------------------------------------------------

dadaFs <- readRDS(file.path(results_dir, "dadaFs.rds"))
dadaRs <- readRDS(file.path(results_dir, "dadaRs.rds"))
sample.names <- readRDS(file.path(results_dir, "sample_names.rds"))

filt_path <- "/home/sam/ENEZA/results/fastp_results/disease_25_2"

filtFs <- sort(list.files(
    filt_path,
    pattern="_F_filt.fastq.gz$",
    full.names=TRUE
))

filtRs <- sort(list.files(
    filt_path,
    pattern="_R_filt.fastq.gz$",
    full.names=TRUE
))

#-----------------------------------------------------------
# Merge reads
#-----------------------------------------------------------

mergers <- vector("list", length(sample.names))

for(i in seq_along(sample.names)){

    cat("Merging", sample.names[i], "\n")

    derepF <- derepFastq(filtFs[i])
    derepR <- derepFastq(filtRs[i])

    names(derepF) <- sample.names[i]
    names(derepR) <- sample.names[i]

    mergers[[i]] <- mergePairs(
        dadaFs[[i]],
        derepF,
        dadaRs[[i]],
        derepR,
        verbose=TRUE
    )

    rm(derepF, derepR)
    gc()
}

names(mergers) <- sample.names

saveRDS(
    mergers,
    file.path(results_dir, "mergers.rds")
)

cat("Finished merging reads.\n")
