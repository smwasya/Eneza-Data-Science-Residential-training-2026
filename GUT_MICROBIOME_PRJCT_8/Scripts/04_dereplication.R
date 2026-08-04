library(dada2)

threads <- 1        # one sample at a time

filt_path <- "/home/sam/ENEZA/results/fastp_results/trimmed_fq/filtered"
results <- "/home/sam/ENEZA/dada2/results"

filtFs <- sort(list.files(
    filt_path,
    pattern="_F_filt.fastq.gz$",
    full.names=TRUE))

filtRs <- sort(list.files(
    filt_path,
    pattern="_R_filt.fastq.gz$",
    full.names=TRUE))

sample.names <- sub("_F_filt.fastq.gz","",basename(filtFs))

errF <- readRDS(file.path(results,"errF.rds"))
errR <- readRDS(file.path(results,"errR.rds"))

dadaFs <- vector("list", length(filtFs))
dadaRs <- vector("list", length(filtRs))

for(i in seq_along(filtFs)){

    cat("\n")
    cat("=====================================\n")
    cat(i,"/",length(filtFs),sample.names[i],"\n")
    cat("=====================================\n")

    derepF <- derepFastq(filtFs[i])
    derepR <- derepFastq(filtRs[i])

    dadaFs[[i]] <- dada(
        derepF,
        err=errF,
        multithread=FALSE
    )

    dadaRs[[i]] <- dada(
        derepR,
        err=errR,
        multithread=FALSE
    )

    rm(derepF,derepR)

    gc()
}

names(dadaFs) <- sample.names
names(dadaRs) <- sample.names

saveRDS(dadaFs,file.path(results,"dadaFs.rds"))
saveRDS(dadaRs,file.path(results,"dadaRs.rds"))
