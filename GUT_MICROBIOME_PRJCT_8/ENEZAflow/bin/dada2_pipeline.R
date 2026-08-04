#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(dada2)
})

###########################################################
## Arguments
###########################################################

option_list <- list(

make_option(
"--samplesheet",
type="character"
),

make_option(
"--silva",
type="character"
),

make_option(
"--output",
type="character"
),

make_option(
"--threads",
type="integer",
default=4
)

)

opt <- parse_args(
OptionParser(option_list=option_list)
)

dir.create(
opt$output,
showWarnings=FALSE,
recursive=TRUE
)

samples <- read.delim(

opt$samplesheet,

stringsAsFactors=FALSE

)

required <- c(
"SampleID",
"forward",
"reverse"
)

missing <- setdiff(
required,
colnames(samples)
)

if(length(missing)>0){

stop(
paste(
"Missing columns:",
paste(missing,collapse=", ")
)
)

}

filtFs <- samples$forward

filtRs <- samples$reverse

sample.names <- samples$SampleID

if(any(!file.exists(filtFs))){

stop("Forward FASTQ missing.")

}

if(any(!file.exists(filtRs))){

stop("Reverse FASTQ missing.")

}

filtered_dir <- file.path(
opt$output,
"filtered"
)

dir.create(
filtered_dir,
showWarnings=FALSE
)

filtFs.out <- file.path(

filtered_dir,

paste0(
sample.names,
"_F_filt.fastq.gz"
)

)

filtRs.out <- file.path(

filtered_dir,

paste0(
sample.names,
"_R_filt.fastq.gz"
)

)

cat("\nFiltering reads...\n")

out <- filterAndTrim(

fwd=filtFs,

filt=filtFs.out,

rev=filtRs,

filt.rev=filtRs.out,

truncLen=c(240,200),

maxN=0,

maxEE=c(2,2),

truncQ=2,

rm.phix=TRUE,

compress=TRUE,

multithread=opt$threads

)

write.csv(

out,

file.path(
opt$output,
"filtering_summary.csv"
)

)

cat("\nLearning forward error model...\n")

errF <- learnErrors(

filtFs.out,

multithread=opt$threads

)

cat("\nLearning reverse error model...\n")

errR <- learnErrors(

filtRs.out,

multithread=opt$threads

)

pdf(

file.path(
opt$output,
"error_rates.pdf"
)

)

plotErrors(
errF,
nominalQ=TRUE
)

plotErrors(
errR,
nominalQ=TRUE
)

dev.off()


dadaFs <- vector(
"list",
length(sample.names)
)

dadaRs <- vector(
"list",
length(sample.names)
)

mergers <- vector(
"list",
length(sample.names)
)

for(i in seq_along(sample.names)){

    cat("---------------------------------\n")
    cat("Sample", i, "/", length(sample.names), "\n")
    cat(sample.names[i], "\n")
    cat("---------------------------------\n")

    ####################################################
    ## Dereplication
    ####################################################

    derepF <- derepFastq(filtFs.out[i])
    derepR <- derepFastq(filtRs.out[i])

    names(derepF) <- sample.names[i]
    names(derepR) <- sample.names[i]

    ####################################################
    ## Denoising
    ####################################################

    dadaFs[[i]] <- dada(
        derepF,
        err = errF,
        multithread = FALSE
    )

    dadaRs[[i]] <- dada(
        derepR,
        err = errR,
        multithread = FALSE
    )

    ####################################################
    ## Merge reads
    ####################################################

    mergers[[i]] <- mergePairs(
        dadaFs[[i]],
        derepF,
        dadaRs[[i]],
        derepR
    )

    rm(derepF, derepR)

    gc()
}

names(dadaFs) <- sample.names
names(dadaRs) <- sample.names
names(mergers) <- sample.names

saveRDS(
    dadaFs,
    file.path(opt$output, "dadaFs.rds")
)

saveRDS(
    dadaRs,
    file.path(opt$output, "dadaRs.rds")
)

saveRDS(
    mergers,
    file.path(opt$output, "mergers.rds")
)

seqtab <- makeSequenceTable(
mergers
)

saveRDS(

seqtab,

file.path(
opt$output,
"seqtab.rds"
)

)

seqtab.nochim <- removeBimeraDenovo(

seqtab,

method="consensus",

multithread=opt$threads,

verbose=TRUE

)

saveRDS(

seqtab.nochim,

file.path(
opt$output,
"seqtab_nochim.rds"
)

)

taxa <- assignTaxonomy(

seqtab.nochim,

opt$silva,

multithread=opt$threads

)

saveRDS(

taxa,

file.path(
opt$output,
"taxonomy.rds"
)

)

write.csv(

taxa,

file.path(
opt$output,
"taxonomy.csv"

)

)


cat("\n=====================================\n")

cat("Samples processed :", length(sample.names), "\n")

cat("ASVs before chimera :", ncol(seqtab), "\n")

cat("ASVs after chimera  :", ncol(seqtab.nochim), "\n")

cat("Results written to :", opt$output, "\n")

cat("=====================================\n")
