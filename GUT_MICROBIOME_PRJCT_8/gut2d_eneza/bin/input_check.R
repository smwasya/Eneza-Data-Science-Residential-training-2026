#!/usr/bin/env Rscript
#install.packages("optparse")
library(optparse)

############################################################
# Arguments
############################################################

option_list <- list(
  
  make_option("--metadata", type = "character"),
  
  make_option("--reads", type = "character"),
  
  make_option("--output", type = "character")
  
)

opt <- parse_args(
  OptionParser(option_list = option_list)
)

############################################################
# Check FASTQ directory
############################################################

if (!dir.exists(opt$reads)) {
  stop("FASTQ directory does not exist.")
}

############################################################
# Read metadata
############################################################

meta <- read.delim(
  opt$metadata,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

############################################################
# Check required columns
############################################################

required <- c("SampleID", "Disease")

missing <- setdiff(required, colnames(meta))

if (length(missing) > 0) {
  
  stop(
    paste(
      "Missing metadata columns:",
      paste(missing, collapse = ", ")
    )
  )
  
}

############################################################
# Check duplicates
############################################################

dup <- duplicated(meta$SampleID)

if (any(dup)) {
  
  stop("Duplicate SampleID values detected.")
  
}

############################################################
# Find FASTQ files
############################################################

fastq_files <- list.files(
  path = opt$reads,
  pattern = "\\.(fastq|fq)\\.gz$",
  full.names = TRUE
)

if (length(fastq_files) == 0) {
  stop("No FASTQ files found.")
}

############################################################
# Extract filenames and sample names
############################################################

files <- basename(fastq_files)

samples <- unique(
  sub("_[12]\\.(fastq|fq)\\.gz$", "", files)
)

############################################################
# Check paired reads
############################################################

for (s in samples) {
  
  fwd_exists <- any(
    basename(fastq_files) %in%
      c(paste0(s, "_1.fastq.gz"),
        paste0(s, "_1.fq.gz"))
  )
  
  rev_exists <- any(
    basename(fastq_files) %in%
      c(paste0(s, "_2.fastq.gz"),
        paste0(s, "_2.fq.gz"))
  )
  
  if (!fwd_exists)
    stop(paste("Missing forward read for", s))
  
  if (!rev_exists)
    stop(paste("Missing reverse read for", s))
  
}

############################################################
# Check metadata consistency
############################################################

missing_meta <- setdiff(
  samples,
  meta$SampleID
)

if (length(missing_meta) > 0) {
  
  stop(
    paste(
      "FASTQ samples missing from metadata:",
      paste(missing_meta, collapse = ", ")
    )
  )
  
}

############################################################
# Check metadata without FASTQs
############################################################

missing_fastq <- setdiff(
  meta$SampleID,
  samples
)

if (length(missing_fastq) > 0) {
  
  stop(
    paste(
      "Metadata samples missing FASTQs:",
      paste(missing_fastq, collapse = ", ")
    )
  )
  
}

############################################################
# Build validated sample sheet
############################################################

samplesheet <- data.frame()

for (s in samples) {
  
  fwd <- fastq_files[
    grepl(
      paste0("^", s, "_1\\.(fastq|fq)\\.gz$"),
      basename(fastq_files)
    )
  ]
  
  rev <- fastq_files[
    grepl(
      paste0("^", s, "_2\\.(fastq|fq)\\.gz$"),
      basename(fastq_files)
    )
  ]
  
  row <- meta[
    meta$SampleID == s,
    ,
    drop = FALSE
  ]
  
  row$forward <- normalizePath(fwd)
  row$reverse <- normalizePath(rev)
  
  samplesheet <- rbind(
    samplesheet,
    row
  )
  
}

############################################################
# Save validated sample sheet
############################################################

write.table(
  samplesheet,
  file = opt$output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nInput validation completed successfully.\n")
cat("Samples validated:", nrow(samplesheet), "\n")
cat("Output written to:", opt$output, "\n")