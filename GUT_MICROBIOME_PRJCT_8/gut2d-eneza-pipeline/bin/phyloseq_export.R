#!/usr/bin/env Rscript

############################################################
# Create phyloseq object and machine learning dataset
############################################################

suppressPackageStartupMessages({
    library(optparse)
    library(phyloseq)
})

############################################################
# Arguments
############################################################

option_list <- list(

    make_option(
        "--seqtab",
        type = "character"
    ),

    make_option(
        "--taxonomy",
        type = "character"
    ),

    make_option(
        "--metadata",
        type = "character"
    ),

    make_option(
        "--output",
        type = "character"
    )

)

opt <- parse_args(
    OptionParser(option_list = option_list)
)

############################################################
# Checks
############################################################

if (is.null(opt$seqtab))
    stop("Missing --seqtab")

if (is.null(opt$taxonomy))
    stop("Missing --taxonomy")

if (is.null(opt$metadata))
    stop("Missing --metadata")

if (is.null(opt$output))
    stop("Missing --output")

dir.create(
    opt$output,
    recursive = TRUE,
    showWarnings = FALSE
)

############################################################
# Load data
############################################################

cat("Loading sequence table...\n")

seqtab <- readRDS(opt$seqtab)

cat("Loading taxonomy...\n")

taxa <- readRDS(opt$taxonomy)

cat("Loading metadata...\n")

metadata <- read.delim(
    opt$metadata,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

############################################################
# Match sample order
############################################################

common.samples <- intersect(
    rownames(seqtab),
    rownames(metadata)
)

seqtab <- seqtab[common.samples, , drop = FALSE]

metadata <- metadata[common.samples, , drop = FALSE]

stopifnot(
    all(
        rownames(seqtab) == rownames(metadata)
    )
)

############################################################
# Build phyloseq object
############################################################

OTU <- otu_table(
    seqtab,
    taxa_are_rows = FALSE
)

TAX <- tax_table(
    as.matrix(taxa)
)

SAM <- sample_data(
    metadata
)

ps <- phyloseq(
    OTU,
    TAX,
    SAM
)

############################################################
# Keep bacteria only
############################################################

ps <- subset_taxa(
    ps,
    Kingdom == "Bacteria"
)

############################################################
# Remove mitochondria
############################################################

ps <- subset_taxa(
    ps,
    Family != "Mitochondria" |
        is.na(Family)
)

############################################################
# Remove chloroplast
############################################################

ps <- subset_taxa(
    ps,
    Order != "Chloroplast" |
        is.na(Order)
)

############################################################
# Relative abundance
############################################################

ps.rel <- transform_sample_counts(

    ps,

    function(x){

        x / sum(x)

    }

)

############################################################
# Aggregate to genus
############################################################

ps.genus <- tax_glom(

    ps.rel,

    taxrank = "Genus",

    NArm = TRUE

)

############################################################
# Remove unknown genera
############################################################

genus <- as.character(
    tax_table(ps.genus)[,"Genus"]
)

keep <-

    !is.na(genus) &

    genus != "" &

    !grepl(

        "uncultured|unclassified|unknown|metagenome",

        genus,

        ignore.case = TRUE

    )

ps.genus <- prune_taxa(

    keep,

    ps.genus

)

############################################################
# Export abundance table
############################################################

ml <- as.data.frame(

    otu_table(ps.genus)

)

if(taxa_are_rows(ps.genus)){

    ml <- t(ml)

    ml <- as.data.frame(ml)

}

############################################################
# Rename columns
############################################################

genus <-

    tax_table(ps.genus)[,"Genus"]

colnames(ml) <-

    make.unique(

        as.character(genus)

    )

############################################################
# Add disease labels
############################################################

ml$Disease <-

    sample_data(ps.genus)$DISEASE

############################################################
# Export taxonomy table
############################################################

tax.export <-

    as.data.frame(

        tax_table(ps.genus)

    )

write.csv(

    tax.export,

    file.path(

        opt$output,

        "taxonomy_table.csv"

    )

)

############################################################
# Export abundance table
############################################################

otu.export <-

    as.data.frame(

        otu_table(ps)

    )

if(taxa_are_rows(ps)){

    otu.export <- t(otu.export)

    otu.export <- as.data.frame(otu.export)

}

write.csv(

    otu.export,

    file.path(

        opt$output,

        "asv_abundance_table.csv"

    )

)

############################################################
# Machine learning dataset
############################################################

write.csv(

    ml,

    file.path(

        opt$output,

        "machine_learning_dataset.csv"

    ),

    row.names = TRUE

)

############################################################
# Save phyloseq objects
############################################################

saveRDS(

    ps,

    file.path(

        opt$output,

        "phyloseq.rds"

    )

)

saveRDS(

    ps.genus,

    file.path(

        opt$output,

        "phyloseq_genus.rds"

    )

)

############################################################
# Representative sequences
############################################################

seqs <- colnames(seqtab)

asv_ids <- paste0(

    "ASV",

    seq_along(seqs)

)

con <- file(

    file.path(

        opt$output,

        "rep_seqs.fasta"

    ),

    "w"

)

for(i in seq_along(seqs)){

    writeLines(

        paste0(

            ">",

            asv_ids[i]

        ),

        con

    )

    writeLines(

        seqs[i],

        con

    )

}

close(con)

############################################################
# Summary
############################################################

cat("\n=====================================\n")

cat("Phyloseq object created successfully\n")

cat("-------------------------------------\n")

cat("Samples :", nsamples(ps), "\n")

cat("ASVs    :", ntaxa(ps), "\n")

cat("Genera  :", ntaxa(ps.genus), "\n")

cat("Output  :", normalizePath(opt$output), "\n")

cat("=====================================\n")

saveRDS(
    seqtab,
    file.path(opt$output, "asv_table.rds")
)

saveRDS(
    taxa,
    file.path(opt$output, "taxonomy_table.rds")
)
