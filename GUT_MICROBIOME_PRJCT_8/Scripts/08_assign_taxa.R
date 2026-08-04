library(dada2)

seqtab.all <- readRDS(
    "/home/sam/ENEZA/dada2/results_2/131_seqtab_nochim_merged.rds"
)

taxa <- assignTaxonomy(
    seqtab.all,
    "/home/sam/ENEZA/dada2/results/silva_nr99_v138.1_train_set.fa.gz",
    multithread = TRUE
)

saveRDS(
    taxa,
    "/home/sam/ENEZA/dada2/results_2/taxonomy_merged.rds"
)


####
#
###

silva_species <- "/home/sam/ENEZA/dada2/results/silva_species_assignment_v138.1.fa.gz"

taxa_species <- addSpecies(
  taxa,
  silva_species,
  verbose = TRUE
)

#-----------------------------------------------------------
# Save
#-----------------------------------------------------------

saveRDS(
  taxa_species,
  file.path(results_dir, "taxonomy_species_merged.rds")
)