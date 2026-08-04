
library(phyloseq)
library(ggplot2)
library(ggpubr)

ps1 <- readRDS(file ="~/ENEZA/physeq_gut_august_bacteria.rds" )




#filtering the unwanted sequences
ps2 <- subset_taxa(ps1, (Order!="Chloroplast") | is.na(Order))
ntaxa(ps2)
ps2 <- subset_taxa(ps2, (Genus!="Chloroflexi") | is.na(Genus))
ntaxa(ps2)
ps2<- subset_taxa(ps2, (Family!="Mitochondria") | is.na(Family))
ntaxa(ps2)
ps2<- subset_taxa(ps2, (Kingdom!="Archaea") | is.na(Kingdom))
ntaxa(ps2)
ps2<- subset_taxa(ps2, (Kingdom!="Eukaryota") | is.na(Kingdom))
ntaxa(ps2)






ps3 <- ps2




ps4 <- prune_samples(sample_sums(ps3) >= 20000, ps3)

nsamples(ps4)



set.seed(42)

ps_rarefied <- rarefy_even_depth(
  ps4,
  sample.size = 20000,
  rngseed = 42
)




alpha <- estimate_richness(
  ps_rarefied,
  measures = c("Observed","Shannon","Simpson","Chao1")
)



alpha$SampleID <- rownames(alpha)

meta <- data.frame(sample_data(ps4))
meta$SampleID <- rownames(meta)

alpha <- merge(alpha, meta, by="SampleID")



library(ggpubr)

ggboxplot(
  alpha,
  x = "DISEASE",
  y = "Shannon",
  color = "DISEASE",
  palette = c(
    "Healthy"="#F8766D",
    "Type 2 diabetes"="#00BFC4"
  ),
  add = "jitter"
) +
  stat_compare_means(method="wilcox.test") +
  theme_classic(base_size=14)


