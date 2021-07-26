# load libraries 
library(tidyverse)
library(phyloseq)
library(vegan)


# load data

cyanosphere <- read.csv("data/autometa_assembly_cyanosphere.csv")
cyanosphere <- cyanosphere %>% dplyr::rename(sample_name = ï..sample_name, OTU_ID = taxid)
cyanosphere.transformed <- cyanosphere %>% pivot_wider(names_from = culture_ID, values_from = contigs)
cyanosphere.transformed[is.na(cyanosphere.transformed)] <-0
metadata <- cyanosphere[c(1:11, 20, 21)]
genus.only <- cyanosphere.transformed[c(20:69)]

# create OTU matrix 
otu <- as.matrix(genus.only)
OTU = otu_table(otu, taxa_are_rows = TRUE)


# create taxonomy matrix
taxmat <- as.matrix(cyanosphere.transformed[c(10:16)])
TAX = tax_table(taxmat)

# create sample data 
sampleData <- sample_data(metadata)
sample_names(sampleData)

# create phyloseq object 
physeq = phyloseq(OTU, TAX, sampleData)

library(ape)
random_tree = rtree(ntaxa(physeq), rooted = TRUE, tip.label = taxa_names(physeq))
plot(random_tree)

physeq1 = merge_phyloseq(physeq, sampleData, random_tree)
physeq1

physeq2 = phyloseq(OTU, TAX, sampleData, random_tree)
physeq2

identical(physeq1, physeq2)

plot_tree(physeq1, color = "Habitat", label.tips = "taxa_names", ladderize = "left", text.size = 1.5)
plot_heatmap(physeq1, taxa.label = "Phylum")
plot_bar(physeq1, fill = "phylum")

plot_ordination(physeq1, ordinate(OTU, "MDS"), color = "Habitat") + geom_point()

genus.only %>% metaMDS(trace = "F") %>% ordiplot(type = "none") %>% text("sites")
data(varespec)
