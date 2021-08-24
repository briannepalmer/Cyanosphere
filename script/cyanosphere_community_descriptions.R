# what organisms are present in the cyanosphere 

# load libraries 
library(tidyverse)
library(vegan)
library(ggrepel)

# load data
cyanosphere <- read.csv("data/autometa_assembly_cyanosphere.csv")

pdf("figures/phylum_plot.pdf", width = 11, height = 8.5)
ggplot(cyanosphere, aes(fill = phylum, y = contigs, x = culture_ID)) + geom_bar(position = "fill", stat = "identity")+ theme_bw()  + labs(y = "% abundance") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size = 8), plot.margin = unit(c(0.2, 0.2, 0.2, 2), "cm"))
dev.off()

heatmap.data <- cyanosphere %>% dplyr::select(taxid, contigs, culture_ID)
heatmap.data <- pivot_wider(heatmap.data, names_from = taxid, values_from = contigs, values_fn = sum)
culture.names <-  heatmap.data$culture_ID
heatmap.data <- heatmap.data[-1]
heatmap.data[c(1:299)] <- ifelse(heatmap.data[c(1:299)]>0, 1, 0)
heatmap.data <- sapply(heatmap.data, as.numeric)
heatmap.data[is.na(heatmap.data)] <-0
rownames(heatmap.data) <- culture.names

pdf("figures/taxid heatmap.pdf")
heatmap(heatmap.data, margins = c(10, 11.5))
dev.off()

# are the communities different between cultures
cyanosphere.transformed <- cyanosphere %>% pivot_wider(names_from = culture_ID, values_from = contigs)

otu.all <- cyanosphere.transformed %>% group_by(Host_Species) %>% summarize_at(vars(`Aetokthonos_hydrillicola_B3-Florida`:`Trichormus_sp._ATA11-4-KO1`),list(sum = sum))
otu.all[is.na(otu.all)] <-0
otu.all[c(2:51)] <- ifelse(otu.all[c(2:51)]>0, 1, 0)
otu.only <- otu.all[c(2:51)]
head(otu.only)

metadata <- otu.all[c(1)]

otu.nmds <- metaMDS(otu.only, trace = FALSE, trymax = 100, "bray")
otu.nmds$stress 
adonis(otu.only ~ Host_Species, data = metadata) # didn't work 


# make heat map

heatmap.data <- cyanosphere %>% dplyr::select(phylum, contigs, culture_ID)
heatmap.data <- pivot_wider(heatmap.data, names_from = phylum, values_from = contigs, values_fn = sum)
culture.names <-  heatmap.data$culture_ID
heatmap.data <- heatmap.data[-1]
heatmap.data[c(1:12)] <- ifelse(heatmap.data[c(1:12)]>0, 1, 0)
heatmap.data <- sapply(heatmap.data, as.numeric)
heatmap.data[is.na(heatmap.data)] <-0
rownames(heatmap.data) <- culture.names

pdf("figures/phylum heatmap.pdf")
heatmap(heatmap.data, margins = c(10, 11.5))
dev.off()

# look just at proteobacteria 

proteo <- cyanosphere %>% filter(phylum == "Proteobacteria")
pdf("figures/proteo_plot.pdf", width = 11, height = 8.5)
ggplot(proteo, aes(fill = order, y = contigs, x = culture_ID)) + geom_bar(position = "fill", stat = "identity")+ theme_bw()  + labs(y = "% abundance") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size = 8), plot.margin = unit(c(0.2, 0.2, 0.2, 2), "cm"))
dev.off()

heatmap.data <- proteo %>% dplyr::select(order, contigs, culture_ID)
heatmap.data <- pivot_wider(heatmap.data, names_from = order, values_from = contigs, values_fn = sum)
culture.names <-  heatmap.data$culture_ID
heatmap.data <- heatmap.data[-1]
heatmap.data[c(1:12)] <- ifelse(heatmap.data[c(1:12)]>0, 1, 0)
heatmap.data <- sapply(heatmap.data, as.numeric)
heatmap.data[is.na(heatmap.data)] <-0
rownames(heatmap.data) <- culture.names

pdf("figures/proteo heatmap.pdf")
heatmap(heatmap.data, margins = c(10, 11.5))
dev.off()



### Based on Substrate ###

pdf("figures/substrate_plot.pdf", width = 11, height = 8.5)
ggplot(cyanosphere, aes(fill = phylum, y = contigs, x = Substrate)) + geom_bar(position = "fill", stat = "identity")+ theme_bw()  + labs(y = "% abundance") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size = 8), plot.margin = unit(c(0.2, 0.2, 0.2, 2), "cm"))
dev.off()

heatmap.data <- cyanosphere %>% dplyr::select(phylum, contigs, Substrate)
heatmap.data <- pivot_wider(heatmap.data, names_from = phylum, values_from = contigs, values_fn = sum)
substrate.names <-  heatmap.data$Substrate
heatmap.data <- heatmap.data[-1]
heatmap.data[c(1:12)] <- ifelse(heatmap.data[c(1:12)]>0, 1, 0)
heatmap.data <- sapply(heatmap.data, as.numeric)
heatmap.data[is.na(heatmap.data)] <-0
rownames(heatmap.data) <- substrate.names

pdf("figures/substrate heatmap.pdf")
heatmap(heatmap.data, margins = c(10, 11.5))
dev.off()

heatmap.data <- proteo %>% dplyr::select(order, contigs, Substrate)
heatmap.data <- pivot_wider(heatmap.data, names_from = order, values_from = contigs, values_fn = sum)
heatmap.data <- heatmap.data[-1]
heatmap.data[c(1:12)] <- ifelse(heatmap.data[c(1:12)]>0, 1, 0)
heatmap.data <- sapply(heatmap.data, as.numeric)
heatmap.data[is.na(heatmap.data)] <-0
rownames(heatmap.data) <- substrate.names

pdf("figures/proteo substrate heatmap.pdf")
heatmap(heatmap.data, margins = c(10, 11.5))
dev.off()

### Based on Continent  ###

pdf("figures/continent_plot.pdf", width = 11, height = 8.5)
ggplot(cyanosphere, aes(fill = phylum, y = contigs, x = Continent)) + geom_bar(position = "fill", stat = "identity")+ theme_bw()  + labs(y = "% abundance") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size = 8), plot.margin = unit(c(0.2, 0.2, 0.2, 2), "cm"))
dev.off()

heatmap.data <- cyanosphere %>% dplyr::select(phylum, contigs, Continent)
heatmap.data <- pivot_wider(heatmap.data, names_from = phylum, values_from = contigs, values_fn = sum)
row.names <-  heatmap.data$Continent 
heatmap.data <- heatmap.data[-1]
heatmap.data[c(1:12)] <- ifelse(heatmap.data[c(1:12)]>0, 1, 0)
heatmap.data <- sapply(heatmap.data, as.numeric)
heatmap.data[is.na(heatmap.data)] <-0
rownames(heatmap.data) <- row.names

pdf("figures/Continent heatmap.pdf")
heatmap(heatmap.data, margins = c(10, 11.5))
dev.off()

# Proeteobacteria only 

heatmap.data <- cyanosphere %>% filter(phylum == "Proteobacteria") %>% dplyr::select(order, contigs, Continent)
heatmap.data <- pivot_wider(heatmap.data, names_from = order, values_from = contigs, values_fn = sum)
row.names <-  heatmap.data$Continent 
heatmap.data <- heatmap.data[-1]
heatmap.data[c(1:12)] <- ifelse(heatmap.data[c(1:12)]>0, 1, 0)
heatmap.data <- sapply(heatmap.data, as.numeric)
heatmap.data[is.na(heatmap.data)] <-0
rownames(heatmap.data) <- row.names

pdf("figures/Continent_proteo heatmap.pdf")
heatmap(heatmap.data, margins = c(10, 11.5))
dev.off()


# based on host Order 

heatmap.data <- cyanosphere %>% dplyr::select(phylum, contigs, Host_Order)
heatmap.data <- pivot_wider(heatmap.data, names_from = phylum, values_from = contigs, values_fn = sum)
row.names <-  heatmap.data$Host_Order 
heatmap.data <- heatmap.data[-1]
heatmap.data[c(1:12)] <- ifelse(heatmap.data[c(1:12)]>0, 1, 0)
heatmap.data <- sapply(heatmap.data, as.numeric)
heatmap.data[is.na(heatmap.data)] <-0
rownames(heatmap.data) <- row.names

pdf("figures/host order heatmap.pdf")
heatmap(heatmap.data, margins = c(10, 11.5))
dev.off()

heatmap.data <- cyanosphere %>% filter(phylum == "Proteobacteria") %>% dplyr::select(order, contigs, Host_Order)
heatmap.data <- pivot_wider(heatmap.data, names_from = order, values_from = contigs, values_fn = sum)
row.names <-  heatmap.data$Host_Order 
heatmap.data <- heatmap.data[-1]
heatmap.data[c(1:12)] <- ifelse(heatmap.data[c(1:12)]>0, 1, 0)
heatmap.data <- sapply(heatmap.data, as.numeric)
heatmap.data[is.na(heatmap.data)] <-0
rownames(heatmap.data) <- row.names

pdf("figures/host order bacteria  heatmap.pdf")
heatmap(heatmap.data, margins = c(12, 11.5))
dev.off()

