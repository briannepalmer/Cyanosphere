# what organisms are present in the cyanosphere 

# load libraries 
library(tidyverse)
library(calecopal)

# load data
cyanosphere <- read.csv("data/autometa_assembly_cyanosphere.csv")

pdf("figures/phylum_plot.pdf", width = 11, height = 8.5)
ggplot(cyanosphere, aes(fill = phylum, y = contigs, x = culture_ID)) + geom_bar(position = "fill", stat = "identity")+ theme_bw()  + labs(y = "% abundance") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size = 8), plot.margin = unit(c(0.2, 0.2, 0.2, 2), "cm"))
dev.off()


# make heat map

heatmap.data <- cyanosphere %>% dplyr::select(phylum, contigs, culture_ID)
heatmap.data <- pivot_wider(heatmap.data, names_from = phylum, values_from = contigs, values_fn = sum)
culture.names <-  heatmap.data$culture_ID
heatmap.data <- heatmap.data[-1]
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
heatmap.data <- sapply(heatmap.data, as.numeric)
heatmap.data[is.na(heatmap.data)] <-0
rownames(heatmap.data) <- culture.names

pdf("figures/proteo heatmap.pdf")
heatmap(heatmap.data, margins = c(10, 11.5))
dev.off()
