# Cyanosphere functions DRAM 
library(tidyverse)
library(reshape)

energy <- read.csv("data/Dram_Use/energy.csv")
energy.long <- energy %>% pivot_longer(cols = Bin.012.fastaAetokthonos_hydrillicola_assembly:Bin.002.fastaTrichormus_sp_assembly, names_to = "genome")

taxonomy <- read.csv("data/cyanosphere_clean_2024_lowCheckM_removed.csv")
taxonomy <- taxonomy %>% filter(Phylum != "p__Cyanobacteria")

data <- merge(taxonomy, energy.long, by = "genome")

grouped <- data %>% dplyr::group_by(Name, module, header) %>% dplyr::summarize(value = sum(value))

nitrogen <- grouped %>% filter(header == "Nitrogen")

nitrogen.wide <- nitrogen %>% dplyr::select(Name, module, value) %>% pivot_wider(names_from = module, values_from = value)

nit.mat <- as.matrix(nitrogen.wide)[,-1]
nit.mat <- apply(nit.mat, 2, as.numeric)
rownames(nit.mat) <- nitrogen.wide$Name
colnames(nit.mat) <- c("Ass. nitrate red.", "Comp. nitrification", "Denitrification", "Diss. nitrate red.", "Nitrate Ass.", "Nitrification", "Nitrogen fixation", "NO2+NH4 to N")

heatmap(nit.mat, scale = "row", margins = c(11,6))

library(pheatmap)

newnames <- lapply(
  rownames(nit.mat),
  function(x) bquote(italic(.(x))))

nit.hm <- pheatmap(nit.mat, 
                   cluster_cols = F, 
                   fontsize_row = 6, 
                   cellheight = 6, 
                   cluster_rows = F, 
                   cellwidth = 12, 
                   color =paletteer_c("ggthemes::Blue-Teal", 30), 
                   labels_row = as.expression(newnames))


save_pheatmap_png <- function(x, filename, width=1200, height=2000, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(nit.hm, "figures/cyanosphere_nitrogen.png")


ps <- grouped %>% filter(header == "Photosynthesis")

ps.wide <- ps %>% dplyr::select(Host, module, value) %>% pivot_wider(names_from = module, values_from = value)

ps.mat <- as.matrix(ps.wide)[,-1]
ps.mat <- apply(ps.mat, 2, as.numeric)
rownames(ps.mat) <- ps.wide$Host
#colnames(ps.mat) <- c("Ass. nitrate red.", "Comp. nitrification", "Denitrification", "Diss. nitrate red.", "Nitrate Ass.", "nitrification", "nitrogen fixation", "NO2+NH4 to N")


library(pheatmap)
pheatmap(ps.mat, cluster_cols = T, fontsize_row = 6, cellheight = 6, cluster_rows = T, cellwidth = 12)


carbon <- read.csv("data/Dram_Use/carbon_utilization.csv")
carbon.long <- carbon %>% pivot_longer(cols = Bin.012.fastaAetokthonos_hydrillicola_assembly:Bin.002.fastaTrichormus_sp_assembly, names_to = "genome")

taxonomy <- read.csv("data/cyanosphere_clean_2024_lowCheckM_removed.csv")
taxonomy <- taxonomy %>% filter(Phylum != "p__Cyanobacteria")

data <- merge(taxonomy, carbon.long, by = "genome")

grouped <- data %>% dplyr::group_by(Name, module, header) %>% dplyr::summarize(value = sum(value))

cazy <- grouped %>% filter(header == "CAZY")

cazy.wide <- cazy %>% dplyr::select(Name, module, value) %>% pivot_wider(names_from = module, values_from = value)

cazy.mat <- as.matrix(cazy.wide)[,-1]
cazy.mat <- apply(cazy.mat, 2, as.numeric)
rownames(cazy.mat) <- cazy.wide$Name
heatmap(cazy.mat, scale = "row", margins = c(11,6))

library(pheatmap)
cazy.hm <- pheatmap(cazy.mat,  
                    cluster_cols = F, 
                    fontsize_row = 6, 
                    cellheight = 6, 
                    cluster_rows = F, 
                    cellwidth = 12, 
                    color = paletteer_c("ggthemes::Blue-Teal", 30), 
                    labels_row = as.expression(newnames))

save_pheatmap_png <- function(x, filename, width=1200, height=2500, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(cazy.hm, "figures/cyanosphere_cazy.png")



carbon <- read.csv("data/Dram_Use/carbon_utilization.csv")
carbon.long <- carbon %>% pivot_longer(cols = Bin.012.fastaAetokthonos_hydrillicola_assembly:Bin.002.fastaTrichormus_sp_assembly, names_to = "genome")

taxonomy <- read.csv("data/cyanosphere_clean_2024_lowCheckM_removed.csv")
taxonomy <- taxonomy %>% filter(Phylum != "p__Cyanobacteria")

data <- merge(taxonomy, carbon.long, by = "genome")

grouped <- data %>% dplyr::group_by(Name, module, header, gene_description, subheader, Genus) %>% dplyr::summarize(value = sum(value))

cazy <- grouped %>% filter(header == "CAZY")

eps <- cazy %>% filter(module == "Polysaccharide Lyases")
sphingomonas <- eps %>% filter(Genus == "g__Sphingomonas") %>% group_by(Name) %>% summarize(value = sum(value))

brevundimonas <- eps %>% filter(Genus == "g__Brevundimonas") %>% group_by(Name) %>% summarize(value = sum(value))

PL_alg <- eps %>% filter(gene_description %in% c("PL7 poly(beta-mannuronate) lyase / M-specific alginate lyase (EC 4.2.2.3); alpha-L-guluronate lyase / G-specific alginate lyase (EC 4.2.2.11); poly-(MG)-lyase / MG-specific alginate lyase (EC 4.2.2.-); endo-beta-1,4-glucuronan lyase (EC 4.2.2.14)", "PL6 alginate lyase (EC 4.2.2.3); chondroitinase B (EC 4.2.2.19); MG-specific alginate lyase (EC 4.2.2.-); poly(alpha-L-guluronate) lyase / G-specific alginate lyase (EC 4.2.2.11);", "PL5 alginate lyase (EC 4.2.2.3).", "PL36 poly(beta-mannuronate) lyase / M-specific alginate lyase (EC 4.2.2.3)", "PL17 alginate lyase (EC 4.2.2.3); oligoalginate lyase (EC 4.2.2.26)", "PL15 oligo-alginate lyase (EC 4.2.2.-); alginate lyase (EC 4.2.2.3)", "PL14 alginate lyase (EC 4.2.2.3); exo-oligoalginate lyase (EC 4.2.2.26); beta-1,4-glucuronan lyase (EC 4.2.2.14)"))

alg.ly <- merge(PL_alg, taxonomy)
alg.ly <- alg.ly %>% group_by(Name, module, gene_description, subheader, Host_Order, value) %>% dplyr::summarise()

ggplot(alg.ly %>% filter(value > 0), aes(x = Name, y = value, fill = Host_Order)) + geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"), text = element_text(size = 18)) + labs(x = "Cyanobacteria Host", y = "Count of Alginate Lyases") + scale_fill_manual(values = pal)

    