# Cyanosphere functions DRAM 
library(tidyverse)
library(reshape)

dram <- read.csv("data/cyanobacteria_functions.csv")
dram.grouped <- dram %>% dplyr::group_by(Host, genome) %>% summarise(across(X3.Hydroxypropionate.bi.cycle:Sulfur.metabolism..thiosulfate....sulfite), sum)


metadata <- read.csv("data/metadata2024.csv")

dram.grouped <- merge(metadata, dram.grouped)


energy <- read.csv("data/Dram_Use/energy.csv")
energy.long <- energy %>% pivot_longer(cols = Bin.012.fastaAetokthonos_hydrillicola_assembly:Bin.002.fastaTrichormus_sp_assembly, names_to = "genome")

data <- merge(dram.grouped, energy.long, by = "genome")

grouped <- data %>% dplyr::group_by(Host, module, header) %>% dplyr::summarize(value = sum(value))

nitrogen <- grouped %>% filter(header == "Nitrogen")

nitrogen.wide <- nitrogen %>% dplyr::select(Host, module, value) %>% pivot_wider(names_from = module, values_from = value)

nit.mat <- as.matrix(nitrogen.wide)[,-1]
nit.mat <- apply(nit.mat, 2, as.numeric)
rownames(nit.mat) <- nitrogen.wide$Host
colnames(nit.mat) <- c("Ass. nitrate red.", "Comp. nitrification", "Denitrification", "Diss. nitrate red.", "Nitrate Ass.", "Nitrification", "Nitrogen fixation", "NO2+NH4 to N")

heatmap(nit.mat, scale = "row", margins = c(11,6))


library(pheatmap)
nit.hm <- pheatmap(nit.mat, cluster_cols = T, fontsize_row = 6, cellheight = 6, cluster_rows = T, cellwidth = 12)

save_pheatmap_png <- function(x, filename, width=1200, height=2000, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


nit.hm <- pheatmap(nit.mat, 
                   cluster_cols = F, 
                   fontsize_row = 6, 
                   cellheight = 6, 
                   cluster_rows = F, 
                   cellwidth = 12, 
                   color =paletteer_c("ggthemes::Classic Green", 30), 
                   labels_row = as.expression(newnames))

save_pheatmap_png(nit.hm, "figures/cyanobacteria_nitrogen.png")
nit.new <- as.data.frame(nitrogen.wide)
nit.new <- nit.new %>% 
  mutate_at(vars(`Assimilatory nitrate reduction, nitrate => ammonia`:`nitrite + ammonia => nitrogen`), as.numeric)
nit.new[c(2:9)] <- ifelse(nit.new[c(2:9)] > 0, 1, 0)

nit.new$Sums <- rowSums(nit.new[c(2:9)])


carbon <- read.csv("data/Dram_Use/carbon_utilization.csv")
carbon.long <- carbon %>% pivot_longer(cols = Bin.012.fastaAetokthonos_hydrillicola_assembly:Bin.002.fastaTrichormus_sp_assembly, names_to = "genome")

data <- merge(dram.grouped, carbon.long, by = "genome")
data <- data %>% filter(Phylum == "p__Cyanobacteria")

grouped <- data %>% dplyr::group_by(Host, module, header) %>% dplyr::summarize(value = sum(value))

cazy <- grouped %>% filter(header == "CAZY")

cazy.wide <- cazy %>% dplyr::select(Host, module, value) %>% pivot_wider(names_from = module, values_from = value)

pls <- cazy.wide %>% filter(`Polysaccharide Lyases` > 0)

cazy.mat <- as.matrix(cazy.wide)[,-1]
cazy.mat <- apply(cazy.mat, 2, as.numeric)
rownames(cazy.mat) <- cazy.wide$Host
#colnames(cazy.mat) <- c("Ass. nitrate red.", "Comp. nitrification", "Denitrification", "Diss. nitrate red.", "Nitrate Ass.", "nitrification", "nitrogen fixation", "NO2+NH4 to N")

heatmap(cazy.mat, scale = "row", margins = c(11,6))


library(pheatmap)
cazy.hm <-pheatmap(cazy.mat, 
                   cluster_cols = F, 
                   fontsize_row = 6, 
                   cellheight = 6, 
                   cluster_rows = F, 
                   cellwidth = 12, 
                   color =paletteer_c("ggthemes::Classic Green", 30), 
                   labels_row = as.expression(newnames))

save_pheatmap_png <- function(x, filename, width=1200, height=2500, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(cazy.hm, "figures/cyanobacteria_cazy.png")

