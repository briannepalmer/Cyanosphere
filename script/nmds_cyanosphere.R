# load libraries 
library(tidyverse)
library(vegan)
library(ggrepel)


# load data

cyanosphere <- read.csv("data/autometa_assembly_cyanosphere.csv")
cyanosphere <- cyanosphere %>% dplyr::rename(sample_name = ï..sample_name, OTU_ID = taxid)
cyanosphere$contigs <- as.numeric(cyanosphere$contigs)
cyanosphere.transformed <- cyanosphere %>% pivot_wider(names_from = culture_ID, values_from = contigs)
#metadata <- cyanosphere[c(2,13, 20,21)]
##metadata <- metadata %>% group_by(phylum) %>% summarize(Habitat = unique(Habitat), Location = unique(Location), culture_ID = unique(culture_ID))

cyano.all <-cyanosphere.transformed
cyanp.all[c(50:99)] <- sapply(cyano.all[c(50:99)],as.numeric)
cyano.all[is.na(cyano.all)] <-0
otu.all <- cyano.all %>% group_by(OTU_ID, Substrate, Host_Genus, Continent, phylum) %>% summarize_at(vars(`Aetokthonos_hydrillicola_B3-Florida`:`Trichormus_sp._ATA11-4-KO1`),list(sum = sum))

otu.only <- otu.all[c(6:55)]
head(otu.only)
#phylum.transposed <- as.data.frame(t(as.matrix(phylum)))
#rownames(phylum.transposed) <- colnames(phylum)
#colnames(phylum.transposed) <- paste("S", 1:77, sep = "")

metadata <- otu.all[-c(6:55)]

otu.nmds <- metaMDS(otu.only, trace = FALSE, trymax = 100, "bray")
otu.nmds$stress

#build a data frame with NMDS coordinates and metadata 
MDS1 = otu.nmds$points[,1] # adds the points of the NMDS 1 dimmension
MDS2 = otu.nmds$points[,2] # adds the points of the NMDS 2 dimmension

NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2)
NMDS = cbind(NMDS, metadata)

ggplot(NMDS, aes(x = MDS1, y = MDS2)) + geom_point(aes(color = Substrate),size = 2, alpha = 0.3) + 
  theme_bw() + theme(legend.position = "top")

as.numeric(cyanosphere.transformed[c(31:49)])

env <- cyanosphere.transformed %>% group_by(OTU_ID, Substrate, Host_Genus, Continent, phylum) %>% summarize_at(vars(bio1:bio19),list(mean = mean))

env.fit <- envfit(otu.nmds, env , permutations = 999, na.rm = TRUE)
env.fit$vectors
arrow <- data.frame(env.fit$vectors$arrows, R = env.fit$vectors$r, P = env.fit$vectors$pvals)
arrow$FG <-  rownames(arrow)
arrow.p <- filter(arrow, P<=0.05)

pdf("figures/nmds with bioclim data.pdf")
ggplot(data=NMDS, aes(x = MDS1, y = MDS2)) + geom_point(data=NMDS, aes(MDS1, MDS2, shape = Continent, color = phylum),position=position_jitter(.1)) + theme_minimal()+ geom_segment(data=arrow.p, aes(x=0, y=0, xend=NMDS1, yend=NMDS2), arrow=arrow(length=unit(.2, "cm")*arrow.p$R)) + geom_label_repel(data=arrow.p, aes(x=NMDS1, y=NMDS2,  label = FG),  label.padding = 0.1, label.size = 0.1, size =3, max.overlaps = 20)
dev.off()

