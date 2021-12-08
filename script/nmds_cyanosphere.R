# load libraries 
library(tidyverse)
library(vegan)
library(ggrepel)
library(pairwiseAdonis)


# load data

cyanosphere <- read.csv("data/autometa_assembly_cyanosphere.csv")
cyanosphere <- cyanosphere %>% dplyr::rename(sample_name = ï..sample_name, OTU_ID = taxid)

cyanosphere.transformed <- cyanosphere %>% pivot_wider(names_from = OTU_ID, values_from = contigs)

cyano.all <-cyanosphere.transformed
cyano.all[c(50:348)] <- sapply(cyano.all[c(50:348)],as.numeric)
cyano.all[is.na(cyano.all)] <-0
otu.all <- cyano.all %>% group_by(culture_ID, Habitat, Location, Continent, Substrate, Host_Phylum, Host_Class, Host_Order, Host_Family, Host_Genus) %>% summarize_at(vars(bio1:`159191`),list(mean = mean))


# create presence absence matrix 
otu.all[c(30:328)] <- ifelse(otu.all[c(30:328)]>0, 1, 0)


otu.only <- otu.all[c(30:328)]
head(otu.only)
metadata <- otu.all[-c(30:328)]

otu.nmds <- metaMDS(otu.only, trace = FALSE, trymax = 100, "bray")
otu.nmds$stress #stress is near zero
stressplot(otu.nmds)

adonis(otu.only ~ Substrate, metadata) 

pairwise.adonis(otu.only, factors = as.vector(metadata$Host_Order))


#build a data frame with NMDS coordinates and metadata 
MDS1 = otu.nmds$points[,1] # adds the points of the NMDS 1 dimmension
MDS2 = otu.nmds$points[,2] # adds the points of the NMDS 2 dimmension

NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2)
NMDS = cbind(NMDS, metadata)

ggplot(NMDS, aes(x = MDS1, y = MDS2)) + geom_point(aes(color = Habitat),size = 2, alpha = 0.3) + 
  theme_bw() + theme(legend.position = "top") #+ stat_ellipse(aes(color = Host_Order))

env <- cyanosphere.transformed %>% group_by(OTU_ID, Substrate, Continent, phylum, Host_Order) %>% summarize_at(vars(Lat,Altitude, bio1:bio19),list(mean = mean))
env.fit <- envfit(otu.nmds, env , permutations = 999, na.rm = TRUE)
env.fit$vectors
arrow <- data.frame(env.fit$vectors$arrows, R = env.fit$vectors$r, P = env.fit$vectors$pvals)
arrow$FG <-  rownames(arrow)
arrow.p <- filter(arrow, P<=0.05)

vcov(env.fit)
adonis <- adonis(otu.only ~ Lat_mean + Altitude_mean + bio1_mean, metadata) 
coef(adonis)

pdf("figures/nmds with mean temp and mean precip data.pdf", height = 10, width = 15)
ggplot(data=NMDS, aes(x = MDS1, y = MDS2)) + geom_point(data=NMDS, aes(MDS1, MDS2, shape = Continent, color = Host_Order),position=position_jitter(.1)) + theme_minimal()+ geom_segment(data=arrow.p, aes(x=0, y=0, xend=NMDS1, yend=NMDS2), arrow=arrow(length=unit(.2, "cm")*arrow.p$R)) + geom_label_repel(data=arrow.p, aes(x=NMDS1, y=NMDS2,  label = FG),  label.padding = 0.1, label.size = 0.1, size =3, max.overlaps = 30)
dev.off()

#####

# look at the significance of the climate variables 

model.data <- cyanosphere %>% dplyr::select(contigs, phylum, Lat, Altitude, Habitat:bio19)
model.data <- model.data[-9]

m1 <- lm(contigs ~ ., data = model.data)
car::Anova(m1)
plot(m1)

m2 <- lm(contigs ~ Lat + Altitude + bio1 + bio2 + bio3 + bio5 + bio6  + bio8 + bio9 + bio10 + bio11+ bio15, model.data)
car::Anova(m2)
plot(m2)

ggplot(cyanosphere, aes(x = bio8, y = contigs)) + geom_point(aes(color = phylum)) + geom_smooth(method = 'glm') + facet_wrap(.~phylum) + theme_bw() + theme(legend.position = "none")
