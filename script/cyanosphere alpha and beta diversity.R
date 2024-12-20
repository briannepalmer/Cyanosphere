library(tidyverse)
library(vegan)
library(ampvis2)
library(ecodist)
library(microeco)
library(magrittr)
library(paletteer)
library(phyloseq)

cyanosphere <- read.csv("data/cyanosphere_clean_2024_lowCheckM_removed.csv")
cyanosphere$Altitude <- as.numeric(cyanosphere$Altitude)
# make ASV column using taxid (this just makes the coding a little smoother later on) 
cyanosphere$MAG <- as.vector(paste0("MAG_", 1:466))
cyanosphere <- cyanosphere %>% filter(Phylum != "p__Cyanobacteria", Phylum != "p__Cyanobacteriota")

pal <- paletteer_d("ggthemes::Superfishel_Stone")
pal <- colorRampPalette(pal)(15)
medpal <- colorRampPalette(pal)(50)
largepal <- colorRampPalette(pal)(75)

cyano.1 <- cyanosphere
cyano.1$value <- c(rep(1,410))
grouped <- cyano.1 %>% dplyr::group_by(Host) %>% summarize(value = sum(value))

transformed.table <- cyanosphere %>% pivot_wider(names_from = Host, values_from = Present) 
transformed.table[c(27:81)][is.na(transformed.table[c(27:81)])] <- 0
transformed.table <- transformed.table %>% group_by(Domain, Phylum, Class, Order, Family, Genus, Species, MAG) %>% summarize_at(vars(Aetokthonos_hydrillicola:Trichotorquatus_andrei), sum)
MAGs <- transformed.table$MAG

otu.table <- transformed.table[c(9:64)]
otu.table[is.na(otu.table)] <-0
otu.table <- as.matrix(otu.table)
rownames(otu.table) <- MAGs

tax.table <- transformed.table[c(1:7)]
tax.table <- as.matrix(tax.table)
rownames(tax.table) <- MAGs

metadata <- read.csv("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/Cyanosphere/data/metadata2024.csv")
md <- as.matrix(metadata)
rownames(md) <- metadata$Host
md <- as.data.frame(md)
env.data <- metadata[,c(17:19)]
env.data$Altitude <- as.numeric(env.data$Altitude)
env.data <- as.matrix(env.data)
rownames(env.data) <- rownames(md)

# write.csv(md, "data/metadata2024.csv")

d <- amp_load(otutable = otu.table,
              metadata = md,
              taxonomy = tax.table)


set.seed(12345)
ps <- phyloseq(otu_table(otu.table, taxa_are_rows = T), tax_table(tax.table), sam_data(md))

melt <- psmelt(ps)
sum = sum(melt$Abundance)
proteobacteria <- melt %>% filter(Phylum %in% c("p__Proteobacteria", "p__Pseudomonadota"), Abundance > 0)
proteobacteria <- sum(proteobacteria$Abundance)
proteobacteria/sum

bacteroidota <- melt %>% filter(Phylum == "p__Bacteroidota", Abundance > 0)
bacteroidota <- sum(bacteroidota$Abundance)
bacteroidota/sum

actinobacteria <- melt %>% filter(Phylum == "p__Actinobacteriota", Abundance > 0)
actinobacteria <- sum(actinobacteria$Abundance)
actinobacteria/sum

brevundomonas <- melt %>% filter(Genus == "g__Brevundimonas", Abundance > 0)
brevundomonas <- sum(brevundomonas$Abundance)
brevundomonas/sum

devosia <- melt %>% filter(Genus == "g__Devosia", Abundance > 0)
devosia <- sum(devosia$Abundance)
devosia/sum

Sphingopyxis <- melt %>% filter(Genus == "g__Sphingopyxis", Abundance > 0)
Sphingopyxis <- sum(Sphingopyxis$Abundance)
Sphingopyxis/sum

48/56

alpha <- amp_alpha_diversity(d)

set.seed(1234)

all.rich <- lm(uniqueOTUs ~ Substrate, data = alpha)
summary(aov(all.rich))

all.rich <- lm(uniqueOTUs ~ Host_Genus, data = alpha)
summary(aov(all.rich))

all.div <- lm(Shannon ~ Host_Genus, data = alpha)
summary(aov(all.div))


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
observed.host <- data_summary(alpha, varname="uniqueOTUs", 
                             groupnames=c("Host_Order"))

library(ggpubr)

ggplot(observed.host, aes(x=Host_Order, y=uniqueOTUs)) + 
  geom_bar(stat="identity", 
           position=position_dodge(), aes(fill = Host_Order)) +
  geom_errorbar(aes(ymin=uniqueOTUs-sd, ymax=uniqueOTUs+sd), width=.2, position=position_dodge(.9)) + theme_bw() + theme(legend.position = "none", text = element_text(size =35)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values = pal) + labs(y = "# of MAGs")

all.div <- lm(Shannon ~ Host_Order, data = alpha)
summary(aov(all.div))
diversity.host <- data_summary(alpha, varname="Shannon", 
                              groupnames=c("Host_Order"))

ggplot(diversity.host, aes(x=Host_Order, y=Shannon)) + 
  geom_bar(stat="identity", 
           position=position_dodge(), aes(fill = Host_Order)) +
  geom_errorbar(aes(ymin=Shannon-sd, ymax=Shannon+sd), width=.2, position=position_dodge(.9)) + theme_bw() + theme(legend.position = "none", text = element_text(size =35)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ scale_fill_manual(values = pal) + labs(y = "Shannon Diversity")


all.rich <- lm(uniqueOTUs ~ Substrate, data = alpha)
summary(aov(all.rich))
all.div <- lm(Shannon ~ Substrate, data = alpha)
summary(aov(all.div))

observed.substrate <- data_summary(alpha, varname="uniqueOTUs", 
                              groupnames=c("Substrate"))

ggplot(observed.substrate, aes(x=Substrate, y=uniqueOTUs)) + 
  geom_bar(stat="identity", 
           position=position_dodge(), aes(fill = Substrate)) +
  geom_errorbar(aes(ymin=uniqueOTUs-sd, ymax=uniqueOTUs+sd), width=.2, position=position_dodge(.9)) + theme_bw() + theme(legend.position = "none", text = element_text(size =35)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values = pal) + labs(x = "Cyanobacteria Habitat", y = "# of MAGs")

shannon.substrate <- data_summary(alpha, varname="Shannon", 
                                   groupnames=c("Substrate"))

ggplot(shannon.substrate, aes(x=Substrate, y=Shannon)) + 
  geom_bar(stat="identity", 
           position=position_dodge(), aes(fill = Substrate)) +
  geom_errorbar(aes(ymin=Shannon-sd, ymax=Shannon+sd), width=.2, position=position_dodge(.9)) + theme_bw() + theme(legend.position = "none", text = element_text(size =35)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values = pal) + labs(x = "Cyanobacteria Habitat")

observed.habitat <- data_summary(alpha, varname="uniqueOTUs", 
                                   groupnames=c("Habitat2"))

ggplot(observed.habitat, aes(x=Habitat2, y=uniqueOTUs)) + 
  geom_bar(stat="identity", 
           position=position_dodge(), aes(fill = Habitat2)) +
  geom_errorbar(aes(ymin=uniqueOTUs-sd, ymax=uniqueOTUs+sd), width=.2, position=position_dodge(.9)) + theme_bw() + theme(legend.position = "none", text = element_text(size =35)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values = pal) + labs(x = "Cyanobacteria Habitat", y = "# of MAGs")


shannon.substrate <- data_summary(alpha, varname="Shannon", 
                                  groupnames=c("Substrate"))

ggplot(shannon.substrate, aes(x=Substrate, y=Shannon)) + 
  geom_bar(stat="identity", 
           position=position_dodge(), aes(fill = Substrate)) +
  geom_errorbar(aes(ymin=Shannon-sd, ymax=Shannon+sd), width=.2, position=position_dodge(.9)) + theme_bw() + theme(legend.position = "none", text = element_text(size =35)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values = pal) + labs(x = "Cyanobacteria Habitat")


# alpha diversity plots 

png("figures/richness_hostorder.png", res = 300, units = "cm", height = 16, width = 16)
ggplot(alpha, aes(x = uniqueOTUs, y = Host, fill = Host_Order)) +
  geom_bar(stat="identity") + theme_bw() +
  facet_grid(Host_Order ~ ., scale = "free", space = "free_y") + 
  theme(panel.spacing = unit(0.1, "lines"), 
        strip.background = element_blank(),
        strip.placement = "inside",
        legend.position = "none", 
        strip.text.y.right = element_text(angle = 0)) + scale_fill_manual(values = pal) + labs(x = "# MAGS", y = "Cyanobacteria")
dev.off()

alpha2 <- alpha                                                 # Replicate original data
alpha2$Name <- factor(alpha2$Name, levels = alpha2$Name[order(alpha2$uniqueOTUs)]) # Factor levels in increasing order

# remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)
richness <- ggplot(data = alpha2, aes(x = Name, y = uniqueOTUs, fill = Host_Order, pattern = Substrate)) +
  geom_bar_pattern(position = "dodge", stat = "identity",
                   color = "grey", 
                   pattern_fill = "grey",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.03,
                   pattern_key_scale_factor = 1)+ 
  scale_pattern_manual(values = c("fresh water" = "stripe", "subaerial" = "circle", "soil" = "none", "hydroterrestrial" = "crosshatch")) +
  scale_fill_manual(values = pal) +
  labs(x = "", y = "# of MAGs", pattern = "Habitat") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"), title = "Cyanobacteria Order")) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1, face = "italic"), text=element_text(size=18), legend.position = "none") + ggtitle("A.") + facet_grid(.~Substrate, scales = "free")

richness

richness <- ggplot(data = alpha2, aes(x = Name, y = uniqueOTUs, fill = Host_Order))+ geom_bar(stat = "identity", color = "gray") +
  scale_fill_manual(values = pal) +
  labs(x = "", y = "# of MAGs") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1, face = "italic"), text=element_text(size=18), legend.position = "none") + ggtitle("A.") + facet_grid(.~Substrate, scales = "free", space = "free")

richness




alpha3 <- alpha       # Replicate original data
alpha3$Shannon <- as.numeric(alpha3$Shannon)
alpha3$Name <- factor(alpha3$Name, levels = alpha3$Name[order(alpha3$Shannon)]) # Factor levels in increasing order

diversity <- ggplot(data = alpha3, aes(x = Name, y = Shannon, fill = Host_Order, pattern = Substrate)) +
  geom_bar_pattern(position = "dodge", stat = "identity",
                   color = "grey", 
                   pattern_fill = "grey",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6)+ 
  scale_pattern_manual(values = c("fresh water" = "stripe", "subaerial" = "circle", "soil" = "none", "hydroterrestrial" = "crosshatch")) +
  scale_fill_manual(values = pal) +
  labs(x = "Cyanobacteria Host", y = "Shannon Diversity", pattern = "Habitat") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"), title = "Cyanobacteria Order")) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1, face = "italic"), text=element_text(size=18)) + ggtitle("B.")


diversity <- ggplot(data = alpha3, aes(x = Name, y = Shannon, fill = Host_Order))+ geom_bar(stat = "identity", color = "gray") +
  scale_fill_manual(values = pal) +
  labs(x = "", y = "Shannon Diversity") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1, face = "italic"), text=element_text(size=18)) + ggtitle("A.") + facet_grid(.~Substrate, scales = "free", space = "free")

diversity


library(patchwork)

png("figures/Manuscript/Richnes_Diversity.png", res = 300, units = "px", height = 3500, width = 5500)
richness/diversity
dev.off()

# beta diversity calculations 
# use Genus level 

ps.genus <- tax_glom(ps, "Genus")
ps.genus <- subset_taxa(ps.genus, Genus != "g__")

library(file2meco)

df <- phyloseq2meco(physeq = ps.genus)
df
df$filter_pollution(taxa = c("mitochondria", "chloroplast"))
df$sample_sums() %>% range
View(df$otu_table)
df$tidy_dataset()

set.seed(12345)
t1 <- trans_alpha$new(dataset = df) 
df$cal_betadiv()
t1 <- trans_beta$new(dataset = df, measure = "bray", group = "Host_Order")
t1$cal_manova(manova_all = TRUE)
t1$res_manova # 0.13
t1 <- trans_beta$new(dataset = df, measure = "bray", group = "Host_Family")
t1$cal_manova(manova_all = TRUE)
t1$res_manova # 0.272
t1 <- trans_beta$new(dataset = df, measure = "bray", group = "Host_Genus")
t1$cal_manova(manova_all = TRUE)
t1$res_manova # 0.04


# habitats differ!
t1 <- trans_beta$new(dataset = df, measure = "bray", group = "Substrate")
t1$cal_manova(manova_all = TRUE)
t1$res_manova # 0.001

t1 <- trans_beta$new(dataset = df, measure = "bray", group = "Habitat")
t1$cal_manova(manova_all = TRUE)
t1$res_manova # 0.0041

t1 <- trans_beta$new(dataset = df, measure = "bray", group = "Habitat2")
t1$cal_manova(manova_all = TRUE)
t1$res_manova # 0.017

t1 <- trans_beta$new(dataset = df, measure = "bray", group = "Habitat3")
t1$cal_manova(manova_all = TRUE)
t1$res_manova # 0.011

df2 <- clone(df)
df2$sample_table <- as.data.frame(df2$sample_table) %>% filter(Parentmaterial != "")
df2$tidy_dataset()
t1 <- trans_beta$new(dataset = df2, measure = "bray", group = "Parentmaterial")
t1$cal_manova(manova_all = TRUE)
t1$res_manova # 0.05

t1 <- trans_beta$new(dataset = df, measure = "bray", group = "Continent")
t1$cal_manova(manova_all = TRUE)
t1$res_manova # 0.004

t1 <- trans_beta$new(dataset = df, measure = "bray", group = "Location")
t1$cal_manova(manova_all = TRUE)
t1$res_manova # 0.002



# PCOA no environmental data 
df$cal_betadiv()
t1 <- trans_beta$new(dataset = df, group = "Habitat2", measure = "bray")
t1$cal_ordination(method = "PCoA")

t1$plot_ordination(plot_color = "Host_Order", plot_type = c("point"), point_size = 6) + theme_bw() + scale_color_manual(values = pal, name = "Cyanobacteria Host") + theme(text = element_text(size = 35)) 

t1$plot_ordination(plot_color = "Habitat2", plot_type = c("ellipse", "point"), centroid_segment_linetype = 1, point_size = 6) + scale_color_manual(values = pal)+ scale_fill_manual(values = pal) + theme_bw()

t1$plot_ordination(plot_color = "Continent", plot_type = c("ellipse", "point"), centroid_segment_linetype = 1, point_size = 6, plot_shape = "Substrate", shape_values = c(16, 17, 15, 18)) + scale_color_manual(values = pal)+ scale_fill_manual(values = pal) + theme_bw()

t1$plot_ordination(plot_color = "Continent", plot_type = c("point", "chull")) + theme_bw()+ scale_color_manual(values = pal)+ scale_fill_manual(values = pal)

df2 <- clone(df)
df2$sample_table <- as.data.frame(df2$sample_table) %>% filter(Parentmaterial != "")
df2$tidy_dataset()
df2$cal_betadiv()
t1 <- trans_beta$new(dataset = df2, group = "Parentmaterial", measure = "bray")
t1$cal_ordination(method = "PCoA")
t1$plot_ordination(plot_color = "Parentmaterial", plot_type = c("point", "chull")) + theme_bw()
t1$cal_manova(manova_all = TRUE)
t1$res_manova # 0.025


# use bray-curtis distance for dbRDA
env.data <- as.data.frame(env.data)
env.data$Altitude <- as.numeric(env.data$Altitude)

set.seed(12345)
## HOST ORDER
t1 <- trans_env$new(dataset = df, add_data = data.frame(env.data), complete_na = TRUE)
t1$cal_diff(group = "Substrate", method = "anova")
t1$cal_ordination(method = "dbRDA", use_measure = "bray")
t1$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.5)
# t1$res_rda_trans is the transformed result for plotting
t1$plot_ordination(plot_color = "Habitat2", point_size = 6) + scale_color_manual(values = pal)+ scale_fill_manual(values = pal) + theme(text = element_text(size = 35))

t1$cal_ordination_anova()
t1$res_ordination_terms # precipitation 0.019
t1$res_ordination_axis
t1$cal_ordination_envfit()
t1$res_ordination_envfit

t1$cal_ordination(method = "RDA", taxa_level = "Genus")
# As the main results of RDA are related with the projection and angles between different arrows,
# we adjust the length of the arrow to show them clearly using several parameters.
t1$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)
# t1$res_rda_trans is the transformed result for plot
t1$plot_ordination(plot_color = "Host_Order")+ scale_color_manual(values = pal)+ scale_fill_manual(values = pal)

t1$cal_ordination(method = "dbRDA", use_measure = "bray")

t1$plot_ordination(plot_color = "Substrate", plot_type = c("ellipse", "point"), centroid_segment_linetype = 1, point_size = 6, env_text_size =5) + scale_color_manual(values = pal)+ scale_fill_manual(values = pal) + theme(text = element_text(size = 35)) 


soil <- clone(df)
soil$sample_table <- as.data.frame(soil$sample_table) %>% filter(Substrate == "soil")
soil$tidy_dataset()
soil$cal_betadiv()
t1 <- trans_beta$new(dataset = soil, group = "Habitat2", measure = "bray")
t1$cal_ordination(method = "PCoA")
t1$plot_ordination(plot_color = "Habitat2", plot_type = c("point", "chull"),centroid_segment_linetype = 1, point_size = 6, env_text_size = 5) + theme_bw()+ scale_color_manual(values = pal)+ scale_fill_manual(values = pal)
t1$cal_manova(manova_all = TRUE)
t1$res_manova # 0.025


dataset1 <- df
dataset1
dataset1$filter_pollution(taxa = c("mitochondria", "chloroplast"))

dataset1 <- dataset1$merge_samples(use_group = "Substrate")
t1 <- trans_venn$new(dataset1, ratio = "numratio")
t1$plot_venn()


df$cal_betadiv()
t1 <- trans_beta$new(dataset = df, group = "Substrate", measure = "bray")
t1$cal_ordination(method = "PCoA")
t1$plot_ordination(plot_type = c("point"), plot_color = "Host_Order",point_size = 5) + theme(panel.grid = element_blank()) + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + theme_bw() + facet_grid(.~Habitat) + scale_color_manual(values = pal)

df$cal_betadiv()
t1 <- trans_beta$new(dataset = df, group = "Substrate", measure = "bray")
t1$cal_ordination(method = "PCoA")
t1$plot_ordination(plot_type = c("point"), plot_color = "Host_Order", point_size = 5) + theme(panel.grid = element_blank()) + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + theme_bw() + facet_grid(.~Substrate) + scale_color_manual(values = pal)

t1<- trans_abund$new(dataset = df, taxrank = "Genus", ntaxa = 3, groupmean = "Name")


png("figures/three_abundant_genera_ppt.png", width = 30, height = 15, units = "cm", res = 300)
t1$plot_bar(others_color = "grey70", xtext_keep = T, legend_text_italic = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust= 1),text = element_text(size = 10)) + scale_fill_manual(values = c("grey70", "#6388B4","#FA9C43", "#8BB464", "#BC8580"))
dev.off()


t1 <- trans_abund$new(dataset = df, taxrank = "Genus", ntaxa = 40)
t1$plot_heatmap(xtext_keep = T, legend_text_italic = FALSE) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1))


png("figures/Manuscript/Heatmap.png", width = 3000, height = 2500, units = "px", res = 300)
t1 <- trans_abund$new(dataset = df, taxrank = "Genus", ntaxa = 50,groupmean = "Name")
t1$plot_heatmap(xtext_keep = T, legend_text_italic = FALSE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1,face = "italic"), axis.text.y = element_text(face = "italic"), text = element_text(size = 18)) +scale_fill_gradient2(low = "#6388B4",high = "#F07366", mid = "white")
dev.off()

t1 <- t1 <- trans_abund$new(dataset = df, taxrank = "Order", ntaxa = 50)
t1$plot_bar(others_color = "grey70", facet = "Substrate", xtext_keep = T, legend_text_italic = FALSE) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1))


t1 <- trans_abund$new(dataset = df, taxrank = "Phylum", groupmean = "Name")

phylum.plot <- t1$plot_bar(others_color = "grey70", xtext_keep = T, legend_text_italic = FALSE) + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust= 1))+ theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "italic"), text = element_text(size = 18)) + scale_fill_manual(values = c("#767676","#C1B24A","#6388B4", "#FA9C43", "#F07366", "#B69EA0", "#80BDBC", "#8BB464", "#C7A061","#BC8580", "#5CB092", "#B8A195","#ADAFA6"))+ ggtitle("A.")

phylum.plot


df.proteo <- clone(df)
df.proteo$tax_table <- as.data.frame(df.proteo$tax_table) %>% filter(Phylum == "p__Proteobacteria")
df.proteo$tidy_dataset()

t1 <- trans_abund$new(dataset = df.proteo, taxrank = "Order", groupmean = "Name")


proteo.plot <- t1$plot_bar(others_color = "grey70", xtext_keep = T, legend_text_italic = FALSE) + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust= 1))+ theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "italic"), text = element_text(size = 18)) + scale_fill_manual(values = c("#767676","#C1B24A","#6388B4", "#FA9C43", "#F07366", "#B69EA0", "#80BDBC", "#8BB464", "#C7A061","#BC8580", "#5CB092", "#B8A195","#ADAFA6"))+ ggtitle("B.")

library(patchwork)
png("figures/Manuscript/Abundance.png", width = 5000, height = 4000, units = "px", res = 300)
phylum.plot/proteo.plot
dev.off()




png("figures/rel_abunda_order.png", width = 4000, height = 2000, units = "px")
t1 <- trans_abund$new(dataset = df, taxrank = "Order", ntaxa = 10)
t1$plot_bar(others_color = "grey70", facet = "Substrate", xtext_keep = T, legend_text_italic = T, xtext_size = 25, ytitle_size = 20, strip_text = 20) + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust= 1),text = element_text(size = 35)) + scale_fill_manual(values = c("grey70", medpal))
dev.off()

t1 <- trans_diff$new(dataset = df, method = "lefse", group = "Substrate", alpha = 0.01, lefse_subgroup = NULL)
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
t1$plot_diff_bar(threshold = 4)


t1 <- trans_diff$new(dataset = df, method = "anova", group = "Host_Order", taxa_level = "Genus", filter_thres = 0.001)
t1$plot_diff_abund(use_number = 1:10, add_sig = T, coord_flip = F, color_values = pal)

df$cal_betadiv(method = "bray")
t1 <- trans_env$new(dataset = df, add_data = env.data, complete_na = TRUE)
# use add_abund_table parameter to add the extra data table
t1$cal_cor(add_abund_table = df$alpha_diversity)
# use pH and bray-curtis distance
# add correlation statistics
t1$plot_scatterfit(
  x = "Annual_Precipitation", 
  y = df$beta_diversity$bray[rownames(t1$data_env), rownames(t1$data_env)], 
  type = "cor",
  point_size = 3, point_alpha = 0.1, 
  label.x.npc = "center", label.y.npc = "bottom", 
  x_axis_title = "Euclidean distance of pH", 
  y_axis_title = "Bray-Curtis distance"
)


t1$cal_cor(use_data = "Genus", p_adjust_method = "fdr", p_adjust_type = "Env")
t1$plot_cor()



