### Oneida eDNA metabarcoding analysis
### This code compares Oneida Lake species inventories from capture gears (electrofishing, fyke netting, gillnetting, and seining) survey data, 
### historical records, and eDNA metabarcoding datasets 
### Created by K. Andres and T. Lambert

#### Clear the working environment and load required packages ####
rm(list = ls())
library(plyr)
library(dplyr)
library(RColorBrewer)
library(stringr)
library(UpSetR)
library(ggplot2)

# Load datasets
historical_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/oneida_historical_species_list.csv", header = TRUE)
eDNA_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/eDNA_dat.csv", header = TRUE)
ef_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/ef_dat.csv", header = TRUE)
fyke_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/fyke_dat.csv", header = TRUE)
gillnet_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/gillnet_dat.csv", header = TRUE)
seine_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/seine_dat.csv", header = TRUE)
all_methods_abundance <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/all_methods_abundance.csv", header = TRUE)
all_methods_presence <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/all_methods_presence.csv", header = TRUE)

# Make modifications to gear files as necessary
eDNA_dat <- subset(eDNA_dat, grepl("G", eDNA_dat$Site)) # Restrict eDNA samples to the "green" (large-volume) samples
ef_dat <- ddply(ef_dat, "Site", numcolwise(sum)) # Combine replicates per site (2)              
fyke_dat <- ddply(fyke_dat, "Site", numcolwise(sum)) # Combine replicates per site (2)

###############################################################################################
############# Part 1: Lake-wide species presence-absence ##########################
###############################################################################################

# Overlap of species presence in all 6 datasets
cols <- c(brewer.pal(8, "Dark2"),"#386CB0","black")
my.cols <- cols[c(10,4,9,1,6,3)]
# pdf("Figures/sp_lists_overlap.pdf", width=10, height=7) 
sp_lists_overlap <- upset(all_methods_presence, nsets = 6, sets.bar.color = my.cols, order.by = "freq",
                          mainbar.y.label = "Number of species", text.scale=c(1.7, 1.5, 1.2,1.5,1.5,1.5),
                          sets.x.label = "Number of species per dataset")
sp_lists_overlap
# dev.off()

# Get species richness totals for table
names(all_methods_presence)
length(all_methods_presence[,5][all_methods_presence[,3]>0]) # historical richness
length(all_methods_presence[,5][all_methods_presence[,4]>0]) # eDNA richness
length(all_methods_presence[,5][all_methods_presence[,5]>0]) # electrofishing richness
length(all_methods_presence[,5][all_methods_presence[,6]>0]) # fyke netting richness
length(all_methods_presence[,5][all_methods_presence[,7]>0]) # gillnetting richness
length(all_methods_presence[,5][all_methods_presence[,8]>0]) # seining richness
length(rowSums(all_methods_presence[,5:8])[rowSums(all_methods_presence[,5:8])]>0) # all capture gears richness

### Species richness per habitat type for capture gears
y <- rowSums(all_methods_presence[,c(5,6,8)]) # subset to nearshore capture gears (electrofishing, fyke netting, seining)
length(y[y>0]) # species richness across all nearshore capture gears
# pelagic species richness for capture gears is equivalent to gillnetting (the only pelagic gear)

### Species richness per habitat type for eDNA
eDNA_metadata <- read.csv("/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/oneida_eDNA_sample_metadata.csv", header = TRUE)
eDNA_dat_x <- eDNA_dat[, -grep(".sp.", colnames(eDNA_dat))][!grepl("Y", eDNA_dat$Site),] # remove subsamples (labeled with "Y")
eDNA_dat_x$Habitat <- c(eDNA_metadata$Location.type[-12]) # remove field blank
eDNA_dat_x$Habitat <- gsub("?","",eDNA_dat_x$Habitat, fixed = TRUE)
eDNA_dat_x$Habitat <- gsub("Inlet","Nearshore",eDNA_dat_x$Habitat, fixed = TRUE)
sp_rich_habitat <- rowsum(eDNA_dat_x[,2:(ncol(eDNA_dat_x)-1)], eDNA_dat_x$Habitat)
rowSums(sp_rich_habitat != 0) # eDNA species richness per habitat type

### Total estimated species richness per method
specpool(eDNA_dat[, -grep(".sp.", colnames(eDNA_dat))][,-1])
specpool(ef_dat[, -grep(".sp.", colnames(ef_dat))][,-1])
specpool(fyke_dat[, -grep(".sp.", colnames(fyke_dat))][,-1])
specpool(gillnet_dat[,-1])
specpool(seine_dat[, -grep(".sp.", colnames(seine_dat))][,-1])

###############################################################################################
################ Part 2: Species richness analyses #############################
###############################################################################################
library(vegan)
library(scales)
library(pals)
library(ggplot2)

### Species accumulation curves
### NOTE: this is now replaced this with code in Oneida_optimal_gear_combos.R
eDNA_dat_spp <- eDNA_dat[, -grep(".sp.", colnames(eDNA_dat))]
sp_acc <- specaccum(eDNA_dat_spp[,-1])
# pdf("/Users/kbja10/Documents/Cornell/Research/Oneida/Figures/spp_accum_curve_2.12.pdf", width=8, height=6)
plot(sp_acc, ci.type="poly", col=my.cols[2], lwd=3, ci.lty=0, ci.col=alpha(my.cols[2], 0.2), 
     ylab="Number of species",xlab="Number of sites", cex.axis=1.5, cex.lab=1.5, ylim = c(0,60))
ef_dat_spp <- ef_dat[, -grep(".sp.", colnames(ef_dat))]
sp_acc_ef <- specaccum(ef_dat_spp[,-1])
plot(sp_acc_ef, ci.type="poly", col=my.cols[3], lwd=3, ci.lty=0, ci.col=alpha(my.cols[3], 0.3), add=T)
fyke_dat_spp <- fyke_dat[, -grep(".sp.", colnames(fyke_dat))]
sp_acc_fyke <- specaccum(fyke_dat_spp[,-1])
plot(sp_acc_fyke, ci.type="poly", col=my.cols[5], lwd=3, ci.lty=0, ci.col=alpha(my.cols[5], 0.3), add=T)
gillnet_dat_spp <- gillnet_dat # no .sp. in gillnet dataset
sp_acc_gillnet <- specaccum(gillnet_dat_spp[,-1])
plot(sp_acc_gillnet, ci.type="poly", col=my.cols[6], lwd=3, ci.lty=0, ci.col=alpha(my.cols[6], 0.3), add=T)
seine_dat_spp <- seine_dat[, -grep(".sp.", colnames(seine_dat))]
sp_acc_seine <- specaccum(seine_dat_spp[,-1])
plot(sp_acc_seine, ci.type="poly", col=my.cols[4], lwd=3, ci.lty=0, ci.col=alpha(my.cols[4], 0.5), add=T)
legend(0, 60, legend=c("eDNA","Electrofishing","Fyke netting","Gillnetting","Seining"), col=my.cols[c(2,3,5,6,4)], lwd=3,bty="n")
# dev.off()

### Comparison of species occurrence for eDNA and capture gears for each habitat type
# Nearshore occurrence: eDNA vs. electrofishing, seining, and fyke netting
nearshore_occur <- bind_rows(eDNA_dat_x[eDNA_dat_x$Habitat=="Nearshore",-(ncol(eDNA_dat_x))], ef_dat_spp, seine_dat_spp, fyke_dat_spp) # nearshore datasets
nearshore_occur[is.na(nearshore_occur)] <- 0 # indicate absence w/ 0
nearshore_occur <- data.frame(capture=colSums(nearshore_occur[23:57,-1] != 0)/nrow(nearshore_occur[23:57,-1]),
                            edna=colSums(nearshore_occur[1:22,-1] != 0)/nrow(nearshore_occur[1:22,-1]))
nearshore_occur <- nearshore_occur[rowSums(nearshore_occur>0)!= 0,]
nearshore_occur$habitat <- rep("Nearshore",nrow(nearshore_occur))

# Pelagic occurrence: eDNA vs. gillnet
pelagic_occur <- bind_rows(eDNA_dat_x[eDNA_dat_x$Habitat=="Midlake",-(ncol(eDNA_dat_x))], gillnet_dat_spp) # nearshore datasets
pelagic_occur[is.na(pelagic_occur)] <- 0 # indicate absence w/ 0
pelagic_occur <- data.frame(capture=colSums(pelagic_occur[4:18,-1] != 0)/nrow(pelagic_occur[4:18,-1]),
                                 edna=colSums(pelagic_occur[1:3,-1] != 0)/nrow(pelagic_occur[1:3,-1]))
pelagic_occur <- pelagic_occur[rowSums(pelagic_occur>0)!= 0,]
pelagic_occur$habitat <- rep("Pelagic",nrow(pelagic_occur))

# Bind into single dataset
sp_occurrence <- rbind(pelagic_occur, nearshore_occur)

# Kendall's rank correlation of species occurrence using eDNA vs. traditional gears per habitat type
cor.test(sp_occurrence[sp_occurrence$habitat=="Pelagic",]$capture,sp_occurrence[sp_occurrence$habitat=="Pelagic",]$edna, method="kendall")
cor.test(sp_occurrence[sp_occurrence$habitat=="Nearshore",]$capture,sp_occurrence[sp_occurrence$habitat=="Nearshore",]$edna, method="kendall")

sp_occurrence_plot <- ggplot(sp_occurrence, aes(x=capture, y=edna, color=habitat)) +
  geom_jitter(size=4, width=0.02, height=0) + # jitter to see overlapping points
  xlab("Capture gears species occurrence") + 
  ylab(expression(atop("eDNA metabarcoding", paste("species occurrence")))) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  facet_wrap(~habitat) +
  theme(legend.position='none') +
  scale_color_manual(name="Habitat type",values=as.vector(stepped(14)[c(9,14)]))
sp_occurrence_plot
# ggsave("Figures/sp_occurrence.pdf", plot = sp_occurrence_plot, dpi=600, width=9, height=4,units="in")

###############################################################################################
################ Part 3: Species relative proportion analyses #################################
###############################################################################################
library(tidyr)
all_methods_prop <- data.frame(species=all_methods_abundance$species,
                               lapply(all_methods_abundance[,3:7], function(x) x / sum(x)))
all_methods_prop <- as.data.frame(all_methods_prop %>%
                                    pivot_longer(!species, names_to = "gear", values_to = "prop"))

# Subset to the top 8 species per dataset
top_spp <- all_methods_prop %>% 
  arrange(desc(prop)) %>% 
  group_by(gear) %>% slice(1:10)

# Define function to replace common names with scientific names
# Read in species names (used as a common-to-scientific name look-up table)
spp_name_table <- all_methods_presence[,c("Scientific.name","Common.name")]

replace_names <- function(i) {
  com_name <- gsub("[.]", " ", i)
  com_name <- gsub(" sp ", " sp.", com_name)
  if(com_name %in% spp_name_table$Common.name) {
    sci_name <- spp_name_table$Scientific.name[match(com_name, spp_name_table$Common.name)]
  } else {
    sci_name <- com_name
  }
  return(sci_name)
}
unique(top_spp$species)

# Stacked barplot of relative proportion (counts or reads) for top 10 species per dataset
top_spp$species <- replace_names(top_spp$species)
species_proportions <- ggplot(top_spp, aes(fill=species, y=prop, x=gear)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=as.vector(stepped(22))) +
  xlab("Sampling gear") + ylab("Species proportion") +
  labs(fill = "Species") +
  scale_x_discrete(labels=c("eDNA","Electrofishing","Fyke netting","Gill netting","Seining")) +
  theme_bw()
species_proportions

# ggsave("Figures/all_methods_species_proportion.pdf", plot = species_proportions, dpi=600)

###############################################################################################
################# Part 4: Habitat comparison: eDNA vs. traditional methods ###################
###############################################################################################
library(ggpubr)

# eDNA NMDS: presence-absence data
eDNA_metadata <- read.csv("/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/oneida_eDNA_sample_metadata.csv", header = TRUE)
nmds_presence_dat <- eDNA_dat_spp
nmds_presence_dat[,-1][nmds_presence_dat[,-1]>0] <- 1 # indicate presence w/ 1
NMDS_dat <- metaMDS(nmds_presence_dat[,-1],distance="bray",trymax = 200) # nmds w/ Sorensen index
plot(NMDS_dat,display = c("sites", "species"))
ordiplot(NMDS_dat,type="n")
orditorp(NMDS_dat,display="species",col="red",air=0.01)
orditorp(NMDS_dat,display="sites",cex=1.25,air=0.01)
data.scores <- as.data.frame(scores(NMDS_dat)) # turn into data frame for plotting
data.scores$Habitat <- eDNA_metadata$Location.type[-12] # all ef surveys are nearshore
data.scores$Habitat <- gsub("?","",data.scores$Habitat, fixed = TRUE)
data.scores$Habitat <- gsub("Inlet","Nearshore",data.scores$Habitat, fixed = TRUE)

# add hulls
grp.a <- data.scores[data.scores$Habitat=="Midlake",][chull(data.scores[data.scores$Habitat=="Midlake", c("NMDS1", "NMDS2")]), ]
grp.b <- data.scores[data.scores$Habitat=="Nearshore",][chull(data.scores[data.scores$Habitat=="Nearshore", c("NMDS1", "NMDS2")]), ]
hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b
hull.data

nmds_edna <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size=5, stroke=2, aes(colour=Habitat))+ 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Habitat,group=Habitat),alpha=0.30) + # add the convex hulls
  theme(axis.text = element_text(colour = "black", size = 12), 
        axis.title = element_text(size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  theme(legend.position="top") +
  theme(legend.title = element_blank()) +
  ggtitle("eDNA metabarcoding") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name="Habitat type", values=as.vector(stepped(14)[c(14,9)]), labels=c("Pelagic","Nearshore")) +
  scale_fill_manual(name="Habitat type", values=as.vector(stepped(14)[c(14,9)]), labels=c("Pelagic","Nearshore")) +
  labs(x = "NMDS1", colour = "Site", y = "NMDS2")
nmds_edna

# traditional gears NMDS: presence-absence data
nmds_presence_dat <- bind_rows(gillnet_dat_spp,ef_dat_spp,fyke_dat_spp,seine_dat_spp) # combine traditional datasets
nmds_presence_dat[is.na(nmds_presence_dat)] <- 0 # indicate absence w/ 0
nmds_presence_dat[,-1][nmds_presence_dat[,-1]>0] <- 1 # indicate presence w/ 1
nmds_presence_dat <- nmds_presence_dat[,c(TRUE, colSums(nmds_presence_dat[,-1])>0)]
nmds_presence_dat$Sampling_gear <- c(rep("Gillnetting",nrow(gillnet_dat_spp)),rep("Electrofishing",nrow(ef_dat_spp)),
                            rep("Fyke netting",nrow(fyke_dat_spp)),rep("Seining",nrow(seine_dat_spp)))
nmds_presence_dat <- nmds_presence_dat %>% select(Sampling_gear, everything()) # reorder columnns
nmds_presence_dat <- nmds_presence_dat[c(rowSums(nmds_presence_dat[,-(1:2)])>0),] # remove samples with no species detected
NMDS_dat <- metaMDS(nmds_presence_dat[,-(1:2)],distance="bray",trymax = 200) # nmds w/ Sorensen index
data.scores <- as.data.frame(scores(NMDS_dat)) # turn into data frame for plotting
data.scores$Sampling_gear <- nmds_presence_dat$Sampling_gear
data.scores$Sampling_gear <- factor(data.scores$Sampling_gear, levels=unique(data.scores$Sampling_gear))

# add hulls
grp.a <- data.scores[data.scores$Sampling_gear=="Gillnetting",][chull(data.scores[data.scores$Sampling_gear=="Gillnetting", c("NMDS1", "NMDS2")]), ]
grp.b <- data.scores[data.scores$Sampling_gear=="Electrofishing",][chull(data.scores[data.scores$Sampling_gear=="Electrofishing", c("NMDS1", "NMDS2")]), ]
grp.c <- data.scores[data.scores$Sampling_gear=="Fyke netting",][chull(data.scores[data.scores$Sampling_gear=="Fyke netting", c("NMDS1", "NMDS2")]), ]
grp.d <- data.scores[data.scores$Sampling_gear=="Seining",][chull(data.scores[data.scores$Sampling_gear=="Seining", c("NMDS1", "NMDS2")]), ]
hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)
hull.data

nmds_traditional <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size=5, stroke=2, aes(colour=Sampling_gear, shape=Sampling_gear))+ 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Sampling_gear,group=Sampling_gear),alpha=0.30) + # add the convex hulls
  theme(axis.text = element_text(colour = "black", size = 12), 
        axis.title = element_text(size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  theme(legend.position="top") +
  theme(legend.title = element_blank()) +
  ggtitle("Capture sampling gears") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name="Sampling gear", values=as.vector(stepped(16)[c(14,9,11,12)])) +
  scale_fill_manual(name="Sampling gear", values=as.vector(stepped(16)[c(14,11,11,10)])) +
  scale_shape_manual(name="Sampling gear",values=c(15:18)) +
  labs(x = "NMDS1", colour = "Sampling gear", y = "NMDS2")
nmds_traditional

nmds_plot <- ggarrange(nmds_edna, nmds_traditional, labels = c("(A)", "(B)"),ncol = 2, nrow = 1)
nmds_plot

# ggsave("Figures/nmds_plot.pdf", plot = nmds_plot, dpi=600, width=12, height=7,units="in")

###############################################################################################
################# Part 5: Site comparison -- eDNA vs. electrofishing ########################
###############################################################################################
site_comparison <- bind_rows(eDNA_dat_spp, ef_dat_spp) # combine eDNA and ef datasets
site_comparison[is.na(site_comparison)] <- 0 # indicate absence w/ 0
site_comparison <- site_comparison[-c(1,12,15),] # remove pelagic eDNA samples
eDNA_metadata_site <- eDNA_metadata$Site[-c(1,12,13,16)] 
site_comparison$Site_name <- c(eDNA_metadata_site, site_comparison$Site[23:30]) # all ef surveys are nearshore
site_comparison$Site_code <- c(NA, "South Bay Sailboat Marina run",NA,"South Bay Sailboat Marina run",NA,NA,
                               "Jewell run","Jewell run","Cleveland run","Cleveland run","Constantia run","Constantia run",
                               "Poddygut Bay run","Poddygut Bay run",NA,"Aero Marina run","Aero Marina run", NA, 
                               "Fisher Bay run","Fisher Bay run","Shackelton Point run","Shackelton Point run",site_comparison$Site[23:30])
site_comparison <- site_comparison %>% select(c(Site,Site_name,Site_code), everything())
site_comparison  <- site_comparison[!(is.na(site_comparison$Site_code)),]

# Get eDNA species richness per site
site_comparison_edna <- rowsum(site_comparison[1:16,-c(1:3)], site_comparison$Site_code[1:16])
site_comparison_edna$Richness <- rowSums(site_comparison_edna != 0) # species richness per site

# Get electrofishing species richness per site
site_comparison_ef <- site_comparison[17:24,-c(1:3)]
rownames(site_comparison_ef) <- site_comparison$Site_code[17:24]
site_comparison_ef$Richness <- rowSums(site_comparison_ef != 0) # species richness per site

# differences/correlations between richness per site
site_comparison_edna$Richness-site_comparison_ef$Richness
mean(site_comparison_edna$Richness) # 19.25
sd(site_comparison_edna$Richness) # 4.80
mean(site_comparison_ef$Richness) # 14.75
sd(site_comparison_ef$Richness) # 3.84
wilcox.test(site_comparison_edna$Richness,site_comparison_ef$Richness)
cor.test(site_comparison_edna$Richness,site_comparison_ef$Richness)
plot(site_comparison_edna$Richness,site_comparison_ef$Richness)

site_comparison <- data.frame(Transect=rep(rownames(site_comparison_edna),2),
                              Sample=rep(c("eDNA","Electrofishing"), each=8),
                              Richness=c(site_comparison_edna$Richness,site_comparison_ef$Richness))

site_richness <- ggplot(site_comparison, aes(fill=Sample, y=Richness, x=Transect)) + 
  geom_bar(position="dodge", stat="identity", width = 0.65) +
  xlab("Transect") + ylab("Species richness") +
  theme_bw() +
  theme(text = element_text(size=15)) +
  scale_fill_manual(name="Sampling method",labels=c("eDNA","Electrofishing"),values=my.cols[c(2:3)])
site_richness
# ggsave("Figures/site_richness.pdf", plot = site_richness, dpi=600)

# Correlation in relative abundance per site: eDNA vs. ef
edna_props <- as.data.frame(lapply(site_comparison_edna, function(x) x / sum(x)))
edna_props[is.na(edna_props)] <- 0

ef_props <- as.data.frame(lapply(site_comparison_ef, function(x) x / sum(x)))
ef_props[is.na(ef_props)] <- 0
c.coeffs <- sapply(1:(ncol(edna_props)-1), FUN = function(x) cor(edna_props[,x], ef_props[,x], method="kendall"))
p.vals <- sapply(1:(ncol(edna_props)-1), FUN = function(x) cor.test(edna_props[,x], ef_props[,x], method="kendall")$p.value)
df <- data.frame(Species = colnames(edna_props[,-(ncol(edna_props))]), Correlation = c.coeffs, p.vals=p.vals)
df <- df[order(df$Correlation, decreasing = TRUE),]
df <- df[!is.na(df$Correlation),]
nrow(df)
df$signif <- ifelse(((df$p.vals<=0.05)),1,0)
summary(df$Correlation, na.rm = TRUE)

sp_correlation <- ggplot(df, aes(x=Correlation)) + geom_histogram(bins=10,binwidth = 0.2) +
  geom_vline(data=df, aes(xintercept=mean(Correlation)),linetype="dashed") +
  xlim(-1,1.1) + ylab("Count") + scale_y_continuous(breaks=c(2,4,6,8,10)) +
  theme_bw()
sp_correlation

df$Species <- replace_names(df$Species)
sp_correlation <- ggplot(df, aes(x=reorder(Species,Correlation), y=Correlation, group=signif)) +
  geom_point(aes(shape=as.factor(signif), size = 1.3)) +
  scale_shape_manual(values=c(16,8)) +
  coord_flip() + labs(y="Rank correlation (Kendall's tau)",x="") +
  theme_bw() +
  theme(legend.position = "none")
sp_correlation
# ggsave("Figures/site_corr_coeffs.pdf", plot = sp_correlation, dpi=600)
