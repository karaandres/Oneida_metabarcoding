### Oneida eDNA metabarcoding analysis
### This code compares Oneida Lake species inventories from traditional (EF, seine, fyke, and gillnet) survey data, 
### historical records, and eDNA metabarcoding datasets 
### Created by K. Andres and T. Lambert

historical_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/oneida_historical_species_list.csv", header = TRUE)
eDNA_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/eDNA_dat.csv", header = TRUE)
ef_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/ef_dat.csv", header = TRUE)
fyke_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/fyke_dat.csv", header = TRUE)
gillnet_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/gillnet_dat.csv", header = TRUE)
seine_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/seine_dat.csv", header = TRUE)
all_methods_abundance <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/all_methods_abundance.csv", header = TRUE)
all_methods_presence <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/all_methods_presence.csv", header = TRUE)

###############################################################################################
############# Part 1: Lake-wide species presence-absence ##########################
###############################################################################################

# Venn diagram of species presence in all 6 datasets
library(RColorBrewer)
library(stringr)
library(UpSetR)
library(ggplot2)

cols <- c(brewer.pal(8, "Dark2"),"#386CB0","black")
my.cols <- cols[c(10,4,9,1,6,3)]
# pdf("/Users/kbja10/Documents/Cornell/Research/Oneida/Figures/sp_lists_overlap.pdf", width=10, height=7) 
sp_lists_overlap <- upset(all_methods_presence, nsets = 6, sets.bar.color = my.cols, order.by = "freq",
                          mainbar.y.label = "Number of overlapping species", text.scale=c(1.7, 1.5, 1.2,1.5,1.5,1.5),
                          sets.x.label = "Number of species per dataset")
sp_lists_overlap
# dev.off()

# Get species richness totals for table
names(all_methods_presence)
length(all_methods_presence[,5][all_methods_presence[,3]>0]) # hist richness
length(all_methods_presence[,5][all_methods_presence[,4]>0]) # eDNA richness
length(all_methods_presence[,5][all_methods_presence[,5]>0]) # ef richness
length(all_methods_presence[,5][all_methods_presence[,6]>0]) # fyke richness
length(all_methods_presence[,5][all_methods_presence[,7]>0]) # gill richness
length(all_methods_presence[,5][all_methods_presence[,8]>0]) # seine richness
length(rowSums(all_methods_presence[,5:8])[rowSums(all_methods_presence[,5:8])]>0) # all traditional gears richness

z <- all_methods_presence$Historical
length(z[z>0]) # species richness of historical dataset
y <- rowSums(all_methods_presence[,c(5,6,8)]) # nearshore gears
length(y[y>0]) # species richness across all nearshore gears

###############################################################################################
################ Part 2: Species richness analyses #############################
###############################################################################################
library(vegan)
library(scales)
library(pals)
library(ggplot2)

### Species accumulation curves
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

legend(0, 60, legend=c("eDNA","Electrofishing","Fyke net","Gill net","Seining"), col=my.cols[c(2,3,5,6,4)], lwd=3,bty="n")
# dev.off()

### Total species richness per method
specpool(eDNA_dat_spp[,-1])
specpool(ef_dat_spp[,-1])
specpool(fyke_dat_spp[,-1])
specpool(seine_dat_spp[,-1])
specpool(gillnet_dat_spp[,-1])

### Species richness per habitat type for eDNA
eDNA_metadata <- read.csv("/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/oneida_eDNA_sample_metadata.csv", header = TRUE)
eDNA_dat_x <- eDNA_dat_spp[!grepl("Y", eDNA_dat_spp$Site),] # remove subsamples (labeled with "Y")
eDNA_dat_x$Habitat <- c(eDNA_metadata$Location.type[-12]) # all ef surveys are nearshore
eDNA_dat_x$Habitat <- gsub("?","",eDNA_dat_x$Habitat, fixed = TRUE)
sp_rich_habitat <- rowsum(eDNA_dat_x[,2:(ncol(eDNA_dat_x)-1)], eDNA_dat_x$Habitat)
rowSums(sp_rich_habitat != 0) # species richness per habitat type

### Comparison of species occurrence for eDNA and traditional gears for each habitat type
# Nearshore occurrence: eDNA vs. ef, seine, fyke
nearshore_occur <- bind_rows(eDNA_dat_x[eDNA_dat_x$Habitat=="Nearshore",-(ncol(eDNA_dat_x))], ef_dat_spp, seine_dat_spp, fyke_dat_spp) # nearshore datasets
nearshore_occur[is.na(nearshore_occur)] <- 0 # indicate absence w/ 0
nearshore_occur <- data.frame(traditional=colSums(nearshore_occur[19:87,-1] != 0)/nrow(nearshore_occur[19:87,-1]),
                            edna=colSums(nearshore_occur[1:18,-1] != 0)/nrow(nearshore_occur[1:18,-1]))
nearshore_occur <- nearshore_occur[rowSums(nearshore_occur>0)!= 0,]
nearshore_occur$habitat <- rep("Nearshore",nrow(nearshore_occur))

# Pelagic occurrence: eDNA vs. gillnet
pelagic_occur <- bind_rows(eDNA_dat_x[eDNA_dat_x$Habitat=="Midlake",-(ncol(eDNA_dat_x))], gillnet_dat_spp) # nearshore datasets
pelagic_occur[is.na(pelagic_occur)] <- 0 # indicate absence w/ 0
pelagic_occur <- data.frame(traditional=colSums(pelagic_occur[4:18,-1] != 0)/nrow(pelagic_occur[4:18,-1]),
                                 edna=colSums(pelagic_occur[1:3,-1] != 0)/3)
pelagic_occur <- pelagic_occur[rowSums(pelagic_occur>0)!= 0,]
pelagic_occur$habitat <- rep("Pelagic",nrow(pelagic_occur))

# Bind into single dataset
sp_occurrence <- rbind(pelagic_occur, nearshore_occur)

# Kendall's rank correlation of species occurrence using eDNA vs. traditional gears per habitat type
cor.test(sp_occurrence[sp_occurrence$habitat=="Pelagic",]$traditional,sp_occurrence[sp_occurrence$habitat=="Pelagic",]$edna, method="kendall")
cor.test(sp_occurrence[sp_occurrence$habitat=="Nearshore",]$traditional,sp_occurrence[sp_occurrence$habitat=="Nearshore",]$edna, method="kendall")

sp_occurrence_plot <- ggplot(sp_occurrence, aes(x=traditional, y=edna, color=habitat)) +
  geom_jitter(size=4, width=0.02, height=0) + # jitter to see overlapping points
  xlab("Traditional methods species occurrence") + 
  ylab(expression(atop("eDNA metabarcoding", paste("species occurrence")))) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  facet_wrap(~habitat) +
  theme(legend.position='none') +
  scale_color_manual(name="Habitat type",values=as.vector(stepped(14)[c(9,14)]))
sp_occurrence_plot

# ggsave("/Users/kbja10/Documents/Cornell/Research/Oneida/Figures/sp_occurrence.pdf", plot = sp_occurrence_plot, dpi=600, width=9, height=4,units="in")
# ggsave("/Users/kbja10/Github/Oneida_metabarcoding/markdown_images/sp_occurrence.png", plot = sp_occurrence_plot, dpi=600, width=9, height=4,units="in")

###############################################################################################
################ Part 3: Species relative proportion analyses #################################
###############################################################################################

all_methods_prop <- data.frame(species=all_methods_abundance$species,
                               lapply(all_methods_abundance[,3:7], function(x) x / sum(x)))
all_methods_prop <- as.data.frame(all_methods_prop %>%
                                    pivot_longer(!species, names_to = "gear", values_to = "prop"))

# Subset to the top 8 species per dataset
top_spp <- all_methods_prop %>% 
  arrange(desc(prop)) %>% 
  group_by(gear) %>% slice(1:8)

# Stacked barplot of relative proportion (counts or reads) for top 8 species per dataset
species_proportions <- ggplot(top_spp, aes(fill=species, y=prop, x=gear)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=as.vector(stepped(21))) +
  xlab("Sampling gear") + ylab("Species proportion") +
  labs(fill = "Species") +
  scale_x_discrete(labels=c("eDNA","Electrofishing","Fyke net","Gill net","Seine")) +
  theme_bw()
species_proportions

# ggsave("/Users/kbja10/Documents/Cornell/Research/Oneida/Figures/all_methods_species_proportion.pdf", plot = species_proportions, dpi=600)
# ggsave("/Users/kbja10/Github/Oneida_metabarcoding/markdown_images/all_methods_species_proportion_xcarp.png", plot = species_proportions, dpi=600)

###############################################################################################
################# Part 4: Habitat comparison: eDNA vs. traditional methods ###################
###############################################################################################
library(ggpubr)

# eDNA NMDS: presence-absence data
eDNA_metadata <- read.csv("/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/oneida_eDNA_sample_metadata.csv", header = TRUE)
nmds_presence_dat <- eDNA_dat_spp
nmds_presence_dat <- nmds_presence_dat[!grepl("Y", nmds_presence_dat$Site),] # remove subsamples (labeled with "Y")
nmds_presence_dat[,-1][nmds_presence_dat[,-1]>0] <- 1 # indicate presence w/ 1
NMDS_dat <- metaMDS(nmds_presence_dat[,-1],distance="bray",trymax = 200) # nmds w/ Sorensen index
data.scores <- as.data.frame(scores(NMDS_dat)) # turn into data frame for plotting
data.scores$Habitat <- eDNA_metadata$Location.type[-12] # all ef surveys are nearshore
data.scores$Habitat <- gsub("?","",data.scores$Habitat, fixed = TRUE)

# add hulls
grp.a <- data.scores[data.scores$Habitat=="Midlake",][chull(data.scores[data.scores$Habitat=="Midlake", c("NMDS1", "NMDS2")]), ]
grp.b <- data.scores[data.scores$Habitat=="Nearshore",][chull(data.scores[data.scores$Habitat=="Nearshore", c("NMDS1", "NMDS2")]), ]
grp.c <- data.scores[data.scores$Habitat=="Inlet",][chull(data.scores[data.scores$Habitat=="Inlet", c("NMDS1", "NMDS2")]), ]
hull.data <- rbind(grp.a, grp.b, grp.c)  #combine grp.a and grp.b
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
  ggtitle("eDNA metabarcoding") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name="Habitat type", values=as.vector(stepped(14)[c(2,14,9)]), labels=c("Tributary","Pelagic","Nearshore")) +
  scale_fill_manual(name="Habitat type", values=as.vector(stepped(14)[c(2,14,9)]), labels=c("Tributary","Pelagic","Nearshore")) +
  labs(x = "NMDS1", colour = "Site", y = "NMDS2")
nmds_edna

# traditional gears NMDS: presence-absence data
nmds_presence_dat <- bind_rows(gillnet_dat_spp,ef_dat_spp,fyke_dat_spp,seine_dat_spp) # combine traditional datasets
nmds_presence_dat[is.na(nmds_presence_dat)] <- 0 # indicate absence w/ 0
nmds_presence_dat[,-1][nmds_presence_dat[,-1]>0] <- 1 # indicate presence w/ 1
nmds_presence_dat <- nmds_presence_dat[,c(TRUE, colSums(nmds_presence_dat[,-1])>0)]
nmds_presence_dat$Sampling_gear <- c(rep("Gill net",nrow(gillnet_dat_spp)),rep("Electrofishing",nrow(ef_dat_spp)),
                            rep("Fyke net",nrow(fyke_dat_spp)),rep("Seine net",nrow(seine_dat_spp)))
nmds_presence_dat <- nmds_presence_dat %>% select(Sampling_gear, everything()) # reorder columnns
nmds_presence_dat <- nmds_presence_dat[c(rowSums(nmds_presence_dat[,-(1:2)])>0),] # remove samples with no species detected
NMDS_dat <- metaMDS(nmds_presence_dat[,-(1:2)],distance="bray",trymax = 200) # nmds w/ Sorensen index
data.scores <- as.data.frame(scores(NMDS_dat)) # turn into data frame for plotting
data.scores$Sampling_gear <- nmds_presence_dat$Sampling_gear
data.scores$Sampling_gear <- factor(data.scores$Sampling_gear, levels=unique(data.scores$Sampling_gear))

# add hulls
grp.a <- data.scores[data.scores$Sampling_gear=="Gill net",][chull(data.scores[data.scores$Sampling_gear=="Gill net", c("NMDS1", "NMDS2")]), ]
grp.b <- data.scores[data.scores$Sampling_gear=="Electrofishing",][chull(data.scores[data.scores$Sampling_gear=="Electrofishing", c("NMDS1", "NMDS2")]), ]
grp.c <- data.scores[data.scores$Sampling_gear=="Fyke net",][chull(data.scores[data.scores$Sampling_gear=="Fyke net", c("NMDS1", "NMDS2")]), ]
grp.d <- data.scores[data.scores$Sampling_gear=="Seine net",][chull(data.scores[data.scores$Sampling_gear=="Seine net", c("NMDS1", "NMDS2")]), ]
hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
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
  ggtitle("Traditional sampling gears") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name="Sampling gear", values=as.vector(stepped(16)[c(14,9,11,12)])) +
  scale_fill_manual(name="Sampling gear", values=as.vector(stepped(16)[c(14,11,11,10)])) +
  scale_shape_manual(name="Sampling gear",values=c(15:18)) +
  labs(x = "NMDS1", colour = "Sampling gear", y = "NMDS2")
nmds_traditional

nmds_plot <- ggarrange(nmds_edna, nmds_traditional, labels = c("(A)", "(B)"),ncol = 2, nrow = 1)
nmds_plot

# ggsave("/Users/kbja10/Documents/Cornell/Research/Oneida/Figures/nmds_plot.pdf", plot = nmds_plot, dpi=600, width=12, height=7,units="in")
# ggsave("/Users/kbja10/Github/Oneida_metabarcoding/markdown_images/nmds_plot.png", plot = nmds_plot, dpi=600, width=12, height=7,units="in")

###############################################################################################
################# Part 5: Site comparison -- eDNA vs. electrofishing ########################
###############################################################################################
site_comparison <- bind_rows(eDNA_dat_spp, ef_dat_spp) # combine eDNA and ef datasets
site_comparison <- site_comparison[!grepl("Y", site_comparison$Site),] # remove subsamples (labeled with "Y")
site_comparison[is.na(site_comparison)] <- 0 # indicate absence w/ 0
rownames(site_comparison) <- c(paste0("eDNA_",c(1:11,13:26)),paste0("ef_",c(1:24)))
site_comparison$Habitat <- c(eDNA_metadata$Location.type[-12], rep("Nearshore",24)) # all ef surveys are nearshore
site_comparison$Habitat <- gsub("?","",site_comparison$Habitat, fixed = TRUE)
site_comparison$Site <- c(eDNA_metadata$Site[-12], site_comparison$Site[26:49]) # all ef surveys are nearshore
site_comparison <- site_comparison[site_comparison$Habitat=="Nearshore",]
site_comparison$Site
site_comparison$Site_code <- c(NA, "South","South",NA,"Jewell","Jewell","Cleveland","Cleveland",
                               "Constantia","Constantia","Poddygut","Poddygut","Aero","Aero",
                               "Fisher","Fisher","Shackelton","Shackelton",rep("Constantia",3),rep("Fisher",3),
                               rep("Cleveland",3),rep("Poddygut",3),rep("Aero",3),rep("Shackelton",3),rep("South",3),rep("Jewell",3))
site_comparison <- site_comparison %>% select(c(Habitat,Site,Site_code), everything())
site_comparison  <- site_comparison[!(is.na(site_comparison$Site_code)),]

# Get eDNA totals per site
site_comparison_edna <- rowsum(site_comparison[1:16,-c(1:3)], site_comparison$Site_code[1:16])
site_comparison_edna$Richness <- rowSums(site_comparison_edna != 0) # species richness per site

# Get ef totals per site
site_comparison_ef <- rowsum(site_comparison[17:40,-c(1:3)], site_comparison$Site_code[17:40])
site_comparison_ef$Richness <- rowSums(site_comparison_ef != 0) # species richness per site

# differences/correlations between richness per site
site_comparison_edna$Richness-site_comparison_ef$Richness
mean(site_comparison_edna$Richness) # 19.25
sd(site_comparison_edna$Richness) # 4.80
mean(site_comparison_ef$Richness) # 17.125
sd(site_comparison_ef$Richness) # 3.482
wilcox.test(site_comparison_edna$Richness,site_comparison_ef$Richness)
cor.test(site_comparison_edna$Richness,site_comparison_ef$Richness)
plot(site_comparison_edna$Richness,site_comparison_ef$Richness)

site_comparison <- data.frame(Transect=rep(rownames(site_comparison_edna),2),
                              Sample=rep(c("eDNA","Electrofishing"), each=8),
                              Richness=c(site_comparison_edna$Richness,site_comparison_ef$Richness))

site_richness <- ggplot(site_comparison, aes(fill=Sample, y=Richness, x=Transect)) + 
  geom_bar(position="stack", stat="identity") +
  xlab("Transect") + ylab("Species richness") +
  theme_bw() +
  theme(text = element_text(size=15)) +
  scale_fill_manual(name="Sampling method",labels=c("eDNA","Electrofishing"),values=my.cols[c(2:3)])
site_richness
# ggsave("/Users/kbja10/Documents/Cornell/Research/Oneida/Figures/site_richness.pdf", plot = site_richness, dpi=600)
# ggsave("/Users/kbja10/Github/Oneida_metabarcoding/markdown_images/site_richness.png", plot = site_richness, dpi=600)

# Correlation in relative abundance per site: eDNA vs. ef
edna_props <- as.data.frame(lapply(site_comparison_edna, function(x) x / sum(x)))
edna_props[is.na(edna_props)] <- 0

ef_props <- as.data.frame(lapply(site_comparison_ef, function(x) x / sum(x)))
ef_props[is.na(ef_props)] <- 0

c.coeffs <- sapply(1:(ncol(edna_props)-1), FUN = function(x) cor(edna_props[,x], ef_props[,x], method="kendall"))
p.vals <- sapply(1:(ncol(edna_props)-1), FUN = function(x) cor.test(edna_props[,x], ef_props[,x], method="kendall")$p.value)
df <- data.frame(Species = colnames(edna_props[,-(ncol(edna_props))]), Correlation = c.coeffs, p.vals=p.vals)
df[order(df$Correlation, decreasing = TRUE),]
df <- df[!is.na(df$Correlation),]
nrow(df)
nrow(df[df$p.vals<0.05,])
summary(df$Correlation, na.rm = TRUE)

sp_correlation <- ggplot(df, aes(x=Correlation)) + geom_histogram(bins=10,binwidth = 0.2) +
  geom_vline(data=df, aes(xintercept=mean(Correlation)),linetype="dashed") +
  xlim(-1,1.1) + ylab("Count") + scale_y_continuous(breaks=c(2,4,6,8,10)) +
  theme_bw()
sp_correlation
# ggsave("/Users/kbja10/Github/Oneida_metabarcoding/markdown_images/sp_correlation.png", plot = sp_correlation, dpi=600, width=9, height=4,units="in")

