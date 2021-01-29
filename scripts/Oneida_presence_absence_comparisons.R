### Oneida eDNA metabarcoding analysis
### This code compares Oneida Lake species inventories from traditional (EF, seine, fyke, and gillnet) survey data, 
### historical records, and eDNA metabarcoding datasets 
### Created by K. Andres and T. Lambert

historical_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/oneida_historical_species_list.csv", header = TRUE)
sp_read_count_by_site <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/sp_read_count_by_site.csv")
ef_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/Oneida_2017_comm_dataset.csv")
fyke_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/2017 oneida fyke.csv")
gillnet_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/2017 oneida gill net.csv")
seine_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/2017 oneida seine.csv")

library(dplyr)
library(tidyr)
library(stringr)

### Reconcile inconsistencies in species names across datasets and create community datasets (vegan)

# eDNA dataset
to_replace <- c(" ","Ambloplites.rupestris","Gobius","Labidesthes","northern.hog.sucker","Percopsis.omiscomaycus","pumpkinseed.sunfish","tessallated.darter",
                "shortnose.gar","Cottidae","Cyprinidae","Cyprinus","Luxilus","Notropis","Percina","blacktip.jumprock")
replace_with <- c("\\.","rock.bass","round.goby","brook.silverside","northern.hogsucker","trout.perch","pumpkinseed","tessellated.darter",
                  "longnose.gar","Cottidae.sp.","Cyprinidae.sp.","Cyprinus.sp.","Luxilus.sp.","Notropis.sp.","Percina.sp.","Moxostoma.sp.")
names(replace_with) <- to_replace
sp_read_count_by_site$scomnames <- str_replace_all(sp_read_count_by_site$scomnames, replace_with)
eDNA_dat <- t(as.matrix(sp_read_count_by_site[,-c(1:4)])) # transpose to create community dataset
colnames(eDNA_dat) <- sp_read_count_by_site$scomnames
eDNA_dat <- data.frame(Site=rownames(eDNA_dat),eDNA_dat)
  
# Electrofishing dataset
to_replace <- c("chubsucker","gizzard.shad","killifish","lepomis","pumpkinseedXgreen.sunfish","redhorse","unknown.shiner")
replace_with <- c("creek.chubsucker","American.gizzard.shad","banded.killifish","Lepomis.1","Lepomis.2","Moxostoma.sp.","Notropis.sp.")
names(replace_with) <- to_replace
colnames(ef_dat) <- str_replace_all(colnames(ef_dat), replace_with)
ef_dat <- ef_dat %>% mutate(Lepomis.sp.=Lepomis.1+Lepomis.2, .keep = "unused")
# Collapse reads/species counts by site (triplicate eDNA samples/multiple EF passes)
ef_dat_sites <- as.data.frame(ef_dat %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), sum)))
ef_dat_sites$Site <- gsub(" [[:digit:]]+", "", oneida_ef$Site) # remove pass number 

# Fyke dataset
fyke_dat <- fyke_dat %>% 
  mutate(largemouth.bass=largemouth.bass+largemouth.bass.YOY, .keep = "unused") %>%
  mutate(smallmouth.bass=smallmouth.bass+smallmouth.bass.YOY, .keep = "unused") %>%
  mutate(walleye=walleye+walleye.YOY, .keep = "unused") %>%
  mutate(yellow.perch=yellow.perch+yellow.perch.YOY, .keep = "unused")
to_replace <- c("american.eel","Cisco","Common.carp","darter","gizzard.shad","redhorse")
replace_with <- c("American.eel","cisco","common.carp","Etheostoma.sp.","American.gizzard.shad","Moxostoma.sp.")
names(replace_with) <- to_replace
colnames(fyke_dat) <- str_replace_all(colnames(fyke_dat), replace_with)
fyke_dat$Other
fyke_dat <- fyke_dat %>% select(Site, American.eel:yellow.perch)
fyke_dat$satinfin.shiner <- c(rep(0,4),1,rep(0,31))
fyke_dat$bluntnose.minnow <- c(rep(0,7),1,0,0,0,1,0,0,0,0,0,2,rep(0,13),2,0,0,0,0)
fyke_dat$creek.chubsucker <- c(rep(0,18),1,0,0,3,rep(0,14))
fyke_dat$PSXGSF <- c(rep(0,26),1,0,1,rep(0,7))
fyke_dat <- fyke_dat %>% mutate(Lepomis.sp.=PSXGSF+Lepomis.YOY, .keep = "unused")

# Gill net dataset
to_replace <- c("gizzard.shad","redhorse.sp.","tiger.muskellunge")
replace_with <- c("American.gizzard.shad","Moxostoma.sp.","tiger.musky")
names(replace_with) <- to_replace
colnames(gillnet_dat) <- str_replace_all(colnames(gillnet_dat), replace_with)
gillnet_dat <- gillnet_dat %>% select(Site, walleye:round.goby)
gillnet_dat[is.na(gillnet_dat)] <- 0
  
# Seine dataset
to_replace <- c(" ","gizzard.shad","shiner.sp.","lepomis.sp.")
replace_with <- c("\\.","American.gizzard.shad","Notropis.sp.","Lepomis.sp.")
names(replace_with) <- to_replace
seine_dat$Species <- str_replace_all(seine_dat$Species, replace_with)
seine_dat <- as.data.frame(seine_dat %>% pivot_wider(names_from = Species, values_from = Number.Caught) %>%
  group_by(Site) %>%
  summarise(across(bluegill:chain.pickerel, sum, na.rm = TRUE)))

### Get species totals across all sites

# eDNA species totals 
eDNA_totals <- data.frame(species = colnames(eDNA_dat[,-1]),
                          number = colSums(eDNA_dat[,-1]))
eDNA_totals <- eDNA_totals[eDNA_totals$number > 0, ]
eDNA_totals <- eDNA_totals[order(eDNA_totals$species),]
eDNA_totals$gear <- rep("eDNA",nrow(eDNA_totals))

# Electrofishing species totals
ef_totals <- data.frame(species = colnames(ef_dat[,-1]),
                        number = colSums(ef_dat[,-1]))
ef_totals <- ef_totals[ef_totals$number > 0, ]
ef_totals <- ef_totals[order(ef_totals$species),]
ef_totals$gear <- rep("ef",nrow(ef_totals))

# Fyke species totals
fyke_totals <- data.frame(species = colnames(fyke_dat[,-1]),
                          number = colSums(fyke_dat[,-1]))
fyke_totals <- fyke_totals[fyke_totals$number > 0, ]
fyke_totals <- fyke_totals[order(fyke_totals$species),]
fyke_totals$gear <- rep("fyke",nrow(fyke_totals))

# Gillnet species totals
gillnet_totals <- data.frame(species = colnames(gillnet_dat[,-1]),
                             number = colSums(gillnet_dat[,-1]))
gillnet_totals <- gillnet_totals[gillnet_totals$number > 0, ]
gillnet_totals <- gillnet_totals[order(gillnet_totals$species),]
gillnet_totals$gear <- rep("gillnet",nrow(gillnet_totals))
  
# Seine species totals
seine_totals <- data.frame(species = colnames(seine_dat[,-1]),
                           number = colSums(seine_dat[,-1]))
seine_totals <- seine_totals[seine_totals$number > 0, ]
seine_totals <- seine_totals[order(seine_totals$species),]
seine_totals$gear <- rep("seine", nrow(seine_totals))

# All methods totals
historical_dat$scomnames <- gsub(" ", ".", historical_dat$scomnames)
historical_totals <- data.frame(species=historical_dat$scomnames,number=historical_dat$since_1990,
                                gear=rep("historical",nrow(historical_dat)))
all_methods_totals <- rbind(historical_totals, eDNA_totals,ef_totals,fyke_totals,
                    gillnet_totals,seine_totals)
all_methods_totals <- as.data.frame(all_methods_totals %>% 
  pivot_wider(names_from = gear, values_from = number))
all_methods_totals[is.na(all_methods_totals)] <- 0

# All methods presence-absence
all_methods_presence <- all_methods_totals[!grepl(".sp", all_methods_totals$species),] # remove genus-level IDs
all_methods_presence[,2:7][all_methods_presence[,2:7]>0] <- 1 # indicate presence w/ 1
all_methods_presence <- left_join(all_methods_presence, historical_dat[,1:2], by=c("species"="scomnames")) %>%
  select(species.y, everything()) %>%  
  arrange(desc(historical), desc(eDNA), desc(ef), desc(fyke), desc(gillnet), desc(seine))
all_methods_presence$species <- gsub("\\.", " ", all_methods_presence$species)
colnames(all_methods_presence) <- c("Scientific name","Common name","Historical","eDNA","Electrofishing",
                                    "Fyke net","Gill net","Seine")
# write.csv(all_methods_presence, "/Users/kbja10/Github/Oneida_metabarcoding/datasets/sp_lists_comparison.csv", row.names=FALSE)

###############################################################################################
############# Part 1: Lake-wide species presence-absence ##########################
###############################################################################################

# Venn diagram of species presence in all 6 datasets
library(RColorBrewer)
library(stringr)
library(UpSetR)

cols <- c(brewer.pal(8, "Dark2"),"#386CB0","black")
my.cols <- cols[c(10,4,9,1,6,3)]
# pdf("/Users/kbja10/Documents/Cornell/Research/Oneida/Figures/sp_lists_overlap.pdf", width=10, height=7) 
sp_lists_overlap <- upset(all_methods_presence, nsets = 6, sets.bar.color = my.cols, order.by = "freq",
      mainbar.y.label = "Number of overlapping species", text.scale=c(1.7, 1.5, 1.2,1.5,1.5,1.5),
      sets.x.label = "Number of species per dataset")
# dev.off()

###############################################################################################
############# Part 2: Species accumulation curves ##########################
###############################################################################################

library(vegan)
library(scales)
sp_acc <- specaccum(eDNA_dat)
# pdf("/Users/kbja10/Documents/Cornell/Research/Oneida/Figures/spp_accum_curve_2.12.pdf", width=8, height=6)
plot(sp_acc, ci.type="poly", col=my.cols[2], lwd=3, ci.lty=0, ci.col=alpha(my.cols[2], 0.2), 
     ylab="Number of species",xlab="Number of sites", cex.axis=1.5, cex.lab=1.5, ylim = c(0,60))
sp_acc_ef <- specaccum(ef_dat[,-1])
plot(sp_acc_ef, ci.type="poly", col=my.cols[3], lwd=3, ci.lty=0, ci.col=alpha(my.cols[3], 0.3), add=T)
sp_acc_fyke <- specaccum(fyke_dat[,-1])
plot(sp_acc_fyke, ci.type="poly", col=my.cols[5], lwd=3, ci.lty=0, ci.col=alpha(my.cols[5], 0.3), add=T)
sp_acc_gillnet <- specaccum(gillnet_dat[,-1])
plot(sp_acc_gillnet, ci.type="poly", col=my.cols[6], lwd=3, ci.lty=0, ci.col=alpha(my.cols[6], 0.3), add=T)
sp_acc_seine <- specaccum(seine_dat[,-1])
plot(sp_acc_seine, ci.type="poly", col=my.cols[4], lwd=3, ci.lty=0, ci.col=alpha(my.cols[4], 0.5), add=T)
legend(0, 60, legend=c("eDNA","Electrofishing","Fyke","Gillnet","Seine"), col=my.cols[c(2,3,5,6,4)], lwd=3,bty="n")
# dev.off()

###############################################################################################
################# Part 3: Species relative proportions #################################
###############################################################################################

library(pals)
library(ggplot2)
all_methods_prop <- data.frame(species=all_methods_totals$species,
                               lapply(all_methods_totals[,3:7], function(x) x / sum(x)))

all_methods_prop <- as.data.frame(all_methods_prop %>%
  pivot_longer(!species, names_to = "gear", values_to = "prop"))
top_spp <- all_methods_prop %>% 
  arrange(desc(prop)) %>% 
  group_by(gear) %>% slice(1:8)

# Stacked + percent
species_proportions <- ggplot(top_spp, aes(fill=species, y=prop, x=gear)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=as.vector(stepped(20))) +
  xlab("Sampling gear") + ylab("Species proportion") +
  labs(fill = "Species") +
  scale_x_discrete(labels=c("eDNA","Electrofishing","Fyke net","Gill net","Seine")) +
  theme_bw()
species_proportions
# ggsave("/Users/kbja10/Documents/Cornell/Research/Oneida/Figures/all_methods_species_proportion.pdf", plot = species_proportions, dpi=600)

###############################################################################################
################# Part 4: Habitat comparison -- eDNA vs. traditional methods ########################
###############################################################################################

eDNA_metadata <- read.csv("/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/oneida_eDNA_sample_metadata.csv", header = TRUE)
nmds_presence_dat <- bind_rows(eDNA_dat, ef_dat) # combine eDNA and ef datasets
nmds_presence_dat <- nmds_presence_dat[!grepl("Y", nmds_presence_dat$Site),] # remove subsamples (labeled with "Y")
nmds_presence_dat[is.na(nmds_presence_dat)] <- 0 # indicate absence w/ 0
nmds_presence_dat[,-1][nmds_presence_dat[,-1]>0] <- 1 # indicate presence w/ 1
rownames(nmds_presence_dat) <- c(paste0("eDNA_",c(1:11,13:26)),paste0("ef_",c(1:24)))

NMDS_dat <- metaMDS(nmds_presence_dat[,-(1:2)],distance="bray",trymax = 200) # nmds w/ Sorensen index
data.scores <- as.data.frame(scores(NMDS_dat)) # turn into data frame for plotting
data.scores$Sample <- factor(rownames(data.scores)) 
data.scores$Method <-  c(rep("eDNA",25),rep("ef",24))
data.scores$Habitat <- c(eDNA_metadata$Location.type[-12], rep("Nearshore",24)) # all ef surveys are nearshore
data.scores$Habitat <- gsub("?","",data.scores$Habitat, fixed = TRUE)

nmds_plot <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size=5, stroke=2, aes(shape=Method, colour=Habitat))+ 
  theme(axis.text = element_text(colour = "black", size = 12), 
        axis.title = element_text(size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  scale_shape_manual(name = "Method", values=c(1,16)) +
  scale_color_manual(values=as.vector(stepped(14)[c(2,14,9)])) +
  labs(x = "NMDS1", colour = "Site", y = "NMDS2")
nmds_plot
# ggsave("/Users/kbja10/Documents/Cornell/Research/Oneida/Figures/nmds_presence_absence.pdf", plot = nmds_plot, dpi=600)
