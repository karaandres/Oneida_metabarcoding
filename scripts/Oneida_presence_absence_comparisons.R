### Oneida eDNA metabarcoding analysis
### This code compares electrofishing data, historical records, and eDNA metabarcoding datasets 
### Created 1.28.2020 by Kara Andres

sp_read_count_by_site <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/sp_read_count_by_site.csv")
oneida_ef <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/Oneida_2017_comm_dataset.csv")
historical_species_dataset <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/oneida_historical_species_list.csv", header = TRUE)

# Reconcile inconsistencies in species names across datasets 
library(dplyr)
# electrofishing dataset
to_replace <- c("\\.","chubsucker","gizzard shad","killifish","lepomis","pumpkinseedXgreen.sunfish","redhorse","unknown shiner")
replace_with <- c(" ","creek chubsucker","American gizzard shad","banded killifish","Lepomis.1","Lepomis.2","Moxostoma","Notropis")
names(replace_with) <- to_replace
colnames(oneida_ef) <- str_replace_all(colnames(oneida_ef), replace_with)
oneida_ef <- oneida_ef %>% mutate(Lepomis=Lepomis.1+Lepomis.2, .keep = "unused")
# eDNA dataset
to_replace <- c("Ambloplites rupestris","Gobius","Labidesthes","northern hog sucker","Percopsis omiscomaycus","pumpkinseed sunfish","tessallated darter")
replace_with <- c("rock bass","round goby","brook silverside","northern hogsucker","trout perch","pumpkinseed","tessellated darter")
names(replace_with) <- to_replace
sp_read_count_by_site$scomnames <- str_replace_all(sp_read_count_by_site$scomnames, replace_with)

############# Part 1: Lake-wide species presence-absence ##########################
# Venn diagram of species presence in all 3 datasets
library(VennDiagram)
library(RColorBrewer)
library(stringr)

# species lists from 3 datasets
hist_spp <- as.character(historical_species_dataset$historical)
ef_spp <- as.character(colnames(oneida_ef[-1])) 
eDNA_spp <- as.character(sp_read_count_by_site$scomnames)


# Make Venn diagram
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(x = list(ef_spp, hist_spp, eDNA_spp),
             category.names = c("Electrofishing" , "Historical" , "eDNA"),
             filename = '/Users/kbja10/Documents/Cornell/Research/Oneida/Figures/venn_diagramm_1.28.png', output=TRUE,
             imagetype="png" , height = 3000, width = 3000, resolution = 300, lwd = 2, lty = 'blank', 
             fill = myCol, cex = 4, cat.cex = 2, fontface = "bold", fontfamily = "sans")

############# Part 2: Species accumulation curves ##########################
library(vegan)
library(scales)
sp_read_count_by_site_t <- t(as.matrix(sp_read_count_by_site[,-c(1:4)])) # transpose to create community dataset
colnames(sp_read_count_by_site_t) <- sp_read_count_by_site[,4]
sp_acc <- specaccum(sp_read_count_by_site_t)
# png("/Users/kbja10/Documents/Cornell/Research/Oneida/Figures/spp_accum_curve_1.28.png", width=8, height=6, units="in", res=600)
plot(sp_acc, ci.type="poly", col="blue", lty=3, lwd=3, ci.lty=0, ci.col=alpha("gray", 0.2), 
     ylab="Number of species",xlab="Number of sites", cex.axis=1.5, cex.lab=1.5, ylim = c(0,60))
plot(sp_acc, ci.type="poly", col="#56B4E9", lty=3, lwd=3, ci.lty=0, ci.col=alpha("gray", 0.2), 
     ylab="Number of species",xlab="Number of sites", cex.axis=1.5, cex.lab=1.5, ylim = c(0,60))
sp_acc_ef <- specaccum(oneida_ef[,-1])
plot(sp_acc_ef, ci.type="poly", col="#009E73", lwd=3, ci.lty=0, ci.col=alpha("gray", 0.5), add=T)
legend(1, 60, legend=c("eDNA samples", "Electrofishing survey"), col=c("#56B4E9","#009E73"), lty=c(3,1), lwd=3, cex=1.2)
# dev.off()

################# Part 3: Site occupancy -- eDNA vs. EF comparison ################################
# Site occupancy (proportion of sites species is detected)
# Collapse reads/species counts by site (triplicate eDNA samples/multiple EF passes)
library(dplyr)
oneida_ef$Site <- gsub(" [[:digit:]]+", "", oneida_ef$Site) # remove pass number 
oneida_ef_sites <- oneida_ef %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), sum))

sp_read_count_sites <- data.frame(sp_read_count_by_site[,1:4], midlake_Buoy113=sp_read_count_by_site$X1G, nearshore_WilsonPoint=sp_read_count_by_site$X2G,
                                  nearshore_SouthBay_EF=sp_read_count_by_site$X3G+sp_read_count_by_site$X5G,
                                  inlet_OneidaCreek=sp_read_count_by_site$X4G, inlet_FishCreek=sp_read_count_by_site$X6G, nearshore_NorthBay=sp_read_count_by_site$X7G,
                                  nearshore_Jewel_EF=sp_read_count_by_site$X8G+sp_read_count_by_site$X9G,
                                  nearshore_Cleveland_EF=sp_read_count_by_site$X10G+sp_read_count_by_site$X11G, midlake_Buoy123=sp_read_count_by_site$X13G,
                                  nearshore_Constantia_EF=sp_read_count_by_site$X14G+sp_read_count_by_site$X15G, midlake_Buoy130=sp_read_count_by_site$X16G,
                                  nearshore_PoddygutBay_EF=sp_read_count_by_site$X17G+sp_read_count_by_site$X18G, inlet_BigBayCreek=sp_read_count_by_site$X19G,
                                  nearshore_LowerSouthBay_EF=sp_read_count_by_site$X20G+sp_read_count_by_site$X21G, inlet_ChittenangoCreek=sp_read_count_by_site$X22G,
                                  nearshore_FisherBay=sp_read_count_by_site$X23G+sp_read_count_by_site$X24G,nearshore_ShackeltonPoint_EF=sp_read_count_by_site$X25G+sp_read_count_by_site$X26G)