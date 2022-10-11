### Oneida eDNA metabarcoding analysis
### This code imports and standardizes all fish survey datasets (eDNA, EF, seine, fyke, gillnet, and historical)
### Outputs: 1) combined gears presence-absenece, and 2) combined gears abundance
### Outputs are community datasets (vegan) for each site as well as lake-wide 
### Created by K. Andres

rm(list = ls())
historical_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/oneida_historical_species_list.csv", header = TRUE)
sp_read_count_by_site <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/sp_read_count_by_site.csv")
ef_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/Oneida_2017_comm_dataset.csv")
fyke_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/2017 oneida fyke.csv")
gillnet_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/2017 oneida gill net.csv")
seine_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/2017 oneida seine.csv")

library(dplyr)
library(tidyr)
library(stringr)

#########################################################################
###############3# Reconcile inconsistencies across datasets ############
#########################################################################

# eDNA dataset
to_replace <- c(" ","Ambloplites.rupestris","northern.hog.sucker","Percopsis.omiscomaycus","pumpkinseed.sunfish","tessallated.darter",
                "Cottidae","Cyprinidae","Cyprinus","Luxilus","Notropis","Percina","Lepisosteus")
replace_with <- c("\\.","rock.bass","northern.hogsucker","trout.perch","pumpkinseed","tessellated.darter",
                  "Cottidae.sp.","Cyprinidae.sp.","Cyprinus.sp.","Luxilus.sp.","Notropis.sp.","Percina.sp.","Lepisosteus.sp.")
names(replace_with) <- to_replace
sp_read_count_by_site$scomnames <- str_replace_all(sp_read_count_by_site$scomnames, replace_with)
eDNA_dat <- t(as.matrix(sp_read_count_by_site[,-c(1:4)])) # transpose to create community dataset
colnames(eDNA_dat) <- sp_read_count_by_site$scomnames
eDNA_dat <- data.frame(Site=rownames(eDNA_dat),eDNA_dat)
# Standardize each row (sample) to 100 individuals (account for differences in sequencing depth among samples)
eDNA_dat <- data.frame(Site=eDNA_dat$Site, eDNA_dat[,-1]/rowSums(eDNA_dat[,-1])*100)
# write.csv(eDNA_dat, "/Users/kbja10/Github/Oneida_metabarcoding/datasets/eDNA_dat.csv", row.names=FALSE)

# Electrofishing dataset
to_replace <- c("chubsucker","gizzard.shad","killifish","lepomis","pumpkinseedXgreen.sunfish","redhorse","unknown.shiner")
replace_with <- c("creek.chubsucker","American.gizzard.shad","banded.killifish","Lepomis.1","Lepomis.2","greater.redhorse","Notropis.sp.")
names(replace_with) <- to_replace
colnames(ef_dat) <- str_replace_all(colnames(ef_dat), replace_with)
ef_dat <- ef_dat %>% mutate(Lepomis.sp.=Lepomis.1+Lepomis.2, .keep = "unused")
# Keep only all-species runs (exclude predator-only runs)
ef_dat <- ef_dat[!grepl(" 2", ef_dat$Site),] # remove run 2 from each site
ef_dat$Site <- gsub(" [[:digit:]]+", "", ef_dat$Site) # remove pass number 
# write.csv(ef_dat, "/Users/kbja10/Github/Oneida_metabarcoding/datasets/ef_dat.csv", row.names=FALSE)

# Fyke dataset
fyke_dat <- fyke_dat %>% 
  mutate(largemouth.bass=largemouth.bass+largemouth.bass.YOY, .keep = "unused") %>%
  mutate(smallmouth.bass=smallmouth.bass+smallmouth.bass.YOY, .keep = "unused") %>%
  mutate(walleye=walleye+walleye.YOY, .keep = "unused") %>%
  mutate(yellow.perch=yellow.perch+yellow.perch.YOY, .keep = "unused")
to_replace <- c("american.eel","Cisco","Common.carp","darter","gizzard.shad","redhorse")
replace_with <- c("American.eel","cisco","common.carp","Etheostoma.sp.","American.gizzard.shad","greater.redhorse")
names(replace_with) <- to_replace
colnames(fyke_dat) <- str_replace_all(colnames(fyke_dat), replace_with)
fyke_dat$Other
fyke_dat <- fyke_dat %>% select(Site, American.eel:yellow.perch)
fyke_dat$satinfin.shiner <- c(rep(0,4),1,rep(0,31))
fyke_dat$bluntnose.minnow <- c(rep(0,7),1,0,0,0,1,0,0,0,0,0,2,rep(0,13),2,0,0,0,0)
fyke_dat$creek.chubsucker <- c(rep(0,18),1,0,0,3,rep(0,14))
fyke_dat$PSXGSF <- c(rep(0,26),1,0,1,rep(0,7))
fyke_dat <- fyke_dat %>% mutate(Lepomis.sp.=PSXGSF+Lepomis.YOY, .keep = "unused")
# Drop columns (species) with no occurrences
cols_to_drop = c(TRUE, colSums(fyke_dat[,-1]) > 0)
fyke_dat <- fyke_dat[,cols_to_drop]
# write.csv(fyke_dat, "/Users/kbja10/Github/Oneida_metabarcoding/datasets/fyke_dat.csv", row.names=FALSE)

# Gill net dataset
to_replace <- c("gizzard.shad","redhorse.sp.","tiger.muskellunge")
replace_with <- c("American.gizzard.shad","greater.redhorse","tiger.musky")
names(replace_with) <- to_replace
colnames(gillnet_dat) <- str_replace_all(colnames(gillnet_dat), replace_with)
gillnet_dat <- gillnet_dat %>% select(Site, walleye:round.goby)
gillnet_dat[is.na(gillnet_dat)] <- 0
# Drop columns (species) with no occurrences
cols_to_drop = c(TRUE, colSums(gillnet_dat[,-1]) > 0)
gillnet_dat <- gillnet_dat[,cols_to_drop]
# write.csv(gillnet_dat, "/Users/kbja10/Github/Oneida_metabarcoding/datasets/gillnet_dat.csv", row.names=FALSE)

# Seine dataset
to_replace <- c(" ","gizzard.shad","shiner.sp.","lepomis.sp.")
replace_with <- c("\\.","American.gizzard.shad","Notropis.sp.","Lepomis.sp.")
names(replace_with) <- to_replace
seine_dat$Species <- str_replace_all(seine_dat$Species, replace_with)
seine_dat <- as.data.frame(seine_dat %>% pivot_wider(names_from = Species, values_from = Number.Caught) %>%
                             group_by(Site) %>%
                             summarise(across(bluegill:chain.pickerel, sum, na.rm = TRUE)))
# write.csv(seine_dat, "/Users/kbja10/Github/Oneida_metabarcoding/datasets/seine_dat.csv", row.names=FALSE)


#########################################################################
######## Lake-wide datasets: get species totals across all sites ########
#########################################################################

# eDNA species read counts normalized per sample (scaled to 100 individuals per sample)
eDNA_totals <- data.frame(species=colnames(eDNA_dat[,-1]),
                          number=colSums(eDNA_dat[,-1]))
# eDNA_totals$number <- eDNA_totals$number/sum(eDNA_totals$number)
eDNA_totals <- eDNA_totals[order(eDNA_totals$species),]
eDNA_totals$gear <- rep("eDNA",nrow(eDNA_totals))

# Electrofishing species totals
ef_totals <- data.frame(species = colnames(ef_dat[,-1]),
                        number = colSums(ef_dat[,-1]))
ef_totals <- ef_totals[ef_totals$number > 0, ]
ef_totals <- ef_totals[order(ef_totals$species),]
ef_totals$gear <- rep("Electrofishing",nrow(ef_totals))

# Fyke species totals
fyke_totals <- data.frame(species = colnames(fyke_dat[,-1]),
                          number = colSums(fyke_dat[,-1]))
fyke_totals <- fyke_totals[fyke_totals$number > 0, ]
fyke_totals <- fyke_totals[order(fyke_totals$species),]
fyke_totals$gear <- rep("Fyke netting",nrow(fyke_totals))

# Gillnet species totals
gillnet_totals <- data.frame(species = colnames(gillnet_dat[,-1]),
                             number = colSums(gillnet_dat[,-1]))
gillnet_totals <- gillnet_totals[gillnet_totals$number > 0, ]
gillnet_totals <- gillnet_totals[order(gillnet_totals$species),]
gillnet_totals$gear <- rep("Gill netting",nrow(gillnet_totals))

# Seine species totals
seine_totals <- data.frame(species = colnames(seine_dat[,-1]),
                           number = colSums(seine_dat[,-1]))
seine_totals <- seine_totals[seine_totals$number > 0, ]
seine_totals <- seine_totals[order(seine_totals$species),]
seine_totals$gear <- rep("Seining", nrow(seine_totals))

# All methods totals
historical_dat$scomnames <- gsub(" ", ".", historical_dat$scomnames)
historical_totals <- data.frame(species=historical_dat$scomnames,number=historical_dat$since_1990,
                                gear=rep("historical",nrow(historical_dat)))
all_methods_abundance <- rbind(historical_totals, eDNA_totals,ef_totals,fyke_totals,
                            gillnet_totals,seine_totals)
all_methods_abundance <- as.data.frame(all_methods_abundance %>% 
                                      pivot_wider(names_from = gear, values_from = number))
all_methods_abundance[is.na(all_methods_abundance)] <- 0
all_methods_abundance <- all_methods_abundance[!grepl(".sp", all_methods_abundance$species),] # remove genus-level IDs
# write.csv(all_methods_abundance, "/Users/kbja10/Github/Oneida_metabarcoding/datasets/all_methods_abundance.csv", row.names=FALSE)

# All methods presence-absence
all_methods_presence <- all_methods_abundance[!grepl(".sp", all_methods_abundance$species),] # remove genus-level IDs
all_methods_presence[,2:7][all_methods_presence[,2:7]>0] <- 1 # indicate presence w/ 1
all_methods_presence <- left_join(all_methods_presence, historical_dat[,1:2], by=c("species"="scomnames")) %>%
  select(species.y, everything()) %>%  
  arrange(desc(all_methods_presence$historical), desc(all_methods_presence$eDNA), desc(all_methods_presence$Electrofishing), desc(all_methods_presence$`Fyke netting`), desc(all_methods_presence$`Gill netting`), desc(all_methods_presence$Seining))
all_methods_presence$species <- gsub("\\.", " ", all_methods_presence$species)
colnames(all_methods_presence) <- c("Scientific name","Common name","Historical","eDNA","Electrofishing",
                                    "Fyke netting","Gill netting","Seining")
# write.csv(all_methods_presence, "/Users/kbja10/Github/Oneida_metabarcoding/datasets/all_methods_presence.csv", row.names=FALSE)
