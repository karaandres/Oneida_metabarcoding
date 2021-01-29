### Oneida eDNA analysis
### This code filters the eDNA dataset by (1) removing ASVs with low read counts
###                                       (2) removing non-target species, and 
###                                       (3) adjusting dataset for any false positives in the controls
### The output is a species-read-count-by-site matrix with controls excluded 
### Created 4.1.2020 by Kara Andres, updated 9.1.2020

# Read in ASV counts by sites: columns are samples, rows are ASVs, taxonomic information in last columns
ASV_count_by_sites <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/ASV_count_by_sites.csv")
ASV_count_by_sites <- subset(ASV_count_by_sites, select=-ASV) # remove unneeded columns

### Remove ASVs with reads below threshold of 0.1% of all reads in a sample ###
ASV_count_by_sites_0.001 <- matrix(NA, nrow = nrow(ASV_count_by_sites), ncol = ncol(ASV_count_by_sites)-5) # create empty matrix to fill
for (i in 1:(ncol(ASV_count_by_sites)-5)){ # for each column (sample)
  for (j in 1:nrow(ASV_count_by_sites)){ # for each row (ASV)
    if (sum(ASV_count_by_sites[1:nrow(ASV_count_by_sites),i])>0) { # if number of reads in sample is >0
      if (ASV_count_by_sites[j,i]/sum(ASV_count_by_sites[1:nrow(ASV_count_by_sites),i])<0.001){
        ASV_count_by_sites_0.001[j,i] <- 0 # replace reads below 0.1% of sample reads with a 0
      } else ASV_count_by_sites_0.001[j,i] <- ASV_count_by_sites[j,i] # otherwise they stay the same
    } else ASV_count_by_sites_0.001[j,i] <- ASV_count_by_sites[j,i]
  }
}

ASV_count_by_sites_0.001 <- cbind(ASV_count_by_sites_0.001, ASV_count_by_sites[,36:40]) # add taxonomic info back in
colnames(ASV_count_by_sites_0.001) <- colnames(ASV_count_by_sites)
ASV_count_by_sites_0.001 <- ASV_count_by_sites_0.001[rowSums(ASV_count_by_sites_0.001[,1:35])>0,] # remove any ASVs with no reads across all samples

### Remove non-target species
unique(ASV_count_by_sites_0.001$scomnames) # remove human, cattle, saithe, red-bellied piranha
ASV_count_by_sites_0.001 <- ASV_count_by_sites_0.001[!(ASV_count_by_sites_0.001$scomnames=="cattle"|ASV_count_by_sites_0.001$scomnames=="human"|ASV_count_by_sites_0.001$scomnames=="saithe"|ASV_count_by_sites_0.001$scomnames=="red-bellied piranha"),]

### Adjust reads for species present in blanks
cumul_control <- data.frame(reads=ASV_count_by_sites_0.001$BL1+ASV_count_by_sites_0.001$BL2+
                              ASV_count_by_sites_0.001$BL3+ASV_count_by_sites_0.001$X1BL+
                              ASV_count_by_sites_0.001$X27G+ASV_count_by_sites_0.001$X12G,
                            family=ASV_count_by_sites_0.001$family,
                            genus=ASV_count_by_sites_0.001$genus,
                            species=ASV_count_by_sites_0.001$scomnames)

# the majority of species in the dataset do not show up in any of the controls
hist(cumul_control$reads, breaks = 75, xlab = "Number of reads", main = "")

# calculating the mean number of non-zero reads in negative controls
nonzero_cumul_control <- cumul_control[cumul_control$reads>0,]
mean(nonzero_cumul_control$reads[nonzero_cumul_control$reads>0], na.rm=T) # mean = 23 reads

### Remove ASVs with reads below mean number of non-zero reads summed across controls ###
ASV_count_by_sites_0.001_adjusted <- matrix(NA, nrow = nrow(ASV_count_by_sites_0.001), ncol = ncol(ASV_count_by_sites_0.001)-5) # create empty matrix to fill
for (i in 1:(ncol(ASV_count_by_sites_0.001)-5)){ # for each column (sample)
  for (j in 1:nrow(ASV_count_by_sites_0.001)){ # for each row (ASV)
    if (ASV_count_by_sites_0.001[j,i]<23) { # if number of reads in sample is <23
      ASV_count_by_sites_0.001_adjusted[j,i] <- 0 # replace reads <23 with a 0
    } else ASV_count_by_sites_0.001_adjusted[j,i] <- ASV_count_by_sites_0.001[j,i] # otherwise they stay the same
  }
}
ASV_count_by_sites_0.001_adjusted <- cbind(ASV_count_by_sites_0.001_adjusted, ASV_count_by_sites_0.001[,36:40]) # add taxonomic info back in
colnames(ASV_count_by_sites_0.001_adjusted) <- colnames(ASV_count_by_sites_0.001)

# Remove negative controls
ASV_count_by_sites_cleaned <- ASV_count_by_sites_0.001_adjusted[ , -which(names(ASV_count_by_sites_0.001_adjusted) %in% 
                                                                            c("BL1","BL2","BL3","X1BL","X12G","X27G"))] # remove negative controls
ASV_count_by_sites_cleaned <- ASV_count_by_sites_cleaned[rowSums(ASV_count_by_sites_cleaned[,1:29])>0,] # remove any ASVs with no reads across all samples

# Assign lowest taxonomic level available to each ASV as the "species"
ASV_count_by_sites_lowest_tax <- ASV_count_by_sites_cleaned
for (i in 1:nrow(ASV_count_by_sites_cleaned)){ # for each ASV
  # if no species match or ambiguous assignment, assign genus level taxonomic info
  if (ASV_count_by_sites_cleaned$species[i]=="NO MATCH"|ASV_count_by_sites_cleaned$species[i]=="AMBIGUOUS"){
    ASV_count_by_sites_lowest_tax$species[i] <- ASV_count_by_sites_cleaned$genus[i]
    ASV_count_by_sites_lowest_tax$scomnames[i] <- ASV_count_by_sites_cleaned$genus[i]
  }
  # if no genus level taxonomic assignment, assign family level taxonomic info
  if (ASV_count_by_sites_cleaned$genus[i]=="AMBIGUOUS"){
    ASV_count_by_sites_lowest_tax$species[i] <- ASV_count_by_sites_cleaned$family[i]
    ASV_count_by_sites_lowest_tax$scomnames[i] <- ASV_count_by_sites_cleaned$family[i]
  }
}

# Sum counts for each species (combine ASVs wiht same taxonomic assigment)
# write.csv(sp_count_by_site, "/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/sp_read_count_by_site.csv") # read counts per site by species excluding blanks
library(dplyr)
sp_count_by_site <- ASV_count_by_sites_lowest_tax %>% 
  group_by(species) %>% 
  summarise(across(where(is.numeric), sum)) %>%
  left_join(distinct(ASV_count_by_sites_lowest_tax[, c("family", "genus", "species", "scomnames")], species,.keep_all = TRUE)) %>%
  select(c(family, genus, scomnames), everything()) %>%
  arrange(scomnames)
# write.csv(sp_count_by_site, "/Users/kbja10/Github/Oneida_metabarcoding/datasets/sp_read_count_by_site.csv", row.names=FALSE) # read counts per site by species excluding blanks
