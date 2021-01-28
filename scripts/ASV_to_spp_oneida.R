### Oneida eDNA analysis
### This code uses BLAST and DADA2 outputs to assign species to each ASV 
### The output is an ASV-read-count-by-site matrix with taxonomic information for each ASV
### Created 4.1.2020 by Tim Lambert, updated 9.1.2020 by Kara Andres

# Read in BLAST output
blast_dat <- read.csv(file = "/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/oneida_output_8.25.2020/blstn_nonchim_fmt10_scinames_Oneida_withseq_taxonimic_matched.csv", header = TRUE)
head(blast_dat)

# Read in ASV counts by sites: columns are samples, rows are ASVs
ASV_count_by_sites <- read.csv(file = "/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/oneida_output_8.25.2020/seqtab.nochim.csv", header = TRUE)
colnames(ASV_count_by_sites)[1] <- "ASV" # rename first column to be ASV
colnames(ASV_count_by_sites)[-1] <- sub("_", "", substring(colnames(ASV_count_by_sites)[-1], 33, 35)) # remove extra characters in sample names
head(ASV_count_by_sites)

### Convert ASVs to species ####
# Select species: (1) subset to the top blast hit(s), provided ≥98%. 
#                 (2) if only one top hit or multiple top hits are the same species, assign species match (common and scientific names)
#                 (3) if multiple top hits are different species, label species "AMBIGUOUS"
#                 (4) use same rules to assign genus and family ranks, using subset not restricted to ≥98%

pident_threshold <- 98 # threshold of % sequence identical above which a match is considered valid
ASV_count_by_sites <- data.frame(ASV_count_by_sites, seqid=unique(blast_dat$seqid),
                                 family=NA, genus=NA, species=NA, scomnames=NA)
for(seqid in ASV_count_by_sites$seqid) { # for each ASV, subset to the associated blast hits
  subdat <- blast_dat[blast_dat$seqid == seqid & blast_dat$pident >= pident_threshold, ] # subset to hits > pident_threshold
  if(nrow(subdat)>0) {subdat <- subdat[subdat$pident==max(subdat$pident),]} # subset to top hit(s)
  subdat2 <- blast_dat[blast_dat$seqid == seqid, ] # create second subset for higher taxonomic matches
  subdat2 <- subdat2[subdat2$pident==max(subdat2$pident),] # subset to top hit(s)
  if (nrow(subdat)==0) { # if no hit exceeds pident_threshold
    ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, c("species", "scomnames")] <- c("NO MATCH", "NO MATCH")
    if (length(unique(subdat2$genus))==1){ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "genus"] <- subdat2[1, "genus"]
    } else {ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "genus"] <- "NO MATCH"}
    if (length(unique(subdat2$family))==1){ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "family"] <- subdat2[1, "family"]
    } else {ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "family"] <- "NO MATCH"}
  } else if (length(unique(subdat$species))==1) { # if single hit or multiple hits > pident_threshold but are the same species
    ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, c("species", "scomnames","genus","family")] <- subdat[1, c("species", "scomnames","genus","family")]
  } else {
    ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, c("species", "scomnames")] <- c("AMBIGUOUS","AMBIGUOUS")
    if (length(unique(subdat$genus))==1){ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "genus"] <- subdat[1, "genus"]
    } else {ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "genus"] <- "AMBIGUOUS"}
    if (length(unique(subdat$family))==1){ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "family"] <- subdat[1, "family"]
    } else {ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "family"] <- "AMBIGUOUS"}
  }
}

# Restrict output ASVs to those with lengths in ASV_length_range
lengths <- sapply(1:nrow(ASV_count_by_sites), FUN = function(x) nchar(as.character(ASV_count_by_sites[x,"ASV"])))
ASV_length_range <- c(160, 192) # lengths of sequences to look at (to eliminate bacteria, etc.)
ASV_count_by_sites <- ASV_count_by_sites[lengths > ASV_length_range[1]-1 & lengths < ASV_length_range[2]+1, ]

nrow(ASV_count_by_sites) # 2982 ASVs
nrow(ASV_count_by_sites[ASV_count_by_sites$scomnames=="AMBIGUOUS",]) # 73 abiguous matches
nrow(ASV_count_by_sites[ASV_count_by_sites$scomnames=="NO MATCH",]) # 725 no match exceeding 98% identity

# Write data files
# write.csv(ASV_count_by_sites, "/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/ASV_count_by_sites.csv", row.names=FALSE)
