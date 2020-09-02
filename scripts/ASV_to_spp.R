### Oneida eDNA analysis
### Assigning species to each ASV in the ASV by site matrix (dada2 output)
### Created 4.1.2020 by Tim Lambert, updated 9.1.2020 by Kara Andres

# Read in BLAST output (in .csv format)
blast_dat <- read.csv(file = "/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/oneida_output_8.25.2020/blstn_nonchim_fmt10_scinames_Oneida_withseq_taxonimic_matched.csv", header = TRUE)
head(blast_dat)

# Read in ASV counts by sites: columns are samples, rows are ASVs
ASV_count_by_sites <- read.csv(file = "/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/oneida_output_8.25.2020/seqtab.nochim.csv", header = TRUE)
head(ASV_count_by_sites)
colnames(ASV_count_by_sites)[1] <- "ASV" # rename first column to be ASV
colnames(ASV_count_by_sites)[-1] <- sub("_", "", substring(colnames(ASV_count_by_sites)[-1], 33, 35)) # remove extra characters in sample names

# PLOT: # reads vs. # samples each ASV occurs in
ASV <- data.frame(ASV_count_by_sites, n_sites = NA, n_reads = NA)
ASV$n_sites <- sapply(1:nrow(ASV_count_by_sites), function(x) sum(ASV_count_by_sites[x,10:36]>0)) # number of G samples ASV is present
ASV$n_reads <- sapply(1:nrow(ASV_count_by_sites), function(x) sum(ASV_count_by_sites[x,10:36])) # total number of reads across all G samples

plot(log(ASV$n_reads, base = 10) ~ jitter(ASV$n_sites, 2), pch = 16, cex = 0.15,
     xlab = "Number of sites", ylab = "log(ASV read count)", main = "ASV read count by number of sites")

# PLOT: Histogram of # reads per ASV (1 site vs. many sites)
hist(log(ASV$n_reads[ASV$n_sites == 1]), breaks = 0:20, col=rgb(1,0,0,0.5),
     xlab = "log(total # of reads)", main = "Histogram of ASV read counts", xlim = c(0,15))
hist(log(ASV$n_reads[ASV$n_sites > 1]), col=rgb(0,0,1,0.5), add=T)
legend("topright", legend = c("1 site", "2+ sites"),
       fill = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))

# PLOT: Histogram of ASVs by number of sites
hist(ASV$n_sites, breaks = 0:27, main = "Histogram of ASVs by number of sites", xlab = "Number of sites", ylab = "Number of ASVs")


#### Convert ASVs to species ####
# Select species: (1) top hit, provided ≥98%. (2) second hit, if == top hit, (3) third hit, if == second hit, etc.
# If two have equal pident values, include them both and assign to highest taxonomic rank
# Add species matches (common and scientific names) as final two columns on the ASV_count_by_sites matrix
pident_threshold <- 98 # threshold of % sequence identical above which a match is considered valid
ASV_count_by_sites$seqid <- unique(blast_dat$seqid)
ASV_count_by_sites$family <- NA
ASV_count_by_sites$genus <- NA
ASV_count_by_sites$sscinames <- NA
ASV_count_by_sites$scomnames <- NA

# for each ASV, subset the associated blast hits
for(seqid in ASV_count_by_sites$seqid) {
    subdat <- blast_dat[blast_dat$seqid == seqid & blast_dat$pident > pident_threshold, ] # subset to hits > pident_threshold
    if(nrow(subdat)>0) {subdat <- subdat[subdat$pident==max(subdat$pident),]} # subset to top hit(s)
    subdat2 <- blast_dat[blast_dat$seqid == seqid, ] # subset to hits > pident_threshold
    subdat2 <- subdat2[subdat2$pident==max(subdat2$pident),]
    if (nrow(subdat)==0) { # if no hit exceeds pident_threshold
      ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, c("sscinames", "scomnames")] <- c("NO MATCH", "NO MATCH")
      if (length(unique(subdat2$genus))==1){ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "genus"] <- subdat2[1, "genus"]
      } else {ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "genus"] <- "NO MATCH"}
      if (length(unique(subdat2$family))==1){ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "family"] <- subdat2[1, "family"]
      } else {ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "family"] <- "NO MATCH"}
    } else if (length(unique(subdat$species))==1) { # if single hit or multiple hits > pident_threshold but are the same species
      ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, c("sscinames", "scomnames","genus","family")] <- subdat[1, c("sscinames", "scomnames","genus","family")]
    } else {
      ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, c("sscinames", "scomnames")] <- c("AMBIGUOUS","AMBIGUOUS")
      if (length(unique(subdat2$genus))==1){ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "genus"] <- subdat2[1, "genus"]
      } else {ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "genus"] <- "AMBIGUOUS"}
      if (length(unique(subdat2$family))==1){ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "family"] <- subdat2[1, "family"]
      } else {ASV_count_by_sites[ASV_count_by_sites$seqid == seqid, "family"] <- "AMBIGUOUS"}
    }
  }
  
# Restrict output ASVs to those with lengths in ASV_length_range
lengths <- sapply(1:nrow(ASV_count_by_sites), FUN = function(x) nchar(as.character(ASV_count_by_sites[x,"ASV"])))
ASV_length_range <- c(160, 192) # lengths of sequences to look at (to eliminate bacteria, etc.)
ASV_count_by_sites <- ASV_count_by_sites[lengths > ASV_length_range[1]-1 & lengths < ASV_length_range[2]+1, ]

nrow(ASV_count_by_sites[ASV_count_by_sites$scomnames=="AMBIGUOUS",]) # 73
nrow(ASV_count_by_sites[ASV_count_by_sites$scomnames=="NO MATCH",]) # 725

# Sum counts for each species (combine multiple hits)
sp_count_by_site <- rowsum(ASV_count_by_sites[,2:36], group = ASV_count_by_sites$scomnames, na.rm = TRUE)

# Write data files
# write.csv(ASV_count_by_sites, "/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/ASV_count_by_sites.csv")
# write.csv(sp_count_by_site, "/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/eDNA_counts_pident_98_inc_blanks.csv") # read counts by species including the blanks
