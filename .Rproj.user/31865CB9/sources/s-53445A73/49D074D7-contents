### Optimal gear combinations for surveying Oneida Lake fishes
### Created by T. Lambert
### This code identifies the Pareto-optimal combination of methods (gears) to survey fish communities. Specifically, it resamples 2017 eDNA, electrofishing, fyke net, gill net, and seine net data from Oneida Lake to estimate how to allocate effort among gears in a way that maximizes the mean number of species detected. The Pareto frontiers are calculated for combinations of all methods as well as for combinations of only traditional (non-eDNA) methods, and both are compared to species accumulation curves for each gear alone.
### Inputs:
###   Standardized 2017 gear files (eDNA_dat.csv, ef_dat.csv, fyke_dat.csv, gillnet_dat.csv, and seine_dat.csv).
### Outputs:
###   (1) The Pareto frontiers for all gears and for all traditional gears, potted alongside species accumulation curves for each gear independently.
###   (2) Stacked area plots indicating the optimal allocation to each survey type, using two alternative approaches: (i) taking discrete samples via resampling without replacement, and (ii) [not yet implemented] an extension to continuous effort by interpolation of resampling with replacement.
###   (3) Heatmap of gear biases, in terms of detection fraction for each species.



#### Clear the working environment and load required packages ####

rm(list = ls())
library(here) # for setting paths
library(dplyr) # for use of bind_rows()
library(arrangements) # for use of combinations()
library(reshape) # for use of melt()
library(secr) # for use of strip.legend()
library(RColorBrewer) # for color palettes



#### Read in and format data ####

# Read in the standardized separate gear files
eDNA_dat <- read.csv(here("datasets", "eDNA_dat.csv"))
ef_dat <- read.csv(here("datasets", "ef_dat.csv"))                  
fyke_dat <- read.csv(here("datasets", "fyke_dat.csv"))                  
gillnet_dat <- read.csv(here("datasets", "gillnet_dat.csv"))                  
seine_dat <- read.csv(here("datasets", "seine_dat.csv"))

# Restrict eDNA samples to the "green" (large-volume) samples
eDNA_dat <- subset(eDNA_dat, grepl("G", eDNA_dat$Site))

# Restrict electrofishing samples to 15-minute all-fish runs (i.e., exclude hour-long predator runs)
ef_dat <- ef_dat[((1:nrow(ef_dat)) %% 3) %in% c(1,0), ] # NOTE: This requires that electrofishing samples be grouped by site, each with 3 samples in the following order: (1) all-fish run, (2) predator run, (3) all-fish run. A more robust approach would require upstream changes in the data preparation to retain information about the sample type.

## No longer used!
# # Aggregate electrofishing samples by site, so each "sample" contains all 3 electrofishing runs per site (2 all-fish samples & 1 hour-long predator sample)
# ef_dat <- aggregate(ef_dat[,-1], by = list(ef_dat$Site), FUN = "sum"); names(ef_dat)[1] <- "Site"

# Add a "Method" column
eDNA_dat$Method <- "eDNA"
ef_dat$Method <- "ef"
fyke_dat$Method <- "fyke"
gillnet_dat$Method <- "gillnet"
seine_dat$Method <- "seine"

# Read in species names (used as a common-to-scientific name look-up table)
all_methods_presence <- read.csv(here("datasets", "all_methods_presence.csv"))
spp_name_table <- all_methods_presence[,c("Scientific.name","Common.name")]
rm(all_methods_presence)

# Combine all data into a single data frame
dat <- bind_rows(eDNA_dat, ef_dat, fyke_dat, gillnet_dat, seine_dat)
dat[is.na(dat)] <- 0 # count for species not detected by a method is 0

# USER CHOICE: Restrict to species-level assignments
# Comment out this line if genus-level assignments are to be included
dat <- dat[,colnames(dat)[!grepl(pattern = "[.]sp[.]", colnames(dat))]]

dat_by_method <- split(dat, dat$Method) # split dat data frame into a list by gear
rm(eDNA_dat, ef_dat, fyke_dat, gillnet_dat, seine_dat) # clean up the workspace



#### Define parameters ####

# effort_per_sample <- c(eDNA = 5, ef = 5, fyke = 1.5, gillnet = 3, seine = 5) # a test scenario that makes things a bit more interesting, i.e., multiple gears contribute to the optimal allocation of effort.
effort_per_sample <- c(eDNA = 3.8, ef = 2.5, fyke = 10.5, gillnet = 12.5, seine = 10.5) # for each gear, the approximate effort (in person-hours of work) required to obtain one sample
## Current estimates:
## eDNA: { (2 people x 12 hours fieldwork) + (12 hours filtration [actually 30 hrs but much less with streamlined protocols???]) + 40 hours lab work } / 20 samples = 3.8 work hours
## Electrofishing: 2-man crew x (15 min. electrofishing + 15 min. travel + 15 min. work up game fish in boat + 30 min. work up all fish on shore) = 2.5 work hours
## Fyke net: 2-man crew x (1 hr. deploy + 1 hr. retrieve + 15 min. travel) + 6 man-hours working up fish on shore [???] = 10.5 work hours
## Gill net: 2-man crew x (1 hr. deploy + 1 hr. retrieve + 15 min. travel) + 4-man crew x 2 hrs. work up fish on shore = 12.5 work hours
## Seine net: 2-man crew x (~2 hours [???] + 15 minutes travel) + 3-man crew x 2 hrs work up fish on shore = 10.5 work hours

# TO DO -- DECISIONS: Should we combine the two 15-minute all-fish electrofishing surveys into a single sample? Note that smaller-effort samples are favored by this algorithm since it's assumed all samples are taken from random locations, so 2 15-minute electrofishing samples may cover much more territory than a single combined 30-minute effort. Questions: Are travel costs--especially between small samples--accurately represented? In real life, are there ways to adjust the "size"/"effort" of a single sample for some/all gear types? If so, this is a separate decision variable that could be optimized.

effort_max <- 100 # maximum effort


## Warnings:
# 1. Print a warning if the names (& order) of effort_per_sample and dat_by_method don't match.
warnings <- character()
if(sum(!names(dat_by_method) == names(effort_per_sample)) > 0) {
  warnings <- c(warnings, "Warning: gear names of effort_per_sample and dat_by_method do not match.")
}
# 2. Print a warning if the maximum effort that can be obtained by resampling the data without replacement falls short of effort_max for one or more gears. 
n_spls <- as.vector(table(dat$Method)) # number of samples per method
indicator <- (n_spls*effort_per_sample) < effort_max
if(sum(indicator) > 0) {
  warnings <- c(warnings, paste("Warning: effort_max exceeds available data from which to sample for the following gear(s): ", paste(names(effort_per_sample)[indicator], collapse = ", "), ".", sep = ""))
}
rm(indicator)
print(warnings)


#### Create list of all possible sample combinations ####
## (with the condition that total effort must be less than effort_max) ##
k_max <- max(floor(effort_max / effort_per_sample)) # maximum number of samples to select per survey/combination (auto-generated so that resampling achieves effort_max for every method and all their possible combinations)
comb_list <- list(as.integer(NULL)) # initialize cumulative list of all sample combinations; first entry is the combination with no samples of any sort
for(k in 1:k_max) { # loop through different values of k (# samples selected per combination)
  combs <- combinations(freq = floor(effort_max/effort_per_sample), k = k, x = 1:5)
  effort_spls <- sapply(1:nrow(combs), function(x) sum(effort_per_sample[combs[x,]])) # identify total effort per sample
  combs <- subset(combs, effort_spls <= effort_max) # cap total effort per sample
  # TO DO: ALSO CAP AT TOTAL NUMBER OF SAMPLES PER METHOD (e.g., >9 seine samples in a combination is not possible given sampling without replacement)
  comb_list <- c(comb_list, as.list(data.frame(t(combs)))) # append new combs to cumulative sample list
}
rm(combs, effort_spls, k_max)

effort_combs <- sapply(1:length(comb_list), function(x) sum(effort_per_sample[comb_list[[x]]])) # recalculate total effort per sample




#### Function definitions ####

### Function definition: pareto() 
## Finds the Pareto front (set of non-dominated solutions).
## Output: A list with two components:
##   front: a data frame of x- and y-values for the Pareto front, with effort ordered from low to high.
##   indices: corresponding indices of the Pareto-optimal solutions (positions in the original x and y vector arguments).
## Inputs:
##   x: x-values (effort)
##   y: y-values (# spp. detected, or more generally an objective to be maximized)
##   condition: integer or logical vector indicating the indices of observations to be kept
pareto = function(x, y, condition = NA) {
  d <- data.frame(x,y)
  rownames(d) <- 1:length(x)
  if(!is.na(condition)[1]) {
    d <- d[condition, ]
  }
  D <- d[order(d$x,d$y,decreasing=c(FALSE,TRUE)),]
  front <- D[which(!duplicated(cummax(D$y))),]
  return(list(front = front, indices = as.integer(rownames(front))))
}




## Function definition: idCombs()
## Identifies the sampling combinations that contain only the specified methods.
## Output:
##   a logical vector indicating which sampling combinations contain only the specified methods (i.e., excluding those that use other methods).
## Inputs:
##   combs: a list of possible sampling combinations, each a vector with length equal to the number of samples, and each element representing the method used for that sample (an integer 1, 2, 3, ...)
##   keep: vector of methods (integers 1, 2, 3, ...) to include; i.e., non-listed methods are excluded
idCombs = function(combs = comb_list, methods = 1:5, keep = 1:5) {
  if(sum(keep %in% methods) == length(methods)) {
    return(rep(TRUE, length(combs)))
  } else if(sum(keep %in% methods) == 0) {
    print("Warning: must keep at least one method.")
    return(list())
  } else {
    exclude <- methods[!(methods %in% keep)]
    logical_vec <- sapply(1:length(combs), function(x) sum(combs[[x]] %in% exclude) == 0)
    return(logical_vec)
  }
}




#### FUNCTION DEFINITION: conv_hull() ####
## (Code adapted from https://stats.stackexchange.com/questions/11919/convex-hull-in-r.)
## Inputs:
##   front: a list returned by front, with first element a dataframe containing the x-y values of the Pareto frontier and second element indicating their indices in the original combination list.
## Outputs:
##   hull: a list analogous to front, but only retaining those points on the convex hull of the front
conv_hull <- function(front) {
  require(grDevices)
  df <- front$front # store data frame of x and y values
  con.hull.pos <- chull(df) # find positions of convex hull
  hull <- list(hull = df[con.hull.pos,], # coordinates of the convex hull
               indices = front$indices[con.hull.pos]) # indices of the convex hull (in the original combination list, not in the front)
  return(hull)
}












# #### METHOD I: simulated sample combinations ####
# # NOTE: For testing purposes only. Use of Method II is preferred. 
# # Perform r simulations per sampling combination and calculate the mean & stdev detected species richness
# r <- 10 # maximum number of replicate simulations per sample combination (FEWER ONLY IF FEWER THAN r SIMULATIONS ARE POSSIBLE; NOT YET ENFORCED -- DO THIS)
# 
# n_spp_detected <- matrix(data = NA, nrow = length(comb_list), ncol = r) # initialize matrix to contain the numbers of species detected for each sample combination (row) and replicate (column)
# 
# # USEFUL OBJECTS: effort_combs, comb_list, n_spls, dat, dat_by_method
# 
# for(i in 1:length(comb_list)) { # for each sampling combination...
#   n_to_select <- sapply(1:5, function(x) sum(comb_list[[i]]==x)) # number of samples per method/gear
#   for(j in 1:r) { # for each of r replicates...
# 
#     for(m in 1:5) { # for each method/gear
#       # Randomly select the appropriate number of samples for each gear
#       if(m==1) {
#         dat_temp <- dat_by_method[[m]][sample(n_spls[m], n_to_select[m]), ] # initialize data frame containing simulated data for the current replicate
#       } else {
#         dat_temp <- rbind(dat_temp, dat_by_method[[m]][sample(n_spls[m], n_to_select[m]), ])
#       }
#     }
#     dat_temp <- dat_temp[ ,!colnames(dat_temp) %in% c("Method","Site")]
#     n_spp_detected[i,j] <- sum(colSums(dat_temp) > 0) # calculate the number of species detected
#     
#   }
# }
# 
# 
# x <- effort_combs # effort for each sample combination
# y <- rowSums(n_spp_detected)/r # average number of species detected




#### METHOD II: Exact calculation of expected number of species detected ####
## This is analogous to the argument <method = "exact"> in specaccum() from package vegan, except it has been extended to account for combinations of multiple methods. It uses resampling without replacement.
# Inputs:
## comb_list
## dat, dat_by_method
## 
spp_names <- colnames(dat)[!colnames(dat) %in% c("Method","Site")]
gear_names <- names(effort_per_sample)
n_spp <- length(spp_names)
n_gears <- length(unique(dat$Method))
n_combs <- length(comb_list)

# A[j] is the total number of samples in the dataset for gear j
A <- n_spls; names(A) <- gear_names

# a[k,j] is the number of samples selected for gear j in combination k
a <- matrix(data = NA, nrow = n_combs, ncol = n_gears)
colnames(a) <- gear_names
for(i in 1:n_combs) {
  a[i,] <- sapply(1:n_gears, function(x) sum(comb_list[[i]] == x))
}

# n[i,j] is the number of samples of gear j that contain species i
n <- matrix(data = NA, nrow = n_spp, ncol = n_gears)
rownames(n) <- spp_names; colnames(n) <- gear_names
for(i in 1:n_spp) {
  sp <- spp_names[i]
  n[i,] <- sapply(1:n_gears, function(x) sum(dat_by_method[[x]][,sp] > 0))
}


# Calculate multiplication factor for 55 spp X 5 methods X number of combinations
# NOTE: This is the most computationally intensive step! May take a few minutes or more depending on parameters.
f <- array(data = 1, dim = c(n_spp, n_gears, n_combs))
for(k in 1:n_combs) {
  for(i in 1:n_spp) {
    for(j in 1:n_gears) {
      f[i,j,k] <- ifelse(A[j] - n[i,j] < a[k,j], # See Methods for explanation of this equation.
                         as.numeric(n[i,j] == 0),
                         exp(lchoose(A[j]-n[i,j], a[k,j]) - lchoose(A[j], a[k,j])))
    }
  }
}

expected_spp_combs <- rep(NA, n_combs) # calculate the expected number of species for each sample combination
for(k in 1:n_combs) {
  expected_spp_combs[k] <- sum(1 - f[,1,k]*f[,2,k]*f[,3,k]*f[,4,k]*f[,5,k])
}

x <- effort_combs # effort for each sample combination
y <- expected_spp_combs # expected number of species detected (exact calculation) for each sample combination






#### CALCULATE & PLOT THE PARETO FRONTIER ####

front_all <- pareto(x, y)
front_trad <- pareto(x, y, condition = idCombs(keep = 2:5))
front_eDNA <- pareto(x, y, condition = idCombs(keep = 1))
front_ef <- pareto(x, y, condition = idCombs(keep = 2))
front_fyke <- pareto(x, y, condition = idCombs(keep = 3))
front_gillnet <- pareto(x, y, condition = idCombs(keep = 4))
front_seine <- pareto(x, y, condition = idCombs(keep = 5))

# PLOT 1(A): the Pareto frontiers both for all gears and for all traditional (non-eDNA) gears are potted alongside species accumulation curves for each gear independently.
cols <- c(brewer.pal(8, "Dark2"),"#386CB0")
my.cols <- cols[c(8,4,9,1,6,3,2)]

plot(front_all$front, xlim = c(0, effort_max), ylim = c(0, max(y)),
     pch = 1, col = my.cols[1], type = "b", lwd = 1,
     xlab = "Effort (work hours)", ylab = "Expected number of species detected",
     main = "Pareto frontiers for combined-method surveys")
points(front_eDNA$front, pch = 16, col = my.cols[2], lwd = 2, type = "b")
points(front_ef$front, pch = 16, col = my.cols[3], lwd = 2, type = "b")
points(front_seine$front, pch = 16, col = my.cols[4], lwd = 2, type = "b")
points(front_fyke$front, pch = 16, col = my.cols[5], lwd = 2, type = "b")
points(front_gillnet$front, pch = 16, col = my.cols[6], lwd = 2, type = "b")
points(front_trad$front, pch = 1, col = my.cols[7], lwd = 1, type = "b")

legend("topleft", legend = c("eDNA", "Electrofishing", "Seine", "Fyke", "Gillnet", "Combined traditional", "Combined eDNA and traditional"),
       pch = c(rep(16,5),1,1),
       col = my.cols[c(2:7,1)], cex = 0.7)



#### CALCULATE AND PLOT THE CONVEX HULLS OF THE PARETO FRONTIERS ####
## Enforce concavity of the Pareto frontier by eliminating points that are dominated by linear combination of other points. The benefits of enforcing concavity of the front include: (1) eliminated points are generally less efficient (in terms of # spp. per work hour), so less likely to be chosen by managers, and (2) noise in the optimal gear allocation plots is reduced. (See continuous version for an alternative remedy to the second problem.)

hull_all <- conv_hull(front_all)
hull_trad <- conv_hull(front_trad)

# PLOT 1(B): the convex hulls of Pareto frontiers -- otherwise same as the last plot
cols <- c(brewer.pal(8, "Dark2"),"#386CB0")
my.cols <- cols[c(8,4,9,1,6,3,2)]

pdf(here("figures","convex_Pareto.pdf"), width = 8, height = 6)
plot(hull_all$hull, xlim = c(0, effort_max), ylim = c(0, max(y)),
     pch = 1, col = my.cols[1], type = "b", lwd = 1,
     xlab = "Effort (work hours)", 
     ylab = "Expected number of species detected")
points(front_eDNA$front, pch = 16, col = my.cols[2], lwd = 2, type = "b")
points(front_ef$front, pch = 16, col = my.cols[3], lwd = 2, type = "b")
points(front_seine$front, pch = 16, col = my.cols[4], lwd = 2, type = "b")
points(front_fyke$front, pch = 16, col = my.cols[5], lwd = 2, type = "b")
points(front_gillnet$front, pch = 16, col = my.cols[6], lwd = 2, type = "b")
points(hull_trad$hull, pch = 1, col = my.cols[7], lwd = 1, type = "b")

legend("topleft", legend = c("eDNA", "Electrofishing", "Seine", "Fyke", "Gillnet", "Combined traditional", "Combined eDNA and traditional"),
       pch = c(rep(16,5),1,1),
       col = my.cols[c(2:7,1)], cex = 0.7)

dev.off()




#### (2) PLOT 2: Stacked area graph -- optimal gear allocation ####
# (i) Discrete samples: for a given sampling effort, the optimal allocation to each survey type.
# https://r-graphics.org/recipe-line-graph-stacked-area
library(ggplot2)

## Function definition: stacked_area_plot() ####
# Finds and plots the optimal method combinations for a given Pareto front.
# Uses the following global variables: comb_list, gear_names, effort_per_sample, n_gears, n_combs, effort_max
### INPUTS:
##     pareto_front
##     gears: vector of gears to include (e.g., leave out eDNA if doing a traditional gear analysis)

stacked_area_plot = function(pareto_front, gears = c("eDNA","ef","fyke","gillnet","seine"), cols=cols) {
  
  # Calculate an effort-by-gear matrix in wide format (object effort_mat)
  effort_mat <- matrix(data = NA, nrow = n_combs, ncol = n_gears)
  for(i in 1:n_combs) {
    if(length(comb_list[[i]]) == 0) {
      effort_mat[i,] <- 0
    } else {
      effort_mat[i,] <- effort_per_sample * colSums(matrix(sapply(1:5, function(x) comb_list[[i]]==x), ncol = n_gears)) # calculate effort as effort_per_sample * number of samples (for each of the 5 gears)
    }
  }
  
  effort_df <- as.data.frame(effort_mat)
  colnames(effort_df) <- gear_names
  effort_df <- effort_df[,gears] # remove undesired gears
  effort_df$Tot_effort <- rowSums(as.matrix(effort_df)) # calculate total effort
  
  # Subset to gear combinations on the Pareto front
  effort_df <- effort_df[pareto_front[["indices"]],]
  names(effort_df)[names(effort_df) == "ef"] <- "Electrofishing" # rename "ef" as "electrofishing" (for legend)
  names(effort_df)[names(effort_df) == "fyke"] <- "Fyke" # rename "ef" as "electrofishing" (for legend)
  names(effort_df)[names(effort_df) == "gillnet"] <- "Gillnet" # rename "ef" as "electrofishing" (for legend)
  names(effort_df)[names(effort_df) == "seine"] <- "Seine" # rename "ef" as "electrofishing" (for legend)
  
  # Convert into long format (three columns: total effort, effort, and method)
  effort_mat_long <- melt(effort_df, id.vars = "Tot_effort")
  colnames(effort_mat_long) <- c("Tot_effort", "Gear", "Effort")
  
  # Make a stacked area plot
  ggplot(effort_mat_long, aes(x = Tot_effort, y = Effort, fill = Gear)) +
    geom_area() + 
    scale_fill_manual(values = cols) +
    xlab("Total effort (work hours)") + ylab("Optimal gear combination (work hours)")
  # TO DO: remove eDNA legend box from traditional samples plot.
  
}


# Stacked area plots for Pareto fronts
cols <- c(brewer.pal(8, "Dark2"),"#386CB0")
my.cols <- cols[c(4,9,1,6,3,2)]
stacked_area_plot(pareto_front = front_all, cols=my.cols)
stacked_area_plot(pareto_front = front_trad, cols=my.cols)

# Stacked area plots for convex hulls of Pareto fronts
# (optimal gear allocation vs. total effort is smoother, but still some noise for reasons described below)
stacked_area_plot(pareto_front = hull_all, cols=my.cols)
stacked_area_plot(pareto_front = hull_trad, gears = c("ef","fyke","gillnet","seine"), cols=my.cols)

## SAVED FIGURE -- for all gears, with convex hull of Pareto front ##
pdf(here("figures","stackedAreaPlot_convexPareto_allGears.pdf"), width = 6, height = 4)
stacked_area_plot(pareto_front = hull_all, cols=my.cols)
dev.off()




#### HEAT MAP OF GEAR BIASES ####
# Species x gear heat map of (1) the fraction of samples in which detection occurred, and (2) the average relative abundance.
#	Motivation: This heat map should help readers to visually see which species are not well sampled by eDNA, and lead naturally to a discussion of why this is the case for certain species (e.g., they don’t have a good reference). The "maximum detection fraction" is a proxy for commonness, and it makes sense to see that most species for which eDNA has a clear advantage over all traditional gears are relatively rare. The figure also may help provide intuition for why certain combinations of gears are effective but others are not.
# SOURCES: https://statisticsglobe.com/heatmap-in-r

detect_frac <- n/matrix(n_spls, nrow = n_spp, ncol = n_gears, byrow = TRUE)
max_detect_frac <- sapply(1:nrow(detect_frac), function(x) {max(detect_frac[x,])})
detect_frac_of_max <- detect_frac / max_detect_frac # normalized between 0 and 1 for each species, where 1 is achieved for the gear(s) attaining the max_detect_frac for that species

# Define function to replace common names with scientific names
replace_names <- function(i) {
  com_name <- gsub("[.]", " ", spp_names[i]); com_name <- gsub(" sp ", " sp.", com_name)
  if(com_name %in% spp_name_table$Common.name) {
    sci_name <- spp_name_table$Scientific.name[which(spp_name_table$Common.name == com_name)]
  } else {
    sci_name <- com_name
  }
  return(sci_name)
}

# Find scientific species names
spp_sci_names <- sapply(1:length(spp_names), function(x) replace_names(x))

# Define function to italicize character vectors
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))


## Plot and save to pdf
pdf(here("figures","gear_bias_heatmap.pdf"), width = 10, height = 8)

colMain <- colorRampPalette(brewer.pal(8, "Reds"))(50)
colSide <- colorRampPalette(brewer.pal(8, "Blues"))(50)
heatmap(detect_frac_of_max, 
        labCol = c("eDNA", "Electrofishing", "Fyke", "Gillnet", "Seine"),
        Colv = NA,
        distfun = function(x) {dist(x, method = "manhattan")},
        scale = "none",
        margins = c(7,21),
        cexRow = 0.9,  cexCol = 1.15,
        labRow = make.italic(spp_sci_names), # Alternatively: gsub("[.]", " ", spp_names)        
        col = colMain, # ... passed to image()
        RowSideColors = colSide[findInterval(max_detect_frac, vec = seq(0, 1, length.out = 51), rightmost.closed = TRUE)] # UNDER CONSTRUCTION!
        )

coords <- strip.legend(xy = "topright", 
                       legend = "    ",  #legend = seq(0,1,by=0.02), 
                       col = colSide, legendtype = c("breaks", "intervals", "other")[1],
             height = 0.25,
             title = "Maximum detection", text.cex = 0.7)
text(x = coords[2], y = coords[3:4]+c(-0.005,0), labels = c("1","0"), pos = 4, cex = 0.6)

coords <- strip.legend(xy = "bottomright", legend = "    ", col = colMain, legendtype = c("breaks", "intervals", "other")[1],
             height = 0.25,
             title = "Relative detection", text.cex = 0.7)
text(x = coords[2], y = coords[3:4]+c(-0.005,0), labels = c("Max", "0"), pos = 4, cex = 0.6)

dev.off()

# Note: Also tried using layout() to organize the plot and custom color bars / legends, but it seems it's not possible to implement when one of the plots is generated by heatmap() because heatmap() utilizes a conflicting call to layout().






#### NOTES ####
# For both single- and multi-method rarefaction curves, >20 samples per method is preferred (see Ch. 4: Estimating species richness by Nicholas J. Gotelli and Robert K. Colwell). We're on the lower end of that for electro-fishing (8 if we aggregate by site, 16 otherwise), seines (8) and gill nets (15).
# TO DO: Look into whether this code can be revised/improved/extended and submitted to R's vegan package, as an extension to vegan:specaccum.



#### UNDER CONSTRUCTION FROM HERE ON OUT ####


#### METHOD III: Continuous effort allocation #### UNDER CONSTRUCTION!!!
## Observation: Optimal method combinations on the Pareto front calculated using Method 2 can be variable across small ranges of total effort. This is due to the discreteness of sample combinations for methods with different effort---e.g., given Method A with sample cost 5 and Method B with sample cost 4, Method B may be Pareto-optimal at 8 and 12 but not 10 simply because its cost per sample allows more effort to be invested in it without exceeding the total effort limit, even if its cost efficiency is lower at all 3 values.
## Remedy: Method 3 eliminates this effect by shifting effort to a continuous scale for all methods. The expected number of species detected is interpolated between discrete numbers of samples for sampling WITH replacement (not without replacement, as in Method 2).
## Note: Taking the convex hull of the Pareto front is another possible remedy that does not require generalizations to continuous effort. See above.
## Additional benefit: Sampling with replacement allows extrapolation to larger sample sizes than we currently have available in the data.

# ## Function definition: expectedRichness() ##
# ## Calculates the expected number of species detected given inputs:
# #    par -- vector of length n_gears-1 indicating the fraction of total effort allocated to the first four gears
# ## And uses global parameters
# #    tot_effort -- the total effort expended
# #    A, n_gears, n_spp, n, effort_per_sample
# expectedRichness = function(par) {
#   b <- tot_effort * c(par, 1 - sum(par))
#   if(sum(b<0) > 0) {
#     return(-100)
#   } else {
#     A_mat <- matrix(A, ncol = n_gears, nrow = n_spp, byrow = TRUE)
#     nondet_frac <- (A_mat-n)/A_mat # (complement of detection_fraction in the heat map)
#     result <- sum(sapply(1:n_spp, function(i) 1 - prod(nondet_frac[i,]^(b/effort_per_sample))))
#     return(result)
#   }
# }
# 
# # USE OPTIM to find optimal allocation to each method, given a fixed total cost.
# n_grid_pts <- 11
# tot_effort_vec <- seq(0, effort_max, length.out = n_grid_pts)
# effort_mat <- matrix(NA, nrow = n_grid_pts, ncol = n_gears); colnames(effort_mat) <- gear_names
# x_cont <- y_cont <- rep(NA, length = n_grid_pts)
# for(i in 1:n_grid_pts) {
#   tot_effort <- tot_effort_vec[i]
#   result <- optim(rep(0.2, 4), expectedRichness,  # PROBLEM : different starts yielding different answers. Improve optim() or even use ga()?
#                   control = list(fnscale = -1), # maxit = 1000, abstol = 0.001
#                   method = "BFGS")
#   result2 <- optim(rep(0, 4), expectedRichness, 
#                    control = list(fnscale = -1), # maxit = 1000, abstol = 0.001
#                    method = "BFGS")
#   if(i>1 & sum((result$par - result2$par) > 0.05) > 0) {print(result, result2); break()}
#   x_cont[i] <- tot_effort
#   y_cont[i] <- result$value
#   effort_mat[i,] <- c(result$par, 1-sum(result$par)) * tot_effort
# }
# 
# plot(x_cont, y_cont)
# ## Goal: produce equivalents of x, y, and comb_list from Method 2
# ## Perhaps just: data frame with total effort, maximum number of species, and optimal allocation to each method (7 total columns).






# Side note: In some cases, gear choices are accumulative, i.e., the solution at one effort is always obtained by adding one sample (of a certain gear) to the optimal gear choices of the next lowest effort. What leads to cases where this property holds? What causes it to be violated? In cases where it holds, the stacked area plot may be simplified to a bar plot, with the y axis indicating the accrued gear choices as you increase from low to high effort.

