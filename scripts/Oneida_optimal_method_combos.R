### Optimal method combinations for surveying Oneida Lake fishes
### This code identifies the Pareto-optimal combination of methods (gears) to survey fish communities. Specifically, it resamples 2017 eDNA, electrofishing, fyke net, gill net, and seine net data from Oneida lake to estimate how to allocate effort among gears in a way that maximizes the mean number of species detected. The Pareto frontiers are calculated for combinations of all methods as well as for combinations of only traditional (non-eDNA) methods, and both are compared to species accumulation curves for each gear alone.
### Inputs: standardized 2017 gear files (eDNA_dat.csv, ef_dat.csv, fyke_dat.csv, gillnet_dat.csv, and seine_dat.csv).
### Outputs: (1) The Pareto frontiers for all gears and for all traditional gears, potted alongside species accumulation curves for each gear independently. (2) Stacked area plots indicating the optimal allocation to each survey type.
### Created by T. Lambert


#### Clear the working environment and load required packages ####

rm(list = ls())
library(here) # for setting paths
library(dplyr) # for use of bind_rows()
library(arrangements) # for use of combinations()
library(reshape) # for use of melt()



#### Read in and format data ####

# Read in the standardized separate gear files
eDNA_dat <- read.csv(here("datasets", "eDNA_dat.csv"))
ef_dat <- read.csv(here("datasets", "ef_dat.csv"))                  
fyke_dat <- read.csv(here("datasets", "fyke_dat.csv"))                  
gillnet_dat <- read.csv(here("datasets", "gillnet_dat.csv"))                  
seine_dat <- read.csv(here("datasets", "seine_dat.csv"))

# Restrict eDNA samples to the "green" (large-volume) samples
eDNA_dat <- subset(eDNA_dat, grepl("G", eDNA_dat$Site))

# Aggregate electrofishing samples by site, so each "sample" contains all 3 electrofishing runs per site (2 predator samples & 1 all-fish sample)
ef_dat <- aggregate(ef_dat[,-1], by = list(ef_dat$Site), FUN = "sum"); names(ef_dat)[1] <- "Site"

# Add a "Method" column
eDNA_dat$Method <- "eDNA"
ef_dat$Method <- "ef"
fyke_dat$Method <- "fyke"
gillnet_dat$Method <- "gillnet"
seine_dat$Method <- "seine"

# Combine all data into a single data frame
dat <- bind_rows(eDNA_dat, ef_dat, fyke_dat, gillnet_dat, seine_dat)
dat[is.na(dat)] <- 0 # count for species not detected by a method is 0


dat_by_method <- split(dat, dat$Method) # split dat data frame into a list by gear
rm(eDNA_dat, ef_dat, fyke_dat, gillnet_dat, seine_dat) # clean up the workspace



#### Define parameters ####

effort_per_sample <- c(eDNA = 5, ef = 5, fyke = 1.5, gillnet = 3, seine = 5) # for each gear, the effort required to obtain one sample (TO DO: OBTAIN REASONABLE ESTIMATES FOR THESE VALUES)

effort_max <- 40 # maximum effort

## Warnings:
# 1. Print a warning if the names (& order) of effort_per_sample and dat_by_method don't match.
if(sum(!names(dat_by_method) == names(effort_per_sample)) > 0) {
  print("Warning: gear names of effort_per_sample and dat_by_method do not match.")
}
# 2. Print a warning if the maximum effort that can be obtained by resampling the data without replacement falls short of effort_max for one or more gears. 
n_spls <- as.vector(table(dat$Method)) # number of samples per method
indicator <- (n_spls*effort_per_sample) < effort_max
if(sum(indicator) > 0) { print(paste("Warning: effort_max exceeds available data from which to sample for the following gear(s): ", paste(names(effort_per_sample)[indicator], collapse = ", "), ".", sep = ""))}
rm(indicator)


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







# #### METHOD #1: simulated sample combinations ####
# # NOTE: For testing purposes only. Use of Method #2 is preferred. 
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




#### Method 2: Exact calculation of expected number of species detected ####
## This is analogous to the argument <method = "exact"> in specaccum() from package vegan, except it has been extended to account for combinations of multiple methods.
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


# Calculate multiplication factor for 56 spp X 5 methods X 377 combs
f <- array(data = 1, dim = c(n_spp, n_gears, n_combs))
for(k in 1:n_combs) {
  for(i in 1:n_spp) {
    for(j in 1:n_gears) {
      f[i,j,k] <- ifelse(A[j] - n[i,j] < a[k,j], 0,
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


# Plots: the Pareto frontiers both for all gears and for all traditional (non-eDNA) gears are potted alongside species accumulation curves for each gear independently.
plot(front_all$front, xlim = c(0, effort_max), ylim = c(0, max(y)),
     pch = 16, col = "black", type = "b",
     xlab = "Effort", ylab = "Number of species detected",
     main = "Pareto frontiers for combined-method surveys")
points(front_eDNA$front, pch = 16, col = "red", type = "b")
points(front_ef$front, pch = 16, col = "lightgreen", type = "b")
points(front_fyke$front, pch = 16, col = "green", type = "b")
points(front_gillnet$front, pch = 16, col = "darkgreen", type = "b")
points(front_seine$front, pch = 16, col = "lightblue", type = "b")
points(front_trad$front, pch = 16, col = "purple", type = "b")

legend("topleft", legend = c("Combined eDNA and traditional", "Only eDNA", "Electrofishing", "Fyke", "Gillnet", "Seine", "Combined traditional"),
       pch = 16,
       col = c("black","red","lightgreen","green","darkgreen","lightblue","purple"), cex = 0.5)




# (2) Stacked area graph: for a given sampling effort, the optimal allocation to each survey type.
# https://r-graphics.org/recipe-line-graph-stacked-area





## Function definition: stacked_area_plot() ####
# Finds and plots the optimal method combinations for a given Pareto front.
# Uses the following global variables: comb_list, gear_names, effort_per_sample, n_gears, n_combs, effort_max
stacked_area_plot = function(pareto_front) {
  
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
  effort_df$Tot_effort <- rowSums(as.matrix(effort_df)) # calculate total effort
  
  # Subset to gear combinations on the Pareto front
  effort_df <- effort_df[pareto_front[["indices"]],]
  
  # Convert into long format (three columns: total effort, effort, and method)
  effort_mat_long <- melt(effort_df, id.vars = "Tot_effort")
  colnames(effort_mat_long) <- c("Tot_effort", "Gear", "Effort")
  
  # Make a stacked area plot
  ggplot(effort_mat_long, aes(x = Tot_effort, y = Effort, fill = Gear)) +
    geom_area()
  # TO DO: remove eDNA legend box from traditional samples plot.
  
}



stacked_area_plot(pareto_front = front_all)
stacked_area_plot(pareto_front = front_trad)






# TO DO: Make it clear that it's discrete, i.e., that Tot_effort increases in jumps, not continuously. One approach would be to draw in the points for each value of total effort that corresponds to a real sampling combination on the Pareto frontier.
# related TO DO: Make sure this is calculated by taking the largest number of species detected not only at that effort but also any lower effort, since in general an effort of exactly E may be impossible to achieve with integer numbers of samples, or only possible with some (possibly non-optimal) gear combinations.
# Try 1: Add a few lines after:
#     effort_df <- effort_df[pareto_front[["indices"]],]
# that are:
#     effort_df <- effort_df[rep(1:nrow(effort_df), each = 2), ]
#     effort_df$Tot_effort[-c(1,nrow(effort_df))] <- effort_df$Tot_effort[-c(1:2)]
#     effort_df$Tot_effort[nrow(effort_df)] <- effort_max
# But that doesn't have desired effect.


# OPTIONAL: Not always the case, but here it appears that cumulative gear choices are accumulative, i.e., the solution at one effort is always obtained by adding one sample (of a certain gear) to the optimal gear choices of the next lowest effort. Check that this still holds for the final choices of gear effort per sample. If so, then the stacked area plot may be simplified to a bar plot, with the y axis indicating the accrued gear choices as you increase from low to high effort.

# Methods -- outline of what to include in manuscript:
# Exhaustive list of all possible sample combinations up to some maximum effort. For each, perform min(# possible combinations given 2017 sample numbers, r) sampling simulations by sampling 2017 data without replacement [UPDATE: as currently implemented, I simply sample r replicates with replacement, even if r > # possible sample combinations], and calculate the average (and sd, if desired) number of species detected.
# 3. Define the Pareto frontier: the set of non-dominated solutions. Find the Pareto front for all traditional gears together, and for all gears including eDNA. Compare to the species accumulation curve for each gear independently (same as Kara's plot).









#### HEAT MAP OF GEAR BIASES ####
# Species x gear heat map of (1) the fraction of samples in which detection occurred, and (2) the average relative abundance.
#	Motivation: This heat map should help readers to visually see which species are not well sampled by eDNA, and lead naturally to a discussion of why this is the case for certain species (e.g., they don’t have a good reference). It also may help provide intuition for why certain combinations of gears are effective but others are not.
# SOURCES:
## https://statisticsglobe.com/heatmap-in-r


detect_frac <- n/matrix(n_spls, nrow = n_spp, ncol = n_gears, byrow = TRUE)

# pdf("Detection_fraction_heatmap.pdf", height = 11, width = 8.5)
heatmap(detect_frac, main = "Fraction of samples in which species was detected",
        cexRow = 0.6)
# dev.off()






# library("ggplot2")
# library(reshape)
# detect_frac_melt <- melt(detect_frac)
# ggp <- ggplot(detect_frac_melt, aes(X1, X2)) +
#   geom_tile(aes(fill = value))
# ggp
# ggp + scale_fill_gradient(low = "white", high = "black")





#### NOTES ####
# For both single- and multi-method rarefaction curves, >20 samples per method is preferred (see Ch. 4: Estimating species richness by Nicholas J. Gotelli and Robert K. Colwell). We're on the lower end of that for electro-fishing (8 if we aggregate by site), seines (8) and gill nets (15).
# TO DO: What happens if total effort is greater than possible effort given limited samples of a certain type? Is this case handled appropriately?

# TO DO (OPTIONAL): Change everything into continuous effort by fitting P(detection) ~ effort for each gear, on a species-by-species basis using the detection probability for a single sample of each gear with fixed effort.

# TO DO: Look into whether this code can be revised/improved/extended and submitted to R's vegan package, as an extension to vegan:specaccum.

