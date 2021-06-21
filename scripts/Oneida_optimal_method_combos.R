### Optimal method combinations for Oneida Lake
### This code identifies the Pareto-optimal combination of methods (gears) to survey fish communities. More specifically, it resamples 2017 eDNA, electrofishing, fyke net, gill net, and seine net data from Oneida lake to estimate how to allocate effort among gears in a way that maximizes the mean number of species detected. The Pareto frontiers are calculated for combinations of all methods as well as for combinations of only traditional (non-eDNA) methods, and both are compared to species accumulation curves for each gear alone.
### Inputs: standardized 2017 gear files (eDNA_dat.csv, ef_dat.csv, fyke_dat.csv, gillnet_dat.csv, and seine_dat.csv).
### Outputs: (1) Plots: the Pareto frontiers both for all gears and for all traditional (non-eDNA) gears are potted alongside species accumulation curves for each gear independently. (2) Tables: for a given sampling effort, the optimal allocation to each survey type.
### Created by T. Lambert


#### Clear the working environment and load required packages ####
rm(list = ls())
library(here) # for setting paths
library(dplyr) # for use of bind_rows()
library(arrangements) # for use of combinations()


#### Read in and format data ####

# Read in the standardized separate gear files
eDNA_dat <- read.csv(here("datasets", "eDNA_dat.csv"))
ef_dat <- read.csv(here("datasets", "ef_dat.csv"))                  
fyke_dat <- read.csv(here("datasets", "fyke_dat.csv"))                  
gillnet_dat <- read.csv(here("datasets", "gillnet_dat.csv"))                  
seine_dat <- read.csv(here("datasets", "seine_dat.csv"))

# Restrict to the "green" (large-volume) samples
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


#### Define parameters ####

effort_per_sample <- c(eDNA = 2, ef = 4, fyke = 4, gillnet = 4, seine = 4) # for each gear, the effort required to obtain one sample (TO DO: OBTAIN REASONABLE ESTIMATES FOR THESE VALUES)

r <- 60 # maximum number of replicate simulations per sample combination (FEWER ONLY IF FEWER THAN r SIMULATIONS ARE POSSIBLE; NOT YET ENFORCED -- DO THIS)

effort_max <- 20 # maximum effort (TO DO: enforce effort_max is below the available sampling effort for each gear that can be sampled without replacement; if not, we are artificially restricting the combinations possible at higher gear)


## Print a warning if the maximum effort that can be obtained by resampling the data without replacement falls short of effort_max for one or more gears. 
n_spls <- as.vector(table(dat$Method)) # number of samples per method
indicator <- (n_spls*effort_per_sample) < effort_max
if(sum(indicator) > 0) { print(paste("Warning: effort_max exceeds available data from which to sample for the following gear(s): ", paste(names(effort_per_sample)[indicator], collapse = ", "), ".", sep = ""))}
rm(indicator)


#### Create list of all possible sample combinations with effort < effort_max ####
k_max <- max(floor(effort_max / effort_per_sample)) # maximum number of samples to select per survey/combination (auto-generated so that resampling achieves effort_max for every method and all their possible combinations)
comb_list <- list() # initialize cumulative list of all sample combinations
for(k in 1:k_max) { # loop through different values of k (# samples selected per combination)
  combs <- combinations(freq = floor(effort_max/effort_per_sample), k = k, x = 1:5)
  effort_spls <- sapply(1:nrow(combs), function(x) sum(effort_per_sample[combs[x,]])) # identify total effort per sample
  combs <- subset(combs, effort_spls <= effort_max) # cap total effort per sample
  # TO DO: ALSO CAP AT TOTAL NUMBER OF SAMPLES PER METHOD (e.g., >9 seine samples in a combination is not possible given sampling without replacement)
  comb_list <- c(comb_list, as.list(data.frame(t(combs)))) # append new combs to cumulative sample list
}
rm(combs)
## TO DO: RENAMING -- perhaps rename as COMBS for combinations rather than sample (e.g., comb_list not comb_list) to help distinguish preliminary from final objects used below.


#### Simulated sample combinations ####
# Perform r simulations per sampling combination and calculate the mean & stdev detected species richness

n_spp_detected <- matrix(data = NA, nrow = length(comb_list), ncol = r) # initialize matrix to contain the numbers of species detected for each sample combination (row) and replicate (column)
effort_combs <- sapply(1:length(comb_list), function(x) sum(effort_per_sample[comb_list[[x]]])) # recalculate total effort per sample -- TO STREAMLINE: if costly, could find a way to use calculations above and subset, rather than recalculating.

# USEFUL OBJECTS: effort_combs, comb_list, n_spls, dat, dat_by_method
dat_by_method <- list(subset(dat, dat$Method == "eDNA"),
                      subset(dat, dat$Method == "ef"),
                      subset(dat, dat$Method == "fyke"),
                      subset(dat, dat$Method == "gillnet"),
                      subset(dat, dat$Method == "seine"))

for(i in 1:length(comb_list)) { # for each sampling combination...
  n_to_select <- sapply(1:5, function(x) sum(comb_list[[i]]==x)) # number of samples per method/gear
  for(j in 1:r) { # for each of r replicates...

    for(m in 1:5) { # for each method/gear
      # Randomly select the appropriate number of samples for each gear
      if(m==1) {
        dat_temp <- dat_by_method[[m]][sample(n_spls[m], n_to_select[m]), ] # initialize data frame containing simulated data for the current replicate
      } else {
        dat_temp <- rbind(dat_temp, dat_by_method[[m]][sample(n_spls[m], n_to_select[m]), ])
      }
    }
    dat_temp <- dat_temp[ ,!colnames(dat_temp) %in% c("Method","Site")]
    n_spp_detected[i,j] <- sum(colSums(dat_temp) > 0) # calculate the number of species detected
    
  }
}



#### Find the Pareto frontiers ####

## Function definition: pareto()
## Finds the Pareto front (set of non-dominated solutions).
## Output: a data frame of x- and y-values for the Pareto front, with effort ordered from low to high.
## Inputs: x: x-values (effort)
##   y: y-values (# spp. detected, or more generally an objective to be maximized)
##   condition: integer or logical vector indicating the indices of observations to be kept
pareto = function(x, y, condition = NA) {
  d <- data.frame(x,y)
  if(!is.na(condition)[1]) {
    d <- d[condition, ]
  }
  D <- d[order(d$x,d$y,decreasing=c(FALSE,TRUE)),]
  front <- D[which(!duplicated(cummax(D$y))),]
  return(front)
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




x <- effort_combs # effort for each sample combination
y <- rowSums(n_spp_detected)/r # average number of species detected

front_all <- pareto(x, y)
front_trad <- pareto(x, y, condition = idCombs(keep = 2:5))
front_eDNA <- pareto(x, y, condition = idCombs(keep = 1))
front_ef <- pareto(x, y, condition = idCombs(keep = 2))
front_fyke <- pareto(x, y, condition = idCombs(keep = 3))
front_gillnet <- pareto(x, y, condition = idCombs(keep = 4))
front_seine <- pareto(x, y, condition = idCombs(keep = 5))






# (1) Plots: the Pareto frontiers both for all gears and for all traditional (non-eDNA) gears are potted alongside species accumulation curves for each gear independently.
plot(front_all, xlim = c(0, effort_max), ylim = c(0, max(y)),
     pch = 16, col = "black", type = "b",
     xlab = "Effort", ylab = "Number of species detected")
points(front_eDNA, pch = 16, col = "red", type = "b")
points(front_ef, pch = 16, col = "lightgreen", type = "b")
points(front_fyke, pch = 16, col = "green", type = "b")
points(front_gillnet, pch = 16, col = "darkgreen", type = "b")
points(front_seine, pch = 16, col = "lightblue", type = "b")
points(front_trad, pch = 16, col = "purple", type = "b")


### PICK UP HERE... UNFINISHED BEYOND THIS!!!!
legend("topleft", legend = c("Combined eDNA and traditional", "Only eDNA", "Only traditional"),
       pch = 16,
       col = c("black","red","green"), cex = 0.8)



# TO DO: compare each method alone to Kara's species accumulation curves (specaccum from package vegan) -- they should match up.
# Storyline:





# (2) Tables: for a given sampling effort, the optimal allocation to each survey type.
# TO DO: Make sure this is calculated by taking the largest number of species detected not only at that effort but also any lower effort, since in general an effort of exactly E may be impossible to achieve with integer numbers of samples, or only possible with some (possibly non-optimal) gear combinations.



# Methods -- outline of what to include in manuscript:
# Exhaustive list of all possible sample combinations up to some maximum effort. For each, perform min(# possible combinations given 2017 sample numbers, r) sampling simulations by sampling 2017 data without replacement [UPDATE: as currently implemented, I simply sample r replicates with replacement, even if r > # possible sample combinations], and calculate the average (and sd, if desired) number of species detected.
# 3. Define the Pareto frontier: the set of non-dominated solutions. Find the Pareto front for all traditional gears together, and for all gears including eDNA. Compare to the species accumulation curve for each gear independently (same as Kara's plot).





#### Method 2: Exact calculation of expected number of species detected ####
## This is analogous to the argument <method = "exact"> in specaccum() from package vegan, except it has been extended to account for combinations of multiple methods.
# Inputs:
## comb_list
## dat, dat_by_method
## 
spp_names <- colnames(dat)[!colnames(dat) %in% c("Method","Site")]
gear_names <- names(effort_per_sample)
n_spp <- length(species_names)
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

res <- rep(NA, n_combs) # calculate the expected number of species for each sample combination
for(k in 1:n_combs) {
  res[k] <- sum(1 - f[,1,k]*f[,2,k]*f[,3,k]*f[,4,k]*f[,5,k])
}


## PLOTTING THE EXACT METHOD -- duplicated code from above... ##

x <- effort_combs # effort for each sample combination
y <- res # expected number of species detected (exact calculation)

front_all <- pareto(x, y)
front_trad <- pareto(x, y, condition = idCombs(keep = 2:5))
front_eDNA <- pareto(x, y, condition = idCombs(keep = 1))
front_ef <- pareto(x, y, condition = idCombs(keep = 2))
front_fyke <- pareto(x, y, condition = idCombs(keep = 3))
front_gillnet <- pareto(x, y, condition = idCombs(keep = 4))
front_seine <- pareto(x, y, condition = idCombs(keep = 5))


# Plots: the Pareto frontiers both for all gears and for all traditional (non-eDNA) gears are potted alongside species accumulation curves for each gear independently.
plot(front_all, xlim = c(0, effort_max), ylim = c(0, max(y)),
     pch = 16, col = "black", type = "b",
     xlab = "Effort", ylab = "Number of species detected",
     main = "Pareto frontiers for combined-method surveys")
points(front_eDNA, pch = 16, col = "red", type = "b")
points(front_ef, pch = 16, col = "lightgreen", type = "b")
points(front_fyke, pch = 16, col = "green", type = "b")
points(front_gillnet, pch = 16, col = "darkgreen", type = "b")
points(front_seine, pch = 16, col = "lightblue", type = "b")
points(front_trad, pch = 16, col = "purple", type = "b")

legend("topleft", legend = c("Combined eDNA and traditional", "Only eDNA", "Electrofishing", "Fyke", "Gillnet", "Seine", "Combined traditional"),
       pch = 16,
       col = c("black","red","lightgreen","green","darkgreen","lightblue","purple"), cex = 0.5)


### Exact method looks to be correct! Now identify the combinations that are on the Pareto front and characterize their composition (e.g., eDNA + a bit of seine is the best possible strategy.)


# TO DO: Look into whether this code can be revised/improved/extended and submitted to R's vegan package, as an extension to vegan:specaccum.



#### HEAT MAP OF GEAR BIASES ####
# Species x gear heat map of (1) the fraction of samples in which detection occurred, and (2) the average relative abundance.
#	Motivation: This heat map should help readers to visually see which species are not well sampled by eDNA, and lead naturally to a discussion of why this is the case for certain species (e.g., they don’t have a good reference). It also may help provide intuition for why certain combinations of gears are effective but others are not.
# SOURCES:
## https://statisticsglobe.com/heatmap-in-r


detect_frac <- n/matrix(n_spls, nrow = n_spp, ncol = n_gears, byrow = TRUE)

pdf("Detection_fraction_heatmap.pdf", height = 11, width = 8.5)
heatmap(detect_frac, main = "Fraction of samples in which species was detected",
        cexRow = 0.6)
dev.off()


library("ggplot2")
library(reshape)
detect_frac_melt <- melt(detect_frac)
ggp <- ggplot(detect_frac_melt, aes(X1, X2)) +
  geom_tile(aes(fill = value))
ggp
ggp + scale_fill_gradient(low = "white", high = "black")





#### NOTES ####
# For both single- and multi-method rarefaction curves, >20 samples per method is preferred (see Ch. 4: Estimating species richness by Nicholas J. Gotelli and Robert K. Colwell). We're on the lower end of that for electro-fishing (8 if we aggregate by site), seines (8) and gill nets (15).
# TO DO: What happens if total effort is greater than possible effort given limited samples of a certain type? Is this case handled appropriately?
