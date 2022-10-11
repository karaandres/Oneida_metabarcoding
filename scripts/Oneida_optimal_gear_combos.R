### Optimal gear combinations for surveying Oneida Lake fishes
### Created by T. Lambert
### This code identifies the Pareto-optimal combination of methods (gears) to survey fish communities. 
### Specifically, it resamples 2017 eDNA, electrofishing, fyke net, gillnet, and seining data from Oneida Lake to estimate how to allocate effort among gears in a way that maximizes the mean number of species detected. 
### The Pareto frontiers are calculated for combinations of all methods as well as for combinations of only traditional (non-eDNA) methods, and both are compared to species accumulation curves for each gear alone.
### Inputs:
###   Standardized 2017 gear files (eDNA_dat.csv, ef_dat.csv, fyke_dat.csv, gillnet_dat.csv, and seine_dat.csv).
### Outputs:
###   (1) The Pareto frontiers for all gears and for all traditional gears, potted alongside species accumulation curves for each gear independently.
###   (2) Stacked area plots indicating the optimal allocation to each survey type, using two alternative approaches: (i) taking discrete samples via resampling without replacement, and (ii) [not yet implemented] an extension to continuous effort by interpolation of resampling with replacement.
###   (3) Heatmap of gear biases, in terms of detection fraction for each species.

#### Clear the working environment and load required packages ####
rm(list = ls())
library(here) # for setting paths
library(plyr) # for use of ddply()
library(dplyr) # for use of bind_rows()
library(arrangements) # for use of combinations()
library(reshape2) # for use of melt()
library(secr) # for use of strip.legend()
library(RColorBrewer) # for color palettes
library(data.table) # for subsetting data

#### Read in and format data ####

# Read in the standardized separate gear files
eDNA_dat <- read.csv("datasets/eDNA_dat.csv")
ef_dat <- read.csv("datasets/ef_dat.csv")                  
fyke_dat <- read.csv("datasets/fyke_dat.csv")                
gillnet_dat <- read.csv("datasets/gillnet_dat.csv")                
seine_dat <- read.csv("datasets/seine_dat.csv")

# Make modifications to gear files as necessary
eDNA_dat <- subset(eDNA_dat, grepl("G", eDNA_dat$Site)) # Restrict eDNA samples to the "green" (large-volume) samples
ef_dat <- ef_dat[((1:nrow(ef_dat)) %% 3) %in% c(1,0), ] # NOTE: This requires that electrofishing samples be grouped by site, each with 3 samples in the following order: (1) all-fish run, (2) predator run, (3) all-fish run.
ef_dat <- ddply(ef_dat, "Site", numcolwise(sum)) # Combine replicates per site (3)              
fyke_dat <- ddply(fyke_dat, "Site", numcolwise(sum)) # Combine replicates per site (2)

# Add a "Method" column
eDNA_dat$Method <- "eDNA"
ef_dat$Method <- "ef"
fyke_dat$Method <- "fyke"
gillnet_dat$Method <- "gillnet"
seine_dat$Method <- "seine"

# Read in species names (used as a common-to-scientific name look-up table)
all_methods_presence <- read.csv("datasets/all_methods_presence.csv")
spp_name_table <- all_methods_presence[,c("Scientific.name","Common.name")]
rm(all_methods_presence)

# Combine all data into a single data frame
dat <- bind_rows(eDNA_dat, ef_dat, fyke_dat, gillnet_dat, seine_dat)
dat[is.na(dat)] <- 0 # count for species not detected by a method is 0

# Restrict to species-level assignments
dat <- dat[,colnames(dat)[!grepl(pattern = "[.]sp[.]", colnames(dat))]]
dat_by_method <- split(dat, dat$Method) # split dat data frame into a list by gear
rm(eDNA_dat, ef_dat, fyke_dat, gillnet_dat, seine_dat) # clean up the workspace

#### Define parameters ####

# EFfort estimates
effort_per_sample <- c(eDNA = 2.0, ef = 4, fyke = 3.6, gillnet = 13, seine = 2.5) # for each gear, the approximate effort (in person-hours of work) required to obtain one sample
effort_max <- 100 # maximum effort

# Low eDNA cost estimates (cost/10 for ease of calculation)
effort_per_sample <- c(eDNA = 15.2, ef = 28, fyke = 9.2, gillnet = 32, seine = 7) # for each gear, the approximate effort (in person-hours of work) required to obtain one sample
effort_max <- 400 # maximum cost

# High cost estimates (cost/10 for ease of calculation)
effort_per_sample <- c(eDNA = 28.4, ef = 28, fyke = 9.2, gillnet = 32, seine = 7) # for each gear, the approximate effort (in person-hours of work) required to obtain one sample
effort_max <- 400 # maximum cost

#### Create list of all possible sample combinations ####
## (with the condition that total effort must be less than effort_max) ##
k_max <- max(floor(effort_max / effort_per_sample)) # maximum number of samples to select per survey/combination (auto-generated so that resampling achieves effort_max for every method and all their possible combinations)
comb_list <- list(as.integer(NULL)) # initialize cumulative list of all sample combinations; first entry is the combination with no samples of any sort
for(k in 1:k_max) { # loop through different values of k (# samples selected per combination)
  combs <- combinations(freq = floor(effort_max/effort_per_sample), k = k, x = 1:5)
  effort_spls <- sapply(1:nrow(combs), function(x) sum(effort_per_sample[combs[x,]])) # identify total effort per sample
  combs <- subset(combs, effort_spls <= effort_max) # cap total effort per sample
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

#### Exact calculation of expected number of species detected ####
## This is analogous to the argument <method = "exact"> in specaccum() from package vegan, except it has been extended to account for combinations of multiple methods. It uses resampling without replacement.
spp_names <- colnames(dat)[!colnames(dat) %in% c("Method","Site")]
gear_names <- names(effort_per_sample)
n_spp <- length(spp_names)
n_gears <- length(unique(dat$Method))
n_combs <- length(comb_list)
n_spls <- as.vector(table(dat$Method)) # number of samples per method
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
f <- array(data=1, dim=c(n_spp, n_gears, n_combs))
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

#### CALCULATE AND PLOT THE CONVEX HULLS OF THE PARETO FRONTIERS ####
## Enforce concavity of the Pareto frontier by eliminating points that are dominated by linear combination of other points. 
hull_all <- conv_hull(front_all)
hull_trad <- conv_hull(front_trad)

cols <- c(brewer.pal(8, "Dark2"),"#386CB0")
my.cols <- cols[c(8,4,9,1,6,3,2)]
# pdf("figures/convex_Pareto_combination_high_cost.pdf", width = 8, height = 6)
plot(hull_all$hull, xlim = c(0, effort_max), ylim = c(0, max(y)),
     pch = 1, col = my.cols[1], type = "b", lwd = 1.5,
     xlab = "Effort (work hours)", 
     ylab = "Expected number of species detected")
points(front_eDNA$front, pch = 16, col = my.cols[2], lwd = 4, type = "o")
points(front_ef$front, pch = 16, col = my.cols[3], lwd = 4, type = "o")
points(front_seine$front, pch = 16, col = my.cols[4], lwd = 4, type = "o")
points(front_fyke$front, pch = 16, col = my.cols[5], lwd = 4, type = "o")
points(front_gillnet$front, pch = 16, col = my.cols[6], lwd = 4, type = "o")
points(hull_trad$hull, pch = 1, col = my.cols[7], lwd = 1.5, type = "b")
legend("topleft", legend = c("eDNA", "Electrofishing", "Seining", "Fyke netting", "Gillnetting", "All capture gears", "All gears"),
       pch = c(rep(16,5),1,1),
       col = my.cols[c(2:7,1)], cex = 1)
# dev.off()

#### (2) PLOT 2: Stacked bar graph -- optimal gear allocation ####
# Discrete samples: for a given sampling effort, the optimal allocation to each survey type.
library(ggplot2)

## Function definition: stacked_area_plot() ####
stacked_bar_plot = function(pareto_front, cols=cols) {
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
  # Rename gears for legend
  supp.labs <- c("eDNA", "Electrofishing", "Fyke netting", "Gillnetting", "Seining")
  names(supp.labs) <- gear_names
  # Convert into long format (three columns: total effort, effort, and method)
  effort_mat_long <- reshape2::melt(effort_df, id.vars = "Tot_effort")
  colnames(effort_mat_long) <- c("Tot_effort", "Gear", "Effort")
  # Subset to effort in multiples of 3 (to make figure look smoother)
  effort_mat_long_sub <- data.frame(Tot_effort=NULL, Gear=NULL, Effort=NULL)
  for (i in seq(10,effort_max,by=effort_max/10)){
    temp <- effort_mat_long[which(effort_mat_long$Tot_effort>i-(effort_max/10/2) & effort_mat_long$Tot_effort<=i+(effort_max/10/2)),]
    if (i %in% temp$Tot_effort) { # if exact amount of effort is estimated
      temp2 <- temp[temp$Tot_effort==i,]
      effort_mat_long_sub <- rbind(effort_mat_long_sub, temp2)
    } else {
      temp3 <- as.data.frame(temp %>% group_by(Gear) %>% summarise_at("Effort", mean, na.rm=TRUE))
      temp4 <- data.frame(Tot_effort=i, temp3)
      effort_mat_long_sub <- rbind(effort_mat_long_sub, temp4)
    }
  }
  # Make a stacked area plot
  ggplot(effort_mat_long_sub, aes(x=Tot_effort, y=Effort, fill=Gear)) +
    geom_bar(position="stack", stat="identity", width=effort_max/10/2) + 
    scale_fill_manual(values=cols, labels=supp.labs) +
    xlab("Total effort (work hours)") + ylab("Optimal gear combination (effort)") +
    scale_x_continuous(breaks=seq(10,effort_max,by=effort_max/10)) +
    scale_y_continuous(breaks=seq(10,effort_max,by=effort_max/10)) +
    theme_bw()
}

# Stacked area plots for Pareto fronts
cols <- c(brewer.pal(8, "Dark2"),"#386CB0")
my.cols <- cols[c(4,9,6,3,1)]
# pdf("figures/stackedAreaPlot_convexPareto_allGears_low_cost.pdf", width = 6, height = 4)
stacked_bar_plot(pareto_front = hull_all, cols=my.cols)
# dev.off()

effort_hull <- hull_all
low_cost_hull <- hull_all
high_cost_hull <- hull_all

### Donut chart (for graphical abstract): effort
data <- data.frame(category=c("A", "B", "C"), count=c(50, 25, 22.5))
data$fraction <- data$count / sum(data$count)
data$ymax <- cumsum(data$fraction)
data$ymin <- c(0, head(data$ymax, n=-1))
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  scale_fill_manual(values=my.cols[c(1,3,5)]) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

# Donut chart (for graphical abstract): external cost
data <- data.frame(category=c("A", "B", "C"), count=c(284, 27,63))
data$fraction <- data$count / sum(data$count)
data$ymax <- cumsum(data$fraction)
data$ymin <- c(0, head(data$ymax, n=-1))
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  scale_fill_manual(values=my.cols[c(1,3,5)]) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

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
# pdf(here("figures","gear_bias_heatmap.pdf"), width = 10, height = 8)
colMain <- colorRampPalette(brewer.pal(8, "Purples"))(50)
colSide <- colorRampPalette(brewer.pal(8, "Blues"))(50)
heatmap(detect_frac_of_max, 
        labCol = c("eDNA", "Electrofishing", "Fyke netting", "Gillnetting", "Seining"),
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
# dev.off()

#### NOTES ####
# For both single- and multi-method rarefaction curves, >20 samples per method is preferred (see Ch. 4: Estimating species richness by Nicholas J. Gotelli and Robert K. Colwell). We're on the lower end of that for electro-fishing (8 if we aggregate by site, 16 otherwise), seines (8) and gillnets (15).
# TO DO: Look into whether this code can be revised/improved/extended and submitted to R's vegan package, as an extension to vegan:specaccum.