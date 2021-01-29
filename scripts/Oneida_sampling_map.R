### January 29, 2021
### Make a map of eDNA sampling sites and electrofishing surveys in Oneida lake

# read in sampling sampling points longitude/latitude
sampling_points <- read.csv("/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/oneida_eDNA_sample_metadata.csv") # not projected
sampling_points <- sampling_points[-12,c("lat","long")] # lat/long, remove field blank (sample 12)

# make map
library(ggmap)
library(ggsn)
mapbox <- c(-76.2, 43.13, -75.7, 43.28) # specify map boundaries
oneida_lake <- get_map(location = mapbox, source = "stamen", maptype = "terrain", zoom = 12)

oneida_lake_map <- ggmap(oneida_lake) +
  geom_point(data = sampling_points,aes(x = long, y = lat), color = "black", size = 3, alpha = 0.7) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        rect = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank()) +
  ggsn::scalebar(location = "bottomleft", x.min = -76.18, x.max = -75.7, 
           y.min = 43.138, y.max = 43.28, 
           dist = 5, dist_unit="km", transform = TRUE, 
           model = "WGS84", height = 0.05, 
           st.dist = 0.05, st.bottom=FALSE)
# ggsave("/Users/kbja10/Documents/Cornell/Research/Oneida/Figures/oneida_lake_map.png", plot = oneida_lake_map)

