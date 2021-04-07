### January 29, 2021
### Make a map of eDNA sampling sites and electrofishing surveys in Oneida lake

# read in sampling sampling points longitude/latitude
eDNA_dat <- read.csv("/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/oneida_eDNA_sample_metadata.csv")
ef_dat <- read.csv("/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/EF_datasets/EF spring std site locations2.csv")
fyke_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/2017 oneida fyke.csv")
gillnet_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/2017 oneida gill net.csv")
seine_dat <- read.csv("/Users/kbja10/Github/Oneida_metabarcoding/datasets/2017 oneida seine.csv")

# Combine lat/longs for all sampling gears
eDNA_sampling_points <- data.frame(eDNA_dat[-12,c("lat","long")],gear=rep("eDNA",nrow(eDNA_dat)-1))
ef_sampling_points <- data.frame(lat=ef_dat$lat, long=ef_dat$long,gear=rep("ef",nrow(ef_dat)))
fyke_sampling_points <- data.frame(lat=fyke_dat$Latitude,long=fyke_dat$Longitude,gear=rep("fyke",nrow(fyke_dat)))
gillnet_sampling_points <- data.frame(lat=gillnet_dat$lat,long=gillnet_dat$long,gear=rep("gillnet",nrow(gillnet_dat)))
seine_sampling_points <- data.frame(lat=unique(na.omit(seine_dat$Latitude)),long=unique(na.omit(seine_dat$Longitude)),gear=rep("seine",length(unique(na.omit(seine_dat$Longitude)))))
all_sampling_points <- rbind(eDNA_sampling_points, ef_sampling_points, fyke_sampling_points,
                             gillnet_sampling_points, seine_sampling_points)

# make map
library(ggmap)
library(ggsn)
library(spData)
library(sf)
library(rgdal)
library(cowplot)
library(RColorBrewer)

lakes <- readOGR("/Users/kbja10/Downloads/ne_10m_lakes") # file path to lake shapefile
lakes <- lakes[complete.cases(lakes$name), ]
great_lakes_map <- subset(lakes, name==c("Lake Ontario","Lake Erie"))
oneida_map <- subset(lakes, name=="Oneida Lake")

mapbox <- c(-76.14, 43.14, -75.7, 43.3) # specify map boundaries
oneida_lake <- get_map(location = mapbox, source = "stamen", maptype = "terrain", zoom = 13)
cols <- c(brewer.pal(8, "Dark2"),"#386CB0")

oneida_lake_map <- ggmap(oneida_lake) +
  geom_point(data=all_sampling_points, mapping=aes(x=long, y=lat, color=gear, shape=gear), size=3, stroke = 1.2) +
  scale_color_manual(values=cols[c(4,9,6,3,1)]) +
  scale_shape_manual(values=c(17,5,3,16,6)) +
  geom_line(data=ef_dat, aes(x=long, y=lat, group=Description), color="white", size=0.5) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
      rect = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank()) +
  ggsn::scalebar(location = "bottomleft", x.min = -76.13, x.max = -75.7, 
           y.min = 43.145, y.max = 43.25, 
           dist = 5, dist_unit="km", transform = TRUE, 
           model = "WGS84", height = 0.05, 
           st.dist = 0.05, st.bottom=FALSE)

region_map <- ggplot() + geom_polygon(data=map_data("state"), aes(x=long, y=lat,group=group), color="white") + 
  geom_polygon(data=great_lakes_map, aes(x=long, y=lat, group=group), color="lightblue", fill="lightblue") + 
  geom_polygon(data=oneida_map, aes(x=long, y=lat, group=group), color="lightblue", fill="lightblue") + 
  coord_cartesian(xlim=c(-79, -70), ylim = c(40.05, 45.5)) +
  geom_rect(aes(xmin=-76.2, xmax=-75.67, ymin=43.05, ymax=43.4),fill = "transparent", color = "red", size = 0.5) +
  theme_map() + theme(panel.background = element_rect(fill = "white"))


inset_map <- ggdraw() + draw_plot(oneida_lake_map) +
  draw_plot(region_map, x = 0.63, y = 0.53, width = 0.26, height = 0.26)
inset_map

ggsave("/Users/kbja10/Documents/Cornell/Research/Oneida/Figures/oneida_lake_map.pdf", plot = inset_map, dpi=600)
