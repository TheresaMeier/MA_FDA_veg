################################################################################
############################ Master's Thesis ###################################
################################################################################

############## Exploratory Analysis: MFPCA - Clustering - Maps #################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/02_FPCA/functions.R")
source("Scripts/MA_FDA_veg/03_MFPCA/MFPCA/MFPCA_calculation.R")
source("Scripts/MA_FDA_veg/01_Description/utils.R")

## Load libraries
library(duckdb)
library(tidyverse)
library(ggplot2)
library(fda)
library(ggmap)
library(gifski)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
library(grid)
library(abind)
library(MFPCA)
library(foreach)
library(funData)

## Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)
scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")
pfts = c("Tundra", "BNE", "IBS", "otherC", "TeBS")
set.seed(1)

cluster_100 = read.table("Scripts/Plots/MFPCA/plot_data/plot_data_all.txt")

locs_picontrol = get_data_fpca("picontrol", 2015, 2040, 1, "Tundra")[[2]]
locs_ssp126 = get_data_fpca("ssp126", 2015, 2040, 1, "Tundra")[[2]]
locs_ssp370 = get_data_fpca("ssp370", 2015, 2040, 1, "Tundra")[[2]]
locs_ssp585 = get_data_fpca("ssp585", 2015, 2040, 1, "Tundra")[[2]]

locs_disturbed = rbind(locs_picontrol[,c("Lon", "Lat")], locs_ssp126[,c("Lon", "Lat")], locs_ssp370[,c("Lon", "Lat")], locs_ssp585[,c("Lon", "Lat")])

cluster_100 = cbind(cluster_100[,c("Cluster", "scenario")], locs_disturbed) %>% 
  mutate(Cluster = if_else(Cluster == 1, "Cluster 1", 
                 if_else(Cluster == 2, "Cluster 2",
                         if_else(Cluster == 3, "Cluster 3",
                                 if_else(Cluster == 4, "Cluster 4", NA))))) %>%
  mutate(Cluster = as.factor(Cluster))

########################## Visualize clustering ################################

# register_google(key = "AIzaSyATIWAZ4gtgJhdH6GP_E8iSubdFh6XQ32Y")
register_google(key = "AIzaSyA_eUpOhj7hoPyzyynWvyMqcGEA1Z_SZVY")

# Get map
worldmap <- get_map(location = c(lon =-4.068561, lat = 58.87355), zoom = 1)

ggmap(worldmap) +
  geom_point(data = cluster_100, 
             aes(x = Lon, y = Lat, fill = Cluster), 
             shape = 21, color = "black", size = 3, stroke = 0.5) + 
  facet_grid(rows = vars(scenario)) + 
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") +
  scale_color_manual(name = "Cluster", drop = TRUE, values = c("Cluster 1" = "#F8766D", "Cluster 2" = "#7CAE00" , "Cluster 3" = "#00BFC4", "Cluster 4" = "#C77CFF"),
                     labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")) +
  theme_bw() + ggtitle("Spatial distribution of the clusters") +
  theme(text = element_text(size = 15), plot.title = element_text(size = 20, face = "bold",hjust = 0.5))


ggsave("Scripts/Plots/MFPCA/Clusters/pdf/Map_clusters.pdf", width = 12, height = 8)
ggsave("Scripts/Plots/MFPCA/Clusters/png/Map_clusters.png", width = 12, height = 8)
