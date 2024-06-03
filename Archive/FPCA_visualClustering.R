################################################################################
############################ Master's Thesis ###################################
################################################################################

############ Exploratory Analysis: FPCA - Analysis of the Results ##############

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/02_FPCA/functions.R")
source("Scripts/MA_FDA_veg/_01_Description/utils.R")

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

## Choose parameters:
start_year = 2015
end_year = 2040
#pft = "TeBS"      # Choose from Tundra, BNE, IBS, otherC, TeBS
pid = 1           # Choose pid (int) or 'all'

## Load plot data generated in file 'FPCA_univ.R'

plot_data_Tundra = read.table(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/plot_data/plot_data_Tundra_",pid))
plot_data_BNE = read.table(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/plot_data/plot_data_BNE_",pid))
plot_data_IBS = read.table(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/plot_data/plot_data_IBS_",pid))
plot_data_otherC = read.table(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/plot_data/plot_data_otherC_",pid))
plot_data_TeBS = read.table(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/plot_data/plot_data_TeBS_",pid))

# Plot Tundra as before
ggplot(plot_data_Tundra) + 
  geom_point(data = plot_data_Tundra, 
             aes(x = PC1, y = PC2, color = region, group = interaction(Lon, Lat))) +
  facet_grid(rows = vars(name)) +
  scale_x_continuous(name = "PC 1") +
  scale_y_continuous(name = "PC 2") +
  scale_color_manual(name = "Region", drop = TRUE,values = c("Europe" = "slateblue4", "Asia" = "green2", "America" = "orange", "Other" = "gray")) +
  ggtitle(paste("Principal Component Scores for one Patch -", long_names_pfts(tolower("Tundra")))) +
  theme_bw() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5)) 

# Disentangle Clusters by visual inspection for Tundra
plot_data_Tundra = plot_data_Tundra %>%
  mutate(Cluster = as.factor(if_else(plot_data_Tundra$PC2 > 0.8 & plot_data_Tundra$PC1 < 3, "Cluster 1", 
                                             if_else(plot_data_Tundra$PC1 > 3, "Cluster 3", "Cluster 2"))),
         Cluster = if_else(plot_data_Tundra$name == "Control", if_else(plot_data_Tundra$PC2 < -1 & plot_data_Tundra$PC1 < 2.5, "Cluster 2", 
                                                                       if_else(plot_data_Tundra$PC1 > 3.1 & plot_data_Tundra$PC1 <= 6, "Cluster 3", 
                                                                               if_else(plot_data_Tundra$PC1 > 6, "Cluster 4", "Cluster 1")))  , Cluster),
         Cluster = if_else(plot_data_Tundra$name == "SSP3-RCP7.0",if_else(plot_data_Tundra$PC2 >-0.9, "Cluster 1", Cluster), Cluster)) %>%
  arrange(name,Lon,Lat)

# Add Clusters to the other PFTs
plot_data_BNE = plot_data_BNE %>% arrange(name,Lon,Lat) %>% mutate(Cluster = plot_data_Tundra$Cluster)
plot_data_IBS = plot_data_IBS %>% arrange(name,Lon,Lat) %>% mutate(Cluster = plot_data_Tundra$Cluster)
plot_data_otherC = plot_data_otherC %>% arrange(name,Lon,Lat) %>% mutate(Cluster = plot_data_Tundra$Cluster)
plot_data_TeBS = plot_data_TeBS %>% arrange(name,Lon,Lat) %>% mutate(Cluster = plot_data_Tundra$Cluster)


# Take Clusters in Tundra and visualize them in the other PFT plots

for (pft in (c("Tundra", "BNE", "IBS", "otherC", "TeBS"))){
  
  plot_data = get(paste0("plot_data_", pft))
  ggplot(plot_data) + 
    geom_point(data = plot_data, 
               aes(x = PC1, y = PC2, color = Cluster, group = interaction(Lon, Lat))) +
    facet_grid(rows = vars(name)) +
    scale_x_continuous(name = "PC 1") +
    scale_y_continuous(name = "PC 2") +
    ggtitle(paste("Principal Component Scores for one Patch -", long_names_pfts(tolower(pft)))) +
    theme_bw() +
    theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5)) 
  
  ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/unrotated/Clustering/visual/PC1_vs_PC2_cluster_",pft,"_",pid,".pdf"), width = 7, height = 4.5)
  ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/unrotated/Clustering/visual/PC1_vs_PC2_cluster_",pft,"_",pid,".png"), width = 7, height = 4.5)
  
}

## Plot the locations of the clusters
register_google(key = "AIzaSyATIWAZ4gtgJhdH6GP_E8iSubdFh6XQ32Y")

# Get map
worldmap <- get_map(location = c(lon =-4.068561, lat = 58.87355), zoom = 1)

ggmap(worldmap)+
  geom_point(data = plot_data, aes(x = Lon, y=Lat, color = Cluster), size = 0.8) +
  facet_grid(rows = vars(name)) +
  xlab("Longitude") +
  ylab("Latitude") +
  ylim(c(30,70)) +
  xlim(c(-180,180)) +
  theme_bw() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5)) +
  ggtitle("Spatial Distribution of the Clusters for Tundra")

ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/unrotated/Clustering/visual/PC1_vs_PC2_cluster_map.pdf"), width = 7, height = 5)
ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/unrotated/Clustering/visual/PC1_vs_PC2_cluster_map.png"), width = 7, height = 5)
