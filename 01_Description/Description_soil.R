################################################################################
############################ Master's Thesis ###################################
################################################################################

########################### Climate and Soil data ##############################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/01_Description/utils.R")

## Load libraries
library(dplyr)
library(data.table)
library(stringr)
library(cowplot)
library(ggplot2)
library(ggmap)

scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")

## Load data
pid = 1

for (scen in scenarios){
  d_scen_2015_2040 = fread(paste0("Data/data_", scen, "_2015_2040.csv"))
  
  d_soil_scen = d_scen_2015_2040 %>%
    distinct(Lon,Lat,sand_fraction,silt_fraction, clay_fraction, bulkdensity_soil, ph_soil, soilcarbon)
  
  assign(paste0("d_soil_", scen), d_soil_scen)
  
}

register_google(key = "AIzaSyATIWAZ4gtgJhdH6GP_E8iSubdFh6XQ32Y")

# Get map
worldmap <- get_map(location = c(lon =-4.068561, lat = 58.87355), zoom = 1)

g_picontrol_sand = ggmap(worldmap) +
  geom_point(data = d_soil_picontrol, 
             aes(x = Lon, y = Lat, color = sand_fraction), 
             size = 1.5) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("Control") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(color = "Sand Fraction") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  #scale_color_gradient(low = "slateblue4", high = "firebrick1")  
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 0.5, limits = c(0, 1))  

g_ssp126_sand = ggmap(worldmap) +
  geom_point(data = d_soil_ssp126, 
             aes(x = Lon, y = Lat, color = sand_fraction), 
             size = 1.5) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("SSP1-RCP2.6") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(color = "Sand Fraction") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  #scale_color_gradient(low = "slateblue4", high = "firebrick1")  
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 0.5, limits = c(0, 1))  

g_ssp370_sand = ggmap(worldmap) +
  geom_point(data = d_soil_ssp370, 
             aes(x = Lon, y = Lat, color = sand_fraction), 
             size = 1.5) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("SSP3-RCP7.0") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(color = "Sand Fraction") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  #scale_color_gradient(low = "slateblue4", high = "firebrick1")  
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 0.5, limits = c(0, 1))  

g_ssp585_sand = ggmap(worldmap) +
  geom_point(data = d_soil_ssp585, 
             aes(x = Lon, y = Lat, color = sand_fraction), 
             size = 1.5) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("SSP5-RCP8.5") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(color = "Sand Fraction") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  #scale_color_gradient(low = "slateblue4", high = "firebrick1")  
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 0.5, limits = c(0, 1))  

plot_grid(g_picontrol_sand, g_ssp126_sand, g_ssp370_sand, g_ssp585_sand, ncol = 1)

ggsave("Scripts/Plots/Descriptive/Soil/pdf/sand.pdf", width = 20, height = 15)
ggsave("Scripts/Plots/Descriptive/Soil/png/sand.png", width = 20, height = 15)

g_picontrol_clay = ggmap(worldmap) +
  geom_point(data = d_soil_picontrol, 
             aes(x = Lon, y = Lat, color = clay_fraction), 
             size = 1.5) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("Control") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(color = "Clay Fraction") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  #scale_color_gradient(low = "slateblue4", high = "firebrick1")  
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 0.5, limits = c(0, 1))  

g_ssp126_clay = ggmap(worldmap) +
  geom_point(data = d_soil_ssp126, 
             aes(x = Lon, y = Lat, color = clay_fraction), 
             size = 1.5) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("SSP1-RCP2.6") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(color = "Clay Fraction") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  #scale_color_gradient(low = "slateblue4", high = "firebrick1")  
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 0.5, limits = c(0, 1))  

g_ssp370_clay = ggmap(worldmap) +
  geom_point(data = d_soil_ssp370, 
             aes(x = Lon, y = Lat, color = clay_fraction), 
             size = 1.5) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("SSP3-RCP7.0") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(color = "Clay Fraction") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  #scale_color_gradient(low = "slateblue4", high = "firebrick1")  
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 0.5, limits = c(0, 1))  

g_ssp585_clay = ggmap(worldmap) +
  geom_point(data = d_soil_ssp585, 
             aes(x = Lon, y = Lat, color = clay_fraction), 
             size = 1.5) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("SSP5-RCP8.5") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(color = "Clay Fraction") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  #scale_color_gradient(low = "slateblue4", high = "firebrick1")  
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 0.5, limits = c(0, 1))  


plot_grid(g_picontrol_clay, g_ssp126_clay, g_ssp370_clay, g_ssp585_clay, ncol = 1)
ggsave("Scripts/Plots/Descriptive/Soil/pdf/clay.pdf", width = 20, height = 15)
ggsave("Scripts/Plots/Descriptive/Soil/png/clay.png", width = 20, height = 15)


g_picontrol_silt = ggmap(worldmap) +
  geom_point(data = d_soil_picontrol, 
             aes(x = Lon, y = Lat, color = silt_fraction), 
             size = 1.5) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("Control") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(color = "Silt Fraction") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  #scale_color_gradient(low = "slateblue4", high = "firebrick1")  
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 0.5, limits = c(0, 1))  

g_ssp126_silt = ggmap(worldmap) +
  geom_point(data = d_soil_ssp126, 
             aes(x = Lon, y = Lat, color = silt_fraction), 
             size = 1.5) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("SSP1-RCP2.6") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(color = "Silt Fraction") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  #scale_color_gradient(low = "slateblue4", high = "firebrick1")  
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 0.5, limits = c(0, 1))  


g_ssp370_silt = ggmap(worldmap) +
  geom_point(data = d_soil_ssp370, 
             aes(x = Lon, y = Lat, color = silt_fraction), 
             size = 1.5) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("SSP3-RCP7.0") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(color = "Silt Fraction") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  #scale_color_gradient(low = "slateblue4", high = "firebrick1")  
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 0.5, limits = c(0, 1))  


g_ssp585_silt = ggmap(worldmap) +
  geom_point(data = d_soil_ssp585, 
             aes(x = Lon, y = Lat, color = silt_fraction), 
             size = 1.5) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("SSP5-RCP8.5") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(color = "Silt Fraction") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  #scale_color_gradient(low = "slateblue4", medium = "white", high = "firebrick1")  
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 0.5, limits = c(0, 1))  

plot_grid(g_picontrol_silt, g_ssp126_silt, g_ssp370_silt, g_ssp585_silt, ncol = 1)

ggsave("Scripts/Plots/Descriptive/Soil/pdf/silt.pdf", width = 20, height = 15)
ggsave("Scripts/Plots/Descriptive/Soil/png/silt.png", width = 20, height = 15)



