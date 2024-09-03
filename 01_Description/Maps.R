################################################################################
############################ Master's Thesis ###################################
################################################################################

############################# Plot world maps ##################################

# Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/01_Description/utils.R")
source("Scripts/MA_FDA_veg/02_FPCA/functions.R")

# Load libraries
library(duckdb)
library(tidyverse)
library(ggplot2)
library(ggplot2)
library(ggmap)
library(dplyr)
library(gridExtra)
library(cowplot)

# Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)

register_google(key = "AIzaSyA_eUpOhj7hoPyzyynWvyMqcGEA1Z_SZVY")
scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")

# Get map
worldmap <- get_map(location = c(lon =-4.068561, lat = 58.87355), zoom = 1)

########################### World Map of Grid Cells ############################

data_loc <- dbGetQuery(con, paste0("SELECT Lon,Lat,PFT,cmass FROM 'picontrol_d150_cmass' WHERE PID == 0 AND Year == 2010")) %>% 
  distinct(Lon,Lat)

ggmap(worldmap) +
  geom_point(data = data_loc, aes(x = Lon, y = Lat), size = 0.5) + ylim(0,80) +
  geom_hline(yintercept = min(data_loc$Lat)-1, lty = 1, color = "red", lwd = 1) +
  geom_hline(yintercept = max(data_loc$Lat)+1, lty = 1, color = "red", lwd = 1) +
  geom_ribbon(data = NULL, aes(ymin = -Inf, ymax = (min(data_loc$Lat)-1)), fill = "grey", alpha = 0.5) +
  geom_ribbon(data = NULL, aes(ymin = (max(data_loc$Lat)+1), ymax = Inf), fill = "grey", alpha = 0.5) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() + 
  theme(text = element_text(size = 10))

ggsave("Scripts/Plots/Descriptive/worldmap.png")
ggsave("Scripts/Plots/Descriptive/worldmap.pdf")


################### Plot disturbed grid cells of patch 1 #######################

start_year = 2015
end_year = 2040
pft = "BNE"      # Choose from Tundra, BNE, IBS, otherC, TeBS
pid = 1           # Choose pid (int) or 'all'

## Get data for all four scenarios in the appropriate shape

d_picontrol = get_data_fpca("picontrol", start_year, end_year,pid,pft)
d_ssp126 = get_data_fpca("ssp126", start_year, end_year,pid,pft)
d_ssp370 = get_data_fpca("ssp370", start_year, end_year,pid,pft)
d_ssp585 = get_data_fpca("ssp585", start_year, end_year,pid,pft)

d_picontrol_loc = d_picontrol[[2]]
d_picontrol = d_picontrol[[1]]
d_ssp126_loc = d_ssp126[[2]]
d_ssp126 = d_ssp126[[1]]
d_ssp370_loc = d_ssp370[[2]]
d_ssp370 = d_ssp370[[1]]
d_ssp585_loc = d_ssp585[[2]]
d_ssp585 = d_ssp585[[1]]

g_picontrol = ggmap(worldmap) +
  geom_point(data = d_picontrol_loc, aes(x = Lon, y = Lat), size = 0.5) +
  xlab("Longitude") + ylim(35,75) + xlim(-200,200) +
  ylab("Latitude") +
  theme_bw() + ggtitle("Control") +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))

g_ssp126 = ggmap(worldmap) +
  geom_point(data = d_ssp126_loc, aes(x = Lon, y = Lat), size = 0.5) +
  xlab("Longitude") + ylim(35,75) + xlim(-200,200) +
  ylab("Latitude") +
  theme_bw() + ggtitle("SSP1-RCP2.6") +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))

g_ssp370 = ggmap(worldmap) +
  geom_point(data = d_ssp370_loc, aes(x = Lon, y = Lat), size = 0.5) +
  xlab("Longitude") + ylim(35,75) + xlim(-200,200) +
  ylab("Latitude") +
  theme_bw() + ggtitle("SSP3-RCP7.0") +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))

g_ssp585 = ggmap(worldmap) +
  geom_point(data = d_ssp585_loc, aes(x = Lon, y = Lat), size = 0.5) +
  xlab("Longitude") + ylim(35,75) + xlim(-200,200) +
  ylab("Latitude") +
  theme_bw() + ggtitle("SSP5-RCP8.5") +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))

plot_grid(g_picontrol, g_ssp126, nrow = 1)

ggsave("Scripts/Plots/FPCA/PCs_2015_2040/Plots_MA/map_disturbed_patches_1.png", width = 10, height = 20)
ggsave("Scripts/Plots/FPCA/PCs_2015_2040/Plots_MA/map_disturbed_patches_1.pdf", width = 10, height = 20)

plot_grid(g_ssp370, g_ssp585, nrow = 1)

ggsave("Scripts/Plots/FPCA/PCs_2015_2040/Plots_MA/map_disturbed_patches_2.png", width = 10, height = 20)
ggsave("Scripts/Plots/FPCA/PCs_2015_2040/Plots_MA/map_disturbed_patches_2.pdf", width = 10, height = 20)
