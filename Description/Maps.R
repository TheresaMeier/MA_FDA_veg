setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/Description/utils.R")

# Load libraries
library(duckdb)
library(tidyverse)
library(ggplot2)
library(ggplot2)
library(ggmap)
library(dplyr)
library(gganimate)
library(gifski)
library(gridExtra)

# Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)

# Replace 'YOUR_API_KEY' with your actual API key
register_google(key = "AIzaSyATIWAZ4gtgJhdH6GP_E8iSubdFh6XQ32Y")

scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")

# Get map
worldmap <- get_map(location = c(lon =-4.068561, lat = 58.87355), zoom = 1)

for (iScen in scenarios){
  for (iYear in c(2015,2100,2200,2296)){
    data_loc <- dbGetQuery(con, paste0("SELECT Lon,Lat,PFT,cmass FROM '", iScen, "_d150_cmass' WHERE PID == 0 AND Year == ", iYear))
    
    # Get most apparent PFT
    max_data_loc <- data_loc %>%
      group_by(Lon, Lat) %>%
      slice(which.max(cmass)) %>%
      ungroup() %>%
      mutate(PFT = long_names_pfts(tolower(PFT)))
    
    
    # Plot the map with points
    
    g = ggmap(worldmap) +
      geom_point(data = max_data_loc, aes(x = Lon, y = Lat, color = PFT), size = 1) +
      ggtitle(paste(long_names_scenarios(iScen), " - ", iYear)) +
      xlab("Longitude") +
      ylab("Latitude") +
      scale_color_manual(name = "Dominant vegetation", drop = TRUE,
                         values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                                    "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
      theme_bw() +
      theme(text = element_text(size = 5),plot.title = element_text(size = 10))
    
    assign(paste0("g_",iScen,"_",iYear), g)
    
    #ggsave(paste0("Scripts/Plots/Descriptive/g_",iScen,"_",iYear,".png"),g)
  }
}


gifski(c("Scripts/Plots/Descriptive/g_picontrol_2015.png","Scripts/Plots/Descriptive/g_picontrol_2100.png","Scripts/Plots/Descriptive/g_picontrol_2200.png","Scripts/Plots/Descriptive/g_picontrol_2296.png"), gif_file = "Scripts/Plots/Descriptive/gif_picontrol.gif", delay = 2)
gifski(c("Scripts/Plots/Descriptive/g_ssp126_2015.png","Scripts/Plots/Descriptive/g_ssp126_2100.png","Scripts/Plots/Descriptive/g_ssp126_2200.png","Scripts/Plots/Descriptive/g_ssp126_2296.png"), gif_file = "Scripts/Plots/Descriptive/gif_ssp126.gif", delay = 2)
gifski(c("Scripts/Plots/Descriptive/g_ssp370_2015.png","Scripts/Plots/Descriptive/g_ssp370_2100.png","Scripts/Plots/Descriptive/g_ssp370_2200.png","Scripts/Plots/Descriptive/g_ssp370_2296.png"), gif_file = "Scripts/Plots/Descriptive/gif_ssp370.gif", delay = 2)
gifski(c("Scripts/Plots/Descriptive/g_ssp585_2015.png","Scripts/Plots/Descriptive/g_ssp585_2100.png","Scripts/Plots/Descriptive/g_ssp585_2200.png","Scripts/Plots/Descriptive/g_ssp585_2296.png"), gif_file = "Scripts/Plots/Descriptive/gif_ssp585.gif", delay = 2)

g_2015 = grid.arrange(g_picontrol_2015,g_ssp126_2015, g_ssp370_2015, g_ssp585_2015, ncol = 2)
ggsave("Scripts/Plots/Descriptive/g_veg_2015.png",g_2015)

g_2100 = grid.arrange(g_picontrol_2100,g_ssp126_2100, g_ssp370_2100, g_ssp585_2100, ncol = 2)
ggsave("Scripts/Plots/Descriptive/g_veg_2100.png",g_2100)

g_2200 = grid.arrange(g_picontrol_2200,g_ssp126_2200, g_ssp370_2200, g_ssp585_2200, ncol = 2)
ggsave("Scripts/Plots/Descriptive/g_veg_2200.png",g_2200)

g_2296 = grid.arrange(g_picontrol_2296,g_ssp126_2296, g_ssp370_2296, g_ssp585_2296, ncol = 2)
ggsave("Scripts/Plots/Descriptive/g_veg_2296.png",g_2296)

gifski(c("Scripts/Plots/Descriptive/g_veg_2015.png","Scripts/Plots/Descriptive/g_veg_2100.png","Scripts/Plots/Descriptive/g_veg_2200.png","Scripts/Plots/Descriptive/g_veg_2296.png"), gif_file = "Scripts/Plots/Descriptive/gif_veg.gif", delay = 3)


########################### World Map of Grid Cells ############################

ggmap(worldmap) +
  geom_point(data = data_loc, aes(x = Lon, y = Lat), size = 0.5) +
  geom_hline(yintercept = min(data_loc$Lat), lty = 1, color = "red", lwd = 1) +
  geom_hline(yintercept = max(data_loc$Lat), lty = 1, color = "red", lwd = 1) +
  geom_ribbon(data = NULL, aes(ymin = -Inf, ymax = min(data_loc$Lat)), fill = "grey", alpha = 0.5) +
  geom_ribbon(data = NULL, aes(ymin = max(data_loc$Lat), ymax = Inf), fill = "grey", alpha = 0.5) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() +
  theme(text = element_text(size = 5))

ggsave("Scripts/Plots/Descriptive/worldmap.png")
ggsave("Scripts/Plots/Descriptive/worldmap.pdf")
