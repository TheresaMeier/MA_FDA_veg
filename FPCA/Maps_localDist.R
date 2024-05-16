################################################################################
############################ Master's Thesis ###################################
################################################################################

#################### Exploratory Analysis: (M)FPCA #############################

# Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/FPCA/functions.R")
source("Scripts/Description/utils.R")

# Load libraries
library(duckdb)
library(tidyverse)
library(ggplot2)
library(fda)
library(ggmap)
library(gifski)
library(RColorBrewer)
library(gridExtra)

# Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)

################################## Get data ####################################

d_picontrol_2015_2040_cmass_1 <- get_data_scenario('picontrol', 2015,2040,'cmass',1)
d_picontrol_2015_2040_cmass_1_tundra <- d_picontrol_2015_2040_cmass_1[d_picontrol_2015_2040_cmass_1$PFT == "Tundra",]

d_picontrol_2015_2040_cmass_1_tundra <- d_picontrol_2015_2040_cmass_1_tundra %>%
  mutate(rel_index = ifelse(relative >0.5, 1,0))

# Plot the trajectories
ggplot() + 
  geom_line(data = d_picontrol_2015_2040_cmass_1_tundra, linewidth = .05, alpha = .25,
            aes(x = age, y = relative, col = PFT,group = interaction(Lon, Lat, PID, PFT))) +
  #geom_line(data = df_mean, aes(x = age, y = relative, color = PFT, group = PFT), linewidth = 2) +
  scale_x_continuous(name = "Year after disturbance", expand = c(0,0), limits = c(0, 100)) +
  scale_y_continuous(name = "Share of aboveground biomass", expand = c(0,0),
                     breaks = c(0.25, 0.50, 0.75, 1.00)) +
  scale_color_manual(name = "Dominant vegetation", drop = TRUE,
                     values = c("Tundra" = "#009E73")) +
  ggtitle("Recovery trajectories with disturbances between 2015 and 2040") +
  add_common_layout() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10))


d_picontrol_2015_2040_cmass_1_tundra_start <- d_picontrol_2015_2040_cmass_1_tundra %>%
  group_by(Lat, Lon) %>%
  filter(age == 10) %>%
  ungroup()

########################## Plot Spatial Dependencies ###########################

register_google(key = "AIzaSyATIWAZ4gtgJhdH6GP_E8iSubdFh6XQ32Y")

# Get map
worldmap <- get_map(location = c(lon =-4.068561, lat = 58.87355), zoom = 1)

scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")

pfts = c("BNE", "IBS", "otherC", "TeBS", "Tundra")

for (iScen in scenarios){
  
  d_2015_2040_cmass_1 <- get_data_scenario(iScen, 2015,2040,'cmass',1)
  
  for (iPFT in pfts){
    
    d_2015_2040_cmass_1_pft <- d_2015_2040_cmass_1 %>% 
      filter(PFT == iPFT)
    
    for (iYear in seq(0,100,by=1)){
      d_2015_2040_cmass_1_pft_iYear <- d_2015_2040_cmass_1_pft %>%
        group_by(Lat, Lon) %>%
        filter(age == iYear) %>%
        ungroup()
      
      g = ggmap(worldmap) +
        geom_point(data = d_2015_2040_cmass_1_pft_iYear, aes(x = Lon, y = Lat, col = relative), size = 1) +
        xlab("Longitude") +
        ylab("Latitude") +
        ylim(c(0,80)) +
        xlim(c(-180,180)) +
        theme_bw() +
        theme(text = element_text(size = 10), plot.title = element_text(size = 12)) +
        #ggtitle(paste(long_names_scenarios(iScen),"-",long_names_pfts(tolower(iPFT)),"- Year", iYear, "after disturbance")) +
        ggtitle(long_names_scenarios(iScen)) +
        scale_color_gradient(name = "Share", limits = c(0, 1), low = brewer.pal(9, "YlOrRd")[2], high = brewer.pal(9, "YlOrRd")[9])
      
      assign(paste0("g_",iScen,"_", iPFT,"_",iYear), g)
      dir_path <- paste0("Scripts/Plots/FPCA/Maps/",iScen,"_",iPFT)
      
      # Create the directory if it doesn't exist
      if (!file.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
      }
      ggsave(paste0(dir_path,"/g_",iScen,"_", iPFT,"_",iYear,".png"),g)
    }
    
    png_files <- list.files(dir_path, pattern = "\\.png$", full.names = TRUE)
    file_numbers <- as.numeric(gsub(".*\\D(\\d+)\\.png", "\\1", png_files))
    ordered_png_files <- png_files[order(file_numbers)]
    
    gifski(ordered_png_files,width = 1000, height = 600, gif_file = paste0("Scripts/Plots/FPCA/Maps/gif_",iScen,"_", iPFT,".gif"), delay = 0.1, loop = FALSE)
  }
}

# 2x2 plots
for (iPFT in pfts){
  
  dir_path <- paste0("Scripts/Plots/FPCA/Maps/",iPFT)
  # Create the directory if it doesn't exist
  if (!file.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  for (iYear in seq(0,100,by=1)){
    g = grid.arrange(get(paste0("g_picontrol_", iPFT,"_", iYear)), get(paste0("g_ssp126_", iPFT,"_", iYear)), get(paste0("g_ssp370_", iPFT,"_", iYear)), get(paste0("g_ssp585_", iPFT,"_", iYear)),
                     ncol = 2,
                     top = textGrob(paste(iPFT, "- Year", iYear, "after disturbance"),gp=gpar(fontsize=15,font=2)))
    assign(paste0("g_", iPFT, "_", iYear),g)
    
    ggsave(paste0(dir_path,"/g_", iPFT,"_",iYear,".png"),g)
  }
  
  png_files <- list.files(dir_path, pattern = "\\.png$", full.names = TRUE)
  file_numbers <- as.numeric(gsub(".*\\D(\\d+)\\.png", "\\1", png_files))
  ordered_png_files <- png_files[order(file_numbers)]
  
  gifski(ordered_png_files, width = 1000, height = 600, gif_file = paste0("Scripts/Plots/FPCA/Maps/gif_", iPFT,".gif"), delay = 0.1, loop = FALSE)
}
