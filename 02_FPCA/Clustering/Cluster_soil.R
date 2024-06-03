################################################################################
############################ Master's Thesis ###################################
################################################################################

########################### CLimate and Soil data ##############################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/Description/utils.R")

## Load libraries
library(dplyr)
library(data.table)
library(stringr)

scenarios = c("picontrol", "ssp126", "ssp585")

## Load data

# Get cluster assignments
plot_data_all = read.table("Scripts/MA_FDA_veg/FPCA/plot_data/plot_data_all.txt")

for (scen in scenarios){
  d_scen_2015_2040 = fread(paste0("Data/data_", scen, "_2015_2040.csv"))
  
  d_soil_scen = d_scen_2015_2040 %>%
    distinct(Lon,Lat,sand_fraction,silt_fraction, clay_fraction, bulkdensity_soil, ph_soil, soilcarbon)
  
  d_scen = as.data.frame(cbind(rownames(plot_data_all[plot_data_all$scenario == long_names_scenarios(scen),]), plot_data_all[plot_data_all$scenario == long_names_scenarios(scen),"Cluster"]))
  colnames(d_scen) = c("Scenario", "Cluster")
  
  # Transform Location to Lon and Lat values
  numeric_values <- as.numeric(unlist(str_extract_all(d_scen$Scenario, "\\d+")))
  
  # Remove the ones from the vector
  numeric_values <- numeric_values[!(numeric_values < 8)]
  
  d_scen$Lon = numeric_values[seq(1, length(numeric_values), by = 4)] + 0.01 * numeric_values[seq(2, length(numeric_values), by = 4)]
  d_scen$Lat = numeric_values[seq(3, length(numeric_values), by = 4)] + 0.01 * numeric_values[seq(4, length(numeric_values), by = 4)]
  
  d_scen[which(!paste(d_scen$Lon, d_scen$Lat) %in% paste(d_soil_scen$Lon, d_soil_scen$Lat)),]$Lon = -d_scen[which(!paste(d_scen$Lon, d_scen$Lat) %in% paste(d_soil_scen$Lon, d_soil_scen$Lat)),]$Lon
  print(nrow(anti_join(d_scen, d_soil_scen, by = c("Lon", "Lat"))))
  
  # Merge data sets
  d_scen_all = inner_join(d_scen, d_soil_scen, by = c("Lon", "Lat")) %>%
    mutate(Scenario = long_names_scenarios(scen))
  
  assign(paste0("d_", scen, "_all"), d_scen_all)
}

# Merge data sets for all scenarios
d_all = rbind(d_picontrol_all, d_ssp126_all, d_ssp585_all)

# Get soil factors means
d_all_means = d_all %>%
  group_by(Cluster) %>%
  summarize(mean_sand_fraction = mean(sand_fraction),
            mean_silt_fraction = mean(silt_fraction),
            mean_clay_fraction = mean(clay_fraction),
            mean_bulkdensity_soil = mean(bulkdensity_soil),
            mean_ph_soil = mean(ph_soil),
            mean_soilcarbon = mean(soilcarbon),
            min_sand_fraction = min(sand_fraction),
            min_silt_fraction = min(silt_fraction),
            min_clay_fraction = min(clay_fraction),
            min_bulkdensity_soil = min(bulkdensity_soil),
            min_ph_soil = min(ph_soil),
            min_soilcarbon = min(soilcarbon),
            max_sand_fraction = max(sand_fraction),
            max_silt_fraction = max(silt_fraction),
            max_clay_fraction = max(clay_fraction),
            max_bulkdensity_soil = max(bulkdensity_soil),
            max_ph_soil = max(ph_soil),
            max_soilcarbon = max(soilcarbon))


