################################################################################
############################ Master's Thesis ###################################
################################################################################

############################ Model: Data Prep ##################################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/02_FPCA/functions.R")
source("Scripts/MA_FDA_veg/01_Description/utils.R")

## Load libraries
library(dplyr)
library(data.table)
library(stringr)
library(cowplot)
library(ggplot2)
library(MFPCA)
library(foreach)
library(funData)
library(abind)
library(tidyverse)
library(RColorBrewer)

scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")
pfts = c("BNE", "IBS", "otherC", "TeBS","Tundra")

# The data provided can be classified into these four groups:
# - PFT-dependent and time-varying: Nuptake
# - PFT-dependent only: initial_recruitment, recruitment_ten_years, previous_state, time_since_dist
# - time-varying only: tas_yearlymax, tas_yearlymin, tas_yearlsmeam, pr_yearlysum, Nuptake_total
# - Location dependent only: sand_fraction, silt_fraction, clay_fraction, bulkdensity_soil, ph_soil, soilcarbon


## Load data
pid = 1

# Get target variable: MFPCA scores
plot_data_all = read.table("Scripts/Plots/MFPCA/plot_data/plot_data_all.txt")

# Get locations
locs_disturbed = read.table("Scripts/Plots/MFPCA/plot_data/locs_disturbed.txt")

# Get other time-constant soil covariates
for (scen in scenarios){
  d_scen_2015_2040 = fread(paste0("Data/data_", scen, "_2015_2040.csv"))
  
  d_soil_scen = d_scen_2015_2040 %>%
    filter(PID == pid) %>%
    distinct(Lon, Lat, sand_fraction, silt_fraction, clay_fraction, bulkdensity_soil, ph_soil, soilcarbon, time_since_dist) %>%
    mutate(Scenario = long_names_scenarios(scen))
  
  d_eco_scen = d_scen_2015_2040 %>%
    filter(PID == pid) %>%
    distinct(Lon,Lat,PFT,Year,initial_recruitment, recruitment_ten_years, previous_state, age, Nuptake) %>%
    mutate(Scenario = long_names_scenarios(scen))
  
  assign(paste0("d_soil_", scen), d_soil_scen)
  assign(paste0("d_eco_", scen), d_eco_scen)
  
}

for (var in c("soil", "eco")){
  
  d_var = rbind(get(paste0("d_", var, "_picontrol")),
                get(paste0("d_", var, "_ssp126")), 
                get(paste0("d_", var, "_ssp370")), 
                get(paste0("d_", var, "_ssp585")))
  
  if (var == "soil"){
    # Reorder d_var to match target
    d_var$key = paste(d_var$Lon, d_var$Lat, d_var$Scenario)
    locs_disturbed$key = paste(locs_disturbed$Lon, locs_disturbed$Lat, locs_disturbed$Scenario)
    
    index = match(locs_disturbed$key, d_var$key)
    d_var = d_var[index,]
    d_var$key = NULL
  }
  
  assign(paste0("d_", var), d_var)
}


######################### Further data pre-processing ##########################

### Ecological PFT-dependent variables

# Cut data to 100 years of resoilery
d_eco = d_eco %>%
  filter(age <= 100)

# Get time-independent variables
d_eco_ind = d_eco %>%
  distinct(Lon,Lat,PFT,initial_recruitment, recruitment_ten_years, previous_state, Scenario) %>%
  group_by(PFT) %>%
  ungroup()

# Spread the data so that each PFT category’s values are in separate columns
d_eco_wide <- d_eco_ind %>%
  pivot_wider(names_from = PFT, values_from = c(initial_recruitment, recruitment_ten_years, previous_state))

# Sort the resulting data according to right order of locations

d_eco_wide$key = paste(d_eco_wide$Lon, d_eco_wide$Lat, d_eco_wide$Scenario)

index = match(locs_disturbed$key, d_eco_wide$key)
d_eco_wide = d_eco_wide[index,]
d_eco_wide$key = NULL

### Get full data set
write.table(d_soil, "Scripts/MA_FDA_veg/04_Model/Data/d_soil.txt")
write.table(d_eco_wide, "Scripts/MA_FDA_veg/04_Model/Data/d_eco_wide.txt")
