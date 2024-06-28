################################################################################
############################ Master's Thesis ###################################
################################################################################

########################### Model: First Trials ################################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/02_FPCA/functions.R")
source("Scripts/MA_FDA_veg/03_MFPCA/MFPCA/MFPCA_calculation.R")
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

scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")

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

# Get other time-constant covariates
for (scen in scenarios){
  d_scen_2015_2040 = fread(paste0("Data/data_", scen, "_2015_2040.csv"))
  
  d_cov_scen = d_scen_2015_2040 %>%
    filter(PID == pid) %>%
    distinct(Lon, Lat, sand_fraction, silt_fraction, clay_fraction, bulkdensity_soil, ph_soil, soilcarbon, time_since_dist) %>%
    mutate(Scenario = long_names_scenarios(scen))
  
  d_PFT_scen = d_scen_2015_2040 %>%
    filter(PID == pid) %>%
    distinct(Lon,Lat,PFT,Year,initial_recruitment, recruitment_ten_years, previous_state, age, Nuptake) %>%
    mutate(Scenario = long_names_scenarios(scen))
  
  d_climate_scen = d_scen_2015_2040 %>%
    filter(PID == pid) %>%
    distinct(Lon,Lat,Year, age, Nuptake_total,
             tas_yearlymax, tas_yearlymin, tas_yearlymeam, pr_yearlysum) %>%
    mutate(Scenario = long_names_scenarios(scen))
  
  assign(paste0("d_cov_", scen), d_cov_scen)
  assign(paste0("d_PFT_", scen), d_PFT_scen)
  assign(paste0("d_climate_", scen), d_climate_scen)
  
}

for (var in c("cov", "PFT", "climate")){
  
  d_var = rbind(get(paste0("d_", var, "_picontrol")),
                get(paste0("d_", var, "_ssp126")), 
                get(paste0("d_", var, "_ssp370")), 
                get(paste0("d_", var, "_ssp585")))
  
  if (var == "cov"){
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

### Climate data: Run MFPCA on mean temperature and precipitation
## Create funData object (1806x100x2)

# Cut data to 100 years of recovery
d_climate = d_climate %>%
  filter(age <= 100)

# store Lon/Lat values
locs_climate = d_climate %>% distinct(Lon,Lat,Scenario)

# Get wider format

d_climate_temp_wide <- d_climate %>%
  distinct(Lon,Lat, age, tas_yearlymeam, Scenario) %>%
  pivot_wider(names_from = c(Lon,Lat,Scenario), values_from = tas_yearlymeam, values_fn = list) %>%
  arrange(age) %>%
  rename_with(~ gsub("_", "/", .), everything()) 


d_climate_precip_wide <- d_climate %>%
  distinct(Lon,Lat, age, pr_yearlysum, Scenario) %>%
  pivot_wider(names_from = c(Lon,Lat,Scenario), values_from = pr_yearlysum, values_fn = list) %>%
  arrange(age) %>%
  rename_with(~ gsub("_", "/", .), everything()) 

# Join data sets 
d_climate_precip_wide = t(d_climate_precip_wide[,-1])
d_climate_temp_wide = t(d_climate_temp_wide[,-1])
d_climate_wide = abind(list(d_climate_temp_wide, d_climate_precip_wide), along = 3)

funData_tmp = funData(argvals = list(0:100, 1:2), X = d_climate_wide)
funData_all = multiFunData(funData_tmp)

## Run MFPCA
MFPCA_climate = MFPCA2(funData_all, M=10, uniExpansions = list(list(type = "given", funData_tmp, ortho = TRUE)), fit = FALSE, approx.eigen = FALSE)

# Variance which is accounted for:
round(MFPCA_climate$values[1:10]/sum(MFPCA_climate$values) *100,3)

# --> only first PC is necessary
# store scores

scores_climate = as.data.frame(MFPCA_climate$scores[,1])

# Reorder scores according to right order of locations
locs_climate$key = paste(locs_climate$Lon, locs_climate$Lat, locs_climate$Scenario)

index = match(locs_disturbed$key, locs_climate$key)
scores_climate = scores_climate[index,]
locs_climate$key = NULL

### Ecological PFT-dependent variables

# Cut data to 100 years of recovery
d_PFT = d_PFT %>%
  filter(age <= 100)

# Get time-independent variables
d_PFT_ind = d_PFT %>%
  distinct(Lon,Lat,PFT,initial_recruitment, recruitment_ten_years, previous_state, Scenario) %>%
  group_by(PFT) %>%
  ungroup()

# Spread the data so that each PFT category’s values are in separate columns
d_PFT_wide <- d_PFT_ind %>%
  pivot_wider(names_from = PFT, values_from = c(initial_recruitment, recruitment_ten_years, previous_state))

# Sort the resulting data according to right order of locations

d_PFT_wide$key = paste(d_PFT_wide$Lon, d_PFT_wide$Lat, d_PFT_wide$Scenario)

index = match(locs_disturbed$key, d_PFT_wide$key)
d_PFT_wide = d_PFT_wide[index,]
d_PFT_wide$key = NULL

### Get full data set
combined_d = cbind(plot_data_all[,c(1,2)], d_cov, d_PFT_wide[,c(4:23)], scores_climate)
combined_d = combined_d[,-c(28:31)]
################################### Model ######################################

# Standardize PC1 and PC2
combined_d <- combined_d %>%
  mutate(
    PC1_standardized = (PC1 - mean(PC1)) / sd(PC1),
    PC2_standardized = (PC2 - mean(PC2)) / sd(PC2)
  )

lm_1 = lm(cbind(PC1_standardized, PC2_standardized) ~., data = combined_d %>%
            select(-PC1, -PC2, -clay_fraction))
summary(lm_1)

predictions = predict(lm_1)
predictions[c(1:10),]
combined_d[c(1:10),] %>% select(PC1_standardized, PC2_standardized)

### Model evaluation
residuals_pc1 <- residuals(lm_1)[, 1]
residuals_pc2 <- residuals(lm_1)[, 2]

# Create residual plots for PC1
par(mfrow = c(2, 2))  # Set up a 2x2 plotting area

plot(fitted(lm_1)[, 1], residuals_pc1, main = "Residuals vs Fitted (PC1)",
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")

qqnorm(residuals_pc1, main = "Normal Q-Q (PC1)")
qqline(residuals_pc1, col = "red")

plot(fitted(lm_1)[, 2], residuals_pc2, main = "Residuals vs Fitted (PC2)",
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")

qqnorm(residuals_pc2, main = "Normal Q-Q (PC2)")
qqline(residuals_pc2, col = "red")
