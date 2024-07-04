################################################################################
############################ Master's Thesis ###################################
################################################################################

###################### Model: (M)FPCA for climate data #########################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/02_FPCA/functions.R")
source("Scripts/MA_FDA_veg/01_Description/utils.R")
source("Scripts/MA_FDA_veg/03_MFPCA/MFPCA/MFPCA_calculation.R")

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

# The data provided can be classified into these four groups:
# - PFT-dependent and time-varying: Nuptake
# - PFT-dependent only: initial_recruitment, recruitment_ten_years, previous_state, time_since_dist
# - time-varying only: tas_yearlymax, tas_yearlymin, tas_yearlsmeam, pr_yearlysum, Nuptake_total
# - Location dependent only: sand_fraction, silt_fraction, clay_fraction, bulkdensity_soil, ph_soil, soilcarbon

#######################
## Get 450 different colors for plotting

palettes <- c("Set1", "Set2", "Set3", "Dark2", "Paired", "Accent")
all_colors <- character()
for (palette_name in palettes) {
  palette <- brewer.pal(20, palette_name)
  all_colors <- c(all_colors, palette)
}
palette_450 <- unique(all_colors[1:470])
pal <- colorRampPalette(palette_450)

################################ Load data #####################################
pid = 1

# Get target variable: MFPCA scores
plot_data_all = read.table("Scripts/Plots/MFPCA/plot_data/plot_data_all.txt")

# Get locations
locs_disturbed = read.table("Scripts/Plots/MFPCA/plot_data/locs_disturbed.txt")

# Get climate covariates
for (scen in scenarios){
  d_scen_2015_2040 = fread(paste0("Data/data_", scen, "_2015_2040.csv"))
  
  d_climate_scen = d_scen_2015_2040 %>%
    filter(PID == pid) %>%
    distinct(Lon,Lat,Year, age, Nuptake_total,
             tas_yearlymax, tas_yearlymin, tas_yearlymeam, pr_yearlysum) %>%
    mutate(Scenario = long_names_scenarios(scen))
  
  assign(paste0("d_climate_", scen), d_climate_scen)
  
}

d_climate = rbind(d_climate_picontrol,d_climate_ssp126, d_climate_ssp370, d_climate_ssp585)

############################### Perform MFPCA ##################################
## Get functional fit
# Create funData object(s)

# Cut data to 100 years of recovery
d_climate = d_climate %>%
  filter(age <= 100 & age > 0)

d_climate = d_climate[!duplicated(d_climate),]

# store Lon/Lat values
locs_climate = d_climate %>% distinct(Lon,Lat,Scenario)

# Get wider format

d_climate_temp_wide <- d_climate %>%
  distinct(Lon,Lat, age, tas_yearlymeam, Scenario) %>%
  pivot_wider(names_from = c(Lon,Lat,Scenario), values_from = tas_yearlymeam) %>%
  arrange(age) %>%
  rename_with(~ gsub("_", "/", .), everything()) 


d_climate_precip_wide <- d_climate %>%
  distinct(Lon,Lat, age, pr_yearlysum, Scenario) %>%
  pivot_wider(names_from = c(Lon,Lat,Scenario), values_from = pr_yearlysum) %>%
  arrange(age) %>%
  rename_with(~ gsub("_", "/", .), everything()) 

# Get funData objects
funData_temp = funData(argvals = 1:100, X = t(d_climate_temp_wide[,-1]))
funData_precip = funData(argvals = 1:100, X = t(d_climate_precip_wide[,-1]))
funData_all = multiFunData(funData_temp, funData_precip)

saveRDS(funData_temp, "Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/funData_temp.rds")
saveRDS(funData_precip, "Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/funData_precip.rds")
saveRDS(funData_all, "Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/funData_all.rds")

## Run MFPCA
uniExpansions <- list(list(type = "uFPCA", npc = 10),
                      list(type = "uFPCA", npc = 10))

MFPCA_climate = MFPCA_2(funData_all, M = 8, uniExpansions = uniExpansions, fit = TRUE, approx.eigen = FALSE)

MFPCA_temp = MFPCA_2(multiFunData(funData_temp), M = 8, uniExpansions = list(list(type = "uFPCA", npc = 10)), fit = TRUE)
MFPCA_precip = MFPCA_2(multiFunData(funData_precip), M = 8, uniExpansions = list(list(type = "uFPCA", npc = 10)), fit = TRUE)

saveRDS(MFPCA_temp, "Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/MFPCA_temp.rds")
saveRDS(MFPCA_precip, "Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/MFPCA_precip.rds")
saveRDS(MFPCA_climate, "Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/MFPCA_climate.rds")

# Variance which is accounted for:
round(MFPCA_climate$values[1:10]/sum(MFPCA_climate$values) *100,3)
round(MFPCA_temp$values[1:10]/sum(MFPCA_temp$values) *100,3)
round(MFPCA_precip$values[1:10]/sum(MFPCA_precip$values) *100,3)

## Plot the functional fits including mean functions
bounds = t(matrix(c(1,434,435,876,877,1338,1339,1803), nrow = 2))

# Precipitation
pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/BasisRep/BasisRep_precip.pdf"), width = 8, height = 6.5)
plot(funData_precip, col = pal(1803)[1:1803],
     xlim = c(1,100), type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = bquote("Yearly precipitation in kg/m"^2), 
     main = paste("Functional fit - Yearly precipitation") )
lines(x=1:100, y = MFPCA_precip$meanFunction@.Data[[1]]@X, col = "black", lwd = 3)
legend("topright", legend = "Mean Function", col = "black", lwd = 3)
dev.off()

# Per scenario
pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/BasisRep/BasisRep_precip_scen.pdf"), width = 10, height = 11)
par(mfrow =(c(2,2)))
plot(funData_precip[1:434], col = pal(434)[1:434],
     xlim = c(1,100), ylim = c(0,4000), type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = bquote("Yearly precipitation in kg/m"^2), 
     main = paste("Control") )

plot(funData_precip[bounds[2,1]:bounds[2,2]], col = pal(500)[1:500],
     xlim = c(1,100), ylim = c(0,4000),type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = bquote("Yearly precipitation in kg/m"^2), 
     main = paste("SSP1-RCP2.6") )

plot(funData_precip[bounds[3,1]:bounds[3,2]], col = pal(500)[1:500],
     xlim = c(1,100), ylim = c(0,4000),type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = bquote("Yearly precipitation in kg/m"^2), 
     main = paste("SSP3-RCP7.0") )

plot(funData_precip[bounds[4,1]:bounds[4,2]], col = pal(500)[1:500],
     xlim = c(1,100), ylim = c(0,4000),type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = bquote("Yearly precipitation in kg/m"^2), 
     main = paste("SSP5-RCP8.5") )
mtext("Functional fit - Yearly precipitation", side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.5)
dev.off()

# Save as png as well
png(paste0("Scripts/Plots/Model/(M)FPCA/png/BasisRep/BasisRep_precip.png"), width = 800, height = 650)
plot(funData_precip, col = pal(1803)[1:1803],
     xlim = c(1,100), type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = bquote("Yearly precipitation in kg/m"^2), 
     main = paste("Functional fit - Yearly precipitation") )
lines(x=1:100, y = MFPCA_precip$meanFunction@.Data[[1]]@X, col = "black", lwd = 3)
legend("topright", legend = "Mean Function", col = "black", lwd = 3)
dev.off()


png(paste0("Scripts/Plots/Model/(M)FPCA/png/BasisRep/BasisRep_precip_scen.png"), width = 1000, height = 1100)
par(mfrow =(c(2,2)))
plot(funData_precip[1:434], col = pal(434)[1:434],
     xlim = c(1,100), ylim = c(0,4000), type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = bquote("Yearly precipitation in kg/m"^2), 
     main = paste("Control") )

plot(funData_precip[bounds[2,1]:bounds[2,2]], col = pal(500)[1:500],
     xlim = c(1,100), ylim = c(0,4000),type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = bquote("Yearly precipitation in kg/m"^2), 
     main = paste("SSP1-RCP2.6") )

plot(funData_precip[bounds[3,1]:bounds[3,2]], col = pal(500)[1:500],
     xlim = c(1,100), ylim = c(0,4000),type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = bquote("Yearly precipitation in kg/m"^2), 
     main = paste("SSP3-RCP7.0") )

plot(funData_precip[bounds[4,1]:bounds[4,2]], col = pal(500)[1:500],
     xlim = c(1,100), ylim = c(0,4000),type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = bquote("Yearly precipitation in kg/m"^2), 
     main = paste("SSP5-RCP8.5") )
mtext("Functional fit - Yearly precipitation", side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.5)
dev.off()

# Mean temperature
pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/BasisRep/BasisRep_mean_temp.pdf"), width = 8, height = 6.5)
plot(funData_temp, col = pal(1803)[1:1803],
     xlim = c(1,100), type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     main = paste("Functional fit - Yearly average temperature") )
lines(x=1:100, y = MFPCA_temp$meanFunction@.Data[[1]]@X, col = "black", lwd = 3)
legend("topright", legend = "Mean Function", col = "black", lwd = 3)
dev.off()

pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/BasisRep/BasisRep_mean_temp_scen.pdf"), width = 10, height = 11)
par(mfrow =(c(2,2)))
plot(funData_temp[1:434], col = pal(434)[1:434],
     xlim = c(1,100), ylim = c(263,284),
     type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     main = paste("Control") )

plot(funData_temp[bounds[2,1]:bounds[2,2]], col = pal(500)[1:500],
     xlim = c(1,100), ylim = c(263,284),
     type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     main = paste("SSP1-RCP2.6") )

plot(funData_temp[bounds[3,1]:bounds[3,2]], col = pal(500)[1:500],
     xlim = c(1,100), ylim = c(263,284),
     type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     main = paste("SSP3-RCP7.0") )

plot(funData_temp[bounds[4,1]:bounds[4,2]], col = pal(500)[1:500],
     xlim = c(1,100), ylim = c(263,284),
     type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     main = paste("SSP5-RCP8.5") )
mtext("Functional fit - Yearly mean temperature", side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.5)
dev.off()


# Save as png as well
png(paste0("Scripts/Plots/Model/(M)FPCA/png/BasisRep/BasisRep_mean_temp.png"), width = 800, height = 650)
plot(funData_temp, col = pal(1803)[1:1803],
     xlim = c(1,100), type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     main = paste("Functional fit - Yearly average temperature") )
lines(x=1:100, y = MFPCA_temp$meanFunction@.Data[[1]]@X, col = "black", lwd = 3)
legend("topright", legend = "Mean Function", col = "black", lwd = 3)
dev.off()


png(paste0("Scripts/Plots/Model/(M)FPCA/png/BasisRep/BasisRep_mean_temp_scen.png"), width = 1000, height = 1100)
par(mfrow =(c(2,2)))
plot(funData_temp[1:434], col = pal(434)[1:434],
     xlim = c(1,100), ylim = c(263,284),
     type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     main = paste("Control") )

plot(funData_temp[bounds[2,1]:bounds[2,2]], col = pal(500)[1:500],
     xlim = c(1,100), ylim = c(263,284),
     type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     main = paste("SSP1-RCP2.6") )

plot(funData_temp[bounds[3,1]:bounds[3,2]], col = pal(500)[1:500],
     xlim = c(1,100), ylim = c(263,284),
     type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     main = paste("SSP3-RCP7.0") )

plot(funData_temp[bounds[4,1]:bounds[4,2]], col = pal(500)[1:500],
     xlim = c(1,100), ylim = c(263,284),
     type = 'l',
     cex.main = 1.5, cex = 1.3, 
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     main = paste("SSP5-RCP8.5") )
mtext("Functional fit - Yearly mean temperature", side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.5)
dev.off()


## Plot the PCs
# mean temperature
# PC 1
pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/PCs/PC1_mean_temp.pdf"), width = 6.5, height = 6.5)
plot(MFPCA_temp, combined = TRUE, plotPCs = 1, xlim = c(1,100), 
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     cex.main = 1.5, cex = 1.3)
dev.off()

png(paste0("Scripts/Plots/Model/(M)FPCA/png/PCs/PC1_mean_temp.png"), width = 650, height = 650)
plot(MFPCA_temp, combined = TRUE, plotPCs = 1, xlim = c(1,100), 
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     cex.main = 1.5, cex = 1.3)
dev.off()

# PC 2
pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/PCs/PC2_mean_temp.pdf"), width = 6.5, height = 6.5)
plot(MFPCA_temp, combined = TRUE, plotPCs = 2, xlim = c(1,100), 
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     cex.main = 1.5, cex = 1.3)
dev.off()

png(paste0("Scripts/Plots/Model/(M)FPCA/png/PCs/PC2_mean_temp.png"), width = 650, height = 650)
plot(MFPCA_temp, combined = TRUE, plotPCs = 2, xlim = c(1,100), 
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     cex.main = 1.5, cex = 1.3)
dev.off()

# precipitation
# PC 1
pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/PCs/PC1_precip.pdf"), width = 6.5, height = 6.5)
plot(MFPCA_precip, combined = TRUE, plotPCs = 1, xlim = c(1,100), 
     xlab = "Year after Disturbance",
     ylab = bquote("Yearly precipitation in kg/m"^2), 
     cex.main = 1.5, cex = 1.3)
dev.off()

png(paste0("Scripts/Plots/Model/(M)FPCA/png/PCs/PC1_precip.png"), width = 650, height = 650)
plot(MFPCA_precip, combined = TRUE, plotPCs = 1, xlim = c(1,100), 
     xlab = "Year after Disturbance",
     ylab = bquote("Yearly precipitation in kg/m"^2), 
     cex.main = 1.5, cex = 1.3)
dev.off()

# PC 2
pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/PCs/PC2_precip.pdf"), width = 6.5, height = 6.5)
plot(MFPCA_precip, combined = TRUE, plotPCs = 2, xlim = c(1,100), 
     xlab = "Year after Disturbance",
     ylab = bquote("Yearly precipitation in kg/m"^2), 
     cex.main = 1.5, cex = 1.3)
dev.off()

png(paste0("Scripts/Plots/Model/(M)FPCA/png/PCs/PC2_precip.png"), width = 650, height = 650)
plot(MFPCA_precip, combined = TRUE, plotPCs = 2, xlim = c(1,100), 
     xlab = "Year after Disturbance",
     ylab = bquote("Yearly precipitation in kg/m"^2), 
     cex.main = 1.5, cex = 1.3)
dev.off()

## Plot reconstructions
# Save as pdf
# precipitation
pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/Reconstruction/FPCA_reconstruct_precip.pdf"), width = 8, height = 6.5)
plot(MFPCA_precip$fit, col = pal(1803)[1:1803],
     xlim = c(1,100), type = 'l',
     cex.main = 1.8, cex = 1.5,
     xlab = "Year after Disturbance",
     ylab = "Yearly precipitation in kg/m^2", 
     main = paste("Reconstructed fit using 8 PCs"))
dev.off()

iScen = 0
for (scen in scenarios){
  pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/Reconstruction/FPCA_reconstruct_", scen, "_precip.pdf"), width = 8, height = 6.5)
  
  iScen = iScen + 1
  nOb = bounds[iScen,2] - bounds[iScen,1] + 1

  plot(MFPCA_precip$fit, obs = c(bounds[iScen,1]: bounds[iScen,2]), col = pal(nOb)[1:nOb],
       xlim = c(1,100), type = 'l',
       cex.main = 1.8, cex = 1.5,
       xlab = "Year after Disturbance",
       ylab = "Yearly precipitation in kg/m^2", 
       main = paste("Reconstructed fit using 8 PCs -", long_names_scenarios(scen)) )
  dev.off()
  
}

# mean temperature
pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/Reconstruction/FPCA_reconstruct_temp.pdf"), width = 8, height = 6.5)
plot(MFPCA_temp$fit, col = pal(1803)[1:1803],
     xlim = c(1,100), type = 'l',
     cex.main = 1.8, cex = 1.5,
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     main = paste("Reconstructed fit using 8 PCs"))
dev.off()

iScen = 0
for (scen in scenarios){
  pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/Reconstruction/FPCA_reconstruct_", scen, "_temp.pdf"), width = 8, height = 6.5)
  
  iScen = iScen + 1
  nOb = bounds[iScen,2] - bounds[iScen,1] + 1
  
  plot(MFPCA_temp$fit, obs = c(bounds[iScen,1]: bounds[iScen,2]), col = pal(nOb)[1:nOb],
       xlim = c(1,100), type = 'l',
       cex.main = 1.8, cex = 1.5,
       xlab = "Year after Disturbance",
       ylab = "Yearly mean temperature in K", 
       main = paste("Reconstructed fit using 8 PCs -", long_names_scenarios(scen)) )
  dev.off()
  
}
# Save as png
png(paste0("Scripts/Plots/Model/(M)FPCA/png/Reconstruction/FPCA_reconstruct_precip.png"), width = 800, height = 650)
plot(MFPCA_precip$fit, col = pal(1803)[1:1803],
     xlim = c(1,100), type = 'l',
     cex.main = 1.8, cex = 1.5,
     xlab = "Year after Disturbance",
     ylab = "Yearly precipitation in kg/m^2", 
     main = paste("Reconstructed fit using 8 PCs"))
dev.off()

iScen = 0
for (scen in scenarios){
  png(paste0("Scripts/Plots/Model/(M)FPCA/png/Reconstruction/FPCA_reconstruct_", scen, "_precip.png"), width = 800, height = 650)
  
  iScen = iScen + 1
  nOb = bounds[iScen,2] - bounds[iScen,1] + 1
  
  plot(MFPCA_precip$fit, obs = c(bounds[iScen,1]: bounds[iScen,2]), col = pal(nOb)[1:nOb],
       xlim = c(1,100), type = 'l',
       cex.main = 1.8, cex = 1.5,
       xlab = "Year after Disturbance",
       ylab = "Yearly precipitation in kg/m^2", 
       main = paste("Reconstructed fit using 8 PCs -", long_names_scenarios(scen)) )
  dev.off()
  
}
# temperature
png(paste0("Scripts/Plots/Model/(M)FPCA/png/Reconstruction/FPCA_reconstruct_temp.png"), width = 800, height = 650)
plot(MFPCA_temp$fit, col = pal(1803)[1:1803],
     xlim = c(1,100), type = 'l',
     cex.main = 1.8, cex = 1.5,
     xlab = "Year after Disturbance",
     ylab = "Yearly mean temperature in K", 
     main = paste("Reconstructed fit using 8 PCs"))
dev.off()

iScen = 0
for (scen in scenarios){
  png(paste0("Scripts/Plots/Model/(M)FPCA/png/Reconstruction/FPCA_reconstruct_", scen, "_temp.png"), width = 800, height = 650)
  
  iScen = iScen + 1
  nOb = bounds[iScen,2] - bounds[iScen,1] + 1
  
  plot(MFPCA_temp$fit, obs = c(bounds[iScen,1]: bounds[iScen,2]), col = pal(nOb)[1:nOb],
       xlim = c(1,100), type = 'l',
       cex.main = 1.8, cex = 1.5,
       xlab = "Year after Disturbance",
       ylab = "Yearly mean temperature in K", 
       main = paste("Reconstructed fit using 8 PCs -", long_names_scenarios(scen)) )
  dev.off()
}


## Store scores

scores_climate = as.data.frame(MFPCA_climate$scores)
scores_temp = as.data.frame(MFPCA_temp$scores)
scores_precip = as.data.frame(MFPCA_precip$scores)

# Reorder scores according to right order of locations
locs_climate$key = paste(locs_climate$Lon, locs_climate$Lat, locs_climate$Scenario)
locs_disturbed$key = paste(locs_disturbed$Lon, locs_disturbed$Lat, locs_disturbed$Scenario)

index = match(locs_disturbed$key, locs_climate$key)
scores_climate = scores_climate[index,]
scores_temp = scores_temp[index,]
scores_precip = scores_precip[index,]

write.table(scores_climate, "Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_climate.txt")
write.table(scores_temp, "Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_temp.txt")
write.table(scores_precip, "Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_precip.txt")

locs_climate$key = NULL


