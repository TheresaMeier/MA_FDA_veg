################################################################################
############################ Master's Thesis ###################################
################################################################################

##################### Model: MFPCA for climate variables #######################

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
library(rlang)
library(mgcv)

scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")

#######################
pid = 1

# Get target variable: MFPCA scores
plot_data_all = read.table("Scripts/Plots/MFPCA/plot_data/plot_data_all.txt")

# Get locations
locs_disturbed = read.table("Scripts/Plots/MFPCA/plot_data/locs_disturbed.txt")


####################### Prepare climate data for MFPCA #########################
# Load the fun data objects
funData_temp = readRDS("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/funData_temp.rds")
funData_temp_min = readRDS("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/funData_temp_min.rds")
funData_temp_max = readRDS("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/funData_temp_max.rds")
funData_precip = readRDS("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/funData_precip.rds")

# Create multiFunData object
multiFun_climate = multiFunData(funData_temp, funData_temp_min, funData_temp_max, funData_precip)
saveRDS(multiFun_climate, "Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/multiFun_climate.rds")

#################################### Run MFPCA #################################

uniExpansions <- list(list(type = "uFPCA", npc = 5),
                      list(type = "uFPCA", npc = 5),
                      list(type = "uFPCA", npc = 5),
                      list(type = "uFPCA", npc = 5))

MFPCA_climate <- MFPCA_2(multiFun_climate, M = 3, fit = TRUE, uniExpansions = uniExpansions)

saveRDS(MFPCA_climate,"Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/MFPCA_climate.rds")

################################ Plot the PCs ##################################
for (iPC in 1:3){
  pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/PCs/PCs_climate_PC", iPC, ".pdf"), width = 12, height = 5)
  plot.MFPCAfit_2(MFPCA_climate, plotPCs = iPC, combined = TRUE, xlab = paste("Year after Disturbance"), 
                  xlim = c(1,100), ylab = "Values", cex.main = 2.3, cex.lab = 1.6,
                  main = c("Mean temperature", "Minimum temperature", "Maximum temperature", "Precipitation"))
  dev.off()
  
  png(paste0("Scripts/Plots/Model/(M)FPCA/png/PCs/PCs_climate_PC", iPC, ".png"), width = 1200, height = 500)
  plot.MFPCAfit_2(MFPCA_climate, plotPCs = iPC, combined = TRUE, xlab = paste("Year after Disturbance"), 
                  xlim = c(1,100), ylab = "Values", cex.main = 2.3, cex.lab = 1.6,
                  main = c("Mean temperature", "Minimum temperature", "Maximum temperature", "Precipitation"))
  dev.off()
  
}

# Variance that is accounted for
round(100*MFPCA_climate$values/(sum(MFPCA_climate$values)),2)

# Save the scores
scores = as.data.frame(MFPCA_climate$scores)
colnames(scores) = c("PC1", "PC2", "PC3")
write.table(scores, "Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_climate.txt")

