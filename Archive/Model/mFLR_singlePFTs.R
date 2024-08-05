################################################################################
############################ Master's Thesis ###################################
################################################################################

######################## Model: Try different PFTs #############################

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
library(RColorBrewer)
library(mgcv)
library(purrr)

scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")
pfts = c("BNE", "IBS", "otherC", "TeBS","Tundra")

######
# The data provided can be classified into these four groups:
# - PFT-dependent and time-varying: Nuptake
# - PFT-dependent only: initial_recruitment, recruitment_ten_years, previous_state, time_since_dist
# - time-varying only: tas_yearlymax, tas_yearlymin, tas_yearlsmeam, pr_yearlysum, Nuptake_total
# - Location dependent only: sand_fraction, silt_fraction, clay_fraction, bulkdensity_soil, ph_soil, soilcarbon

############################ Get the prepared data #############################
# Get soil and ecological variables
d_soil = read.table("Scripts/MA_FDA_veg/04_Model/Data/d_soil.txt")
d_eco_wide = read.table("Scripts/MA_FDA_veg/04_Model/Data/d_eco_wide.txt")

# Get FPCA results for temp and precip
scores_temp = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_temp.txt")
scores_precip = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_precip.txt")

for (iPFT in c(1:5)){
  run = paste0("_",pfts[iPFT])
  
  print(paste0("Start with PFT ", paste(pfts[iPFT], collapse = ", ")))
  # multiFun_pft = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_1803.rds")
  multiFun_pft = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_cmass.rds")
  #for (i in iPFT:5)multiFun_pft@.Data <- multiFun_pft@.Data[-(2:iPFT)]
  
  multiFun_pft@.Data <- multiFun_pft@.Data[-setdiff(c(1:5), iPFT)]
  # uniExpansions <- replicate(iPFT-1, list(type = "uFPCA", npc = 10), simplify = FALSE)
  uniExpansions <- list(list(type = "uFPCA", npc = 10))
  
  MFPCA_iPFT <- MFPCA(multiFun_pft, M = 10, fit = TRUE, uniExpansions = uniExpansions)
  
  plot_data_all = as.data.frame(MFPCA_iPFT$scores) %>% 
    mutate(scenario = c(rep("Control", 434),rep("SSP1-RCP2.6", 442),rep("SSP3-RCP7.0", 462),rep("SSP5-RCP8.5", 465)))
  
  combined_d = cbind(plot_data_all[,c(1:2)], d_soil, d_eco_wide[,c(4:18)], scores_temp[,1], scores_precip[,1])
  
  colnames(combined_d)[c(1:2,28:29)] = c("PC1", "PC2", "PC1_temp", "PC1_precip")
  
  # Produce train-/test- split
  set.seed(1)
  index = sample.int(nrow(combined_d), size = ceiling(0.8 * nrow(combined_d)))
  d_train = combined_d[index,]
  d_test = combined_d[-index,]
  
  # Fit the model
  lm_1 = lm(cbind(PC1, PC2) ~ ., data = d_train %>%
              select(-c(clay_fraction)))
  summary(lm_1)
  saveRDS(lm_1, paste0("Scripts/MA_FDA_veg/04_Model/ModelRDS/lm_1", run, ".rds"))
  
  ### Model evaluation
  residuals_pc1 <- residuals(lm_1)[, 1]
  residuals_pc2 <- residuals(lm_1)[, 2]
  
  # Create residual plots for PC1
  pdf(paste0("Scripts/Plots/Model/Evaluation/pdf/residuals", run, ".pdf"), width = 8, height = 8)
  par(mfrow = c(2, 2)) 
  plot(fitted(lm_1)[, 1], residuals_pc1, main = paste0("Residuals vs. Fitted PC 1 scores - PFTs ", paste(pfts[iPFT], collapse = ", ")),
       xlab = "Fitted values", ylab = "Residuals", pch = 16)
  abline(h = 0, col = "red")
  
  qqnorm(residuals_pc1, main = "Normal Q-Q (PC 1)", pch = 16)
  qqline(residuals_pc1, col = "red")
  
  plot(fitted(lm_1)[, 2], residuals_pc2, main = paste0("Residuals vs. Fitted PC 2 scores - PFTs ", paste(pfts[iPFT], collapse = ", ")),
       xlab = "Fitted values", ylab = "Residuals", pch = 16)
  abline(h = 0, col = "red")
  
  qqnorm(residuals_pc2, main = "Normal Q-Q (PC 2)", pch = 16)
  qqline(residuals_pc2, col = "red")
  dev.off()
  
  # Save as png as well
  png(paste0("Scripts/Plots/Model/Evaluation/png/residuals", run, ".png"), width = 800, height = 800)
  par(mfrow = c(2, 2)) 
  plot(fitted(lm_1)[, 1], residuals_pc1, main = paste0("Residuals vs. Fitted PC 1 scores - PFTs ", paste(pfts[iPFT], collapse = ", ")),
       xlab = "Fitted values", ylab = "Residuals", pch = 16)
  abline(h = 0, col = "red")
  
  qqnorm(residuals_pc1, main = "Normal Q-Q (PC 1)", pch = 16)
  qqline(residuals_pc1, col = "red")
  
  plot(fitted(lm_1)[, 2], residuals_pc2, main = paste0("Residuals vs. Fitted PC 2 scores - PFTs ", paste(pfts[iPFT], collapse = ", ")),
       xlab = "Fitted values", ylab = "Residuals", pch = 16)
  abline(h = 0, col = "red")
  
  qqnorm(residuals_pc2, main = "Normal Q-Q (PC 2)", pch = 16)
  qqline(residuals_pc2, col = "red")
  dev.off()
}


