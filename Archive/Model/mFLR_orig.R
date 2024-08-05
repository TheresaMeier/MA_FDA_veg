################################################################################
############################ Master's Thesis ###################################
################################################################################

###################### Model: Original data as input ###########################

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

scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")
pfts = c("BNE", "IBS", "otherC", "TeBS","Tundra")

run = "_orig" # Choose "" (non-transformed target PCs), 

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

# Get target variable: Try with original data for different time points and PFT BNE
funData_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_1803.rds")
combined_d = cbind(funData_all@.Data[[1]]@X[,20], funData_all@.Data[[1]]@X[,40],
                   funData_all@.Data[[1]]@X[,60], funData_all@.Data[[1]]@X[,80],
                   funData_all@.Data[[1]]@X[,100],
                   d_soil, d_eco_wide[,c(4:18)], scores_temp[,1], scores_precip[,1])

colnames(combined_d)[c(1:5,31:32)] = c("Orig_20","Orig_40","Orig_60","Orig_80",
                                       "Orig_100", "PC1_temp", "PC1_precip")

################################### Model ######################################

# Produce train-/test- split
set.seed(1)
index = sample.int(nrow(combined_d), size = ceiling(0.8 * nrow(combined_d)))
d_train = combined_d[index,]
d_test = combined_d[-index,]

# Fit the model
lm_1 = lm(cbind(Orig_20, Orig_40, Orig_60, Orig_80, Orig_100) ~ ., data = d_train %>%
            select(-c(clay_fraction)))
summary(lm_1)
#saveRDS(lm_1, paste0("Scripts/MA_FDA_veg/04_Model/ModelRDS/lm_1", run, ".rds"))

### Model evaluation

# Create residual plots for PC1
pdf(paste0("Scripts/Plots/Model/Evaluation/pdf/residuals", run, ".pdf"), width = 8, height = 10)
par(mfrow = c(5, 2)) 
plot(fitted(lm_1)[, 1], residuals(lm_1)[, 1], main = "Residuals vs. Fitted - 20 years post-disturbance",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals(lm_1)[, 1], main = "Normal Q-Q (20 years)", pch = 16)
qqline(residuals(lm_1)[, 1], col = "red")

plot(fitted(lm_1)[, 2], residuals(lm_1)[, 2], main = "Residuals vs. Fitted - 40 years post-disturbance",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals(lm_1)[, 2], main = "Normal Q-Q (40 years)", pch = 16)
qqline(residuals(lm_1)[, 2], col = "red")

plot(fitted(lm_1)[, 3], residuals(lm_1)[, 3], main = "Residuals vs. Fitted - 60 years post-disturbance",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals(lm_1)[, 3], main = "Normal Q-Q (60 years)", pch = 16)
qqline(residuals(lm_1)[, 3], col = "red")

plot(fitted(lm_1)[, 4], residuals(lm_1)[, 4], main = "Residuals vs. Fitted - 80 years post-disturbance",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals(lm_1)[, 4], main = "Normal Q-Q (80 years)", pch = 16)
qqline(residuals(lm_1)[, 4], col = "red")

plot(fitted(lm_1)[, 5], residuals(lm_1)[, 5], main = "Residuals vs. Fitted - 100 years post-disturbance",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals(lm_1)[, 5], main = "Normal Q-Q (100 years)", pch = 16)
qqline(residuals(lm_1)[, 5], col = "red")
dev.off()

# Save as png as well
png(paste0("Scripts/Plots/Model/Evaluation/png/residuals", run, ".png"), width = 800, height = 1000)
par(mfrow = c(5, 2)) 
plot(fitted(lm_1)[, 1], residuals(lm_1)[, 1], main = "Residuals vs. Fitted - 20 years post-disturbance",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals(lm_1)[, 1], main = "Normal Q-Q (20 years)", pch = 16)
qqline(residuals(lm_1)[, 1], col = "red")

plot(fitted(lm_1)[, 2], residuals(lm_1)[, 2], main = "Residuals vs. Fitted - 40 years post-disturbance",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals(lm_1)[, 2], main = "Normal Q-Q (40 years)", pch = 16)
qqline(residuals(lm_1)[, 2], col = "red")

plot(fitted(lm_1)[, 3], residuals(lm_1)[, 3], main = "Residuals vs. Fitted - 60 years post-disturbance",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals(lm_1)[, 3], main = "Normal Q-Q (60 years)", pch = 16)
qqline(residuals(lm_1)[, 3], col = "red")

plot(fitted(lm_1)[, 4], residuals(lm_1)[, 4], main = "Residuals vs. Fitted - 80 years post-disturbance",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals(lm_1)[, 4], main = "Normal Q-Q (80 years)", pch = 16)
qqline(residuals(lm_1)[, 4], col = "red")

plot(fitted(lm_1)[, 5], residuals(lm_1)[, 5], main = "Residuals vs. Fitted - 100 years post-disturbance",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals(lm_1)[, 5], main = "Normal Q-Q (100 years)", pch = 16)
qqline(residuals(lm_1)[, 5], col = "red")
dev.off()

