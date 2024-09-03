################################################################################
############################ Master's Thesis ###################################
################################################################################

########################### Model: First Trials ################################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/02_FPCA/functions.R")
source("Scripts/MA_FDA_veg/03_MFPCA/MFPCA/MFPCA_calculation.R")
source("Scripts/MA_FDA_veg/01_Description/utils.R")
source("Scripts/MA_FDA_veg/04_Model/functions_mFLR.R")

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

run = "" # Choose "" (non-transformed target PCs), 
         # "_log" (simple log transformation)
         # "_cmass" (absolute values of cmass)

######
# The data provided can be classified into these four groups:
# - PFT-dependent and time-varying: Nuptake
# - PFT-dependent only: initial_recruitment, recruitment_ten_years, previous_state, time_since_dist
# - time-varying only: tas_yearlymax, tas_yearlymin, tas_yearlsmeam, pr_yearlysum, Nuptake_total
# - Location dependent only: sand_fraction, silt_fraction, clay_fraction, bulkdensity_soil, ph_soil, soilcarbon

############################ Get the prepared data #############################

# Get target variable: MFPCA scores
plot_data_all = read.table("Scripts/Plots/MFPCA/plot_data/plot_data_all.txt")
plot_data_log = read.table("Scripts/Plots/Model/plot_data/plot_data_log.txt")
plot_data_cmass = read.table("Scripts/Plots/Model/plot_data/plot_data_cmass.txt")

# Get soil and ecological variables
d_soil = read.table("Scripts/MA_FDA_veg/04_Model/Data/d_soil.txt")
d_eco_wide = read.table("Scripts/MA_FDA_veg/04_Model/Data/d_eco_wide.txt")

# Get FPCA results for temp and precip
scores_temp = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_temp.txt")
scores_precip = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_precip.txt")

if (run == "_cmass") combined_d = cbind(plot_data_cmass[,c(1,2)], d_soil, d_eco_wide[,c(4:18)], scores_temp[,1], scores_precip[,1])
if (run == "_log") combined_d = cbind(plot_data_log[,c(1,2)], d_soil, d_eco_wide[,c(4:18)], scores_temp[,1], scores_precip[,1])
if (run == "") combined_d = cbind(plot_data_all[,c(1:2)], d_soil, d_eco_wide[,c(4:18)], scores_temp[,1], scores_precip[,1])

colnames(combined_d)[c(1:2,28:29)] = c("PC1", "PC2", "PC1_temp", "PC1_precip")

# Filter for North-America
# locs_disturbed = read.table("Scripts/Plots/MFPCA/plot_data/locs_disturbed.txt")
# index_NA = locs_disturbed$Lon < -40
# combined_d = combined_d[index_NA,]
################################### Model ######################################

# Produce train-/test- split
set.seed(1)
index = sample.int(nrow(combined_d), size = ceiling(0.8 * nrow(combined_d)))
d_train = combined_d[index,]
d_test = combined_d[-index,]

# Fit the model
lm_1 = lm(cbind(PC1, PC2) ~ ., data = d_train %>%
            dplyr::select(-c(clay_fraction)))
summary(lm_1)
#saveRDS(lm_1, paste0("Scripts/MA_FDA_veg/04_Model/ModelRDS/lm_1", run, ".rds"))

### Model evaluation
# AIC
multivariate_AIC(lm_1)

residuals_pc1 <- residuals(lm_1)[, 1]
residuals_pc2 <- residuals(lm_1)[, 2]

# Create residual plots for PC1
pdf(paste0("Scripts/Plots/Model/Evaluation/pdf/residuals", run, ".pdf"), width = 8, height = 8)
par(mfrow = c(2, 2)) 
plot(fitted(lm_1)[, 1], d_train$PC1, main = "Residuals vs. Fitted PC 1 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals_pc1, main = "Normal Q-Q (PC 1)", pch = 16)
qqline(residuals_pc1, col = "red")

plot(fitted(lm_1)[, 2], residuals_pc2, main = "Residuals vs. Fitted PC 2 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals_pc2, main = "Normal Q-Q (PC 2)", pch = 16)
qqline(residuals_pc2, col = "red")
dev.off()

# Save as png as well
png(paste0("Scripts/Plots/Model/Evaluation/png/residuals", run, ".png"), width = 800, height = 800)
par(mfrow = c(2, 2)) 
plot(fitted(lm_1)[, 1], residuals_pc1, main = "Residuals vs. Fitted PC 1 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals_pc1, main = "Normal Q-Q (PC 1)", pch = 16)
qqline(residuals_pc1, col = "red")

plot(fitted(lm_1)[, 2], residuals_pc2, main = "Residuals vs. Fitted PC 2 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals_pc2, main = "Normal Q-Q (PC 2)", pch = 16)
qqline(residuals_pc2, col = "red")
dev.off()

## Predict on the test set
predictions = predict(lm_1, newdata = d_test %>% select(-clay_fraction))

if(run == "_cmass"){
  MFPCA_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_all_cmass.rds")
  funData_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_cmass.rds")
}

if(run == "_log"){
  MFPCA_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_log.rds")
  funData_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_pft_log.rds")
}

if(run == ""){
  MFPCA_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_all_1803_100y.rds")
  funData_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_1803.rds")
}

uniExpansions <- list(list(type = "uFPCA", npc = 10),
                      list(type = "uFPCA", npc = 10),
                      list(type = "uFPCA", npc = 10),
                      list(type = "uFPCA", npc = 10),
                      list(type = "uFPCA", npc = 10))


MFPCA_2PCs <- MFPCA(funData_all, M = 2, fit = TRUE, uniExpansions = uniExpansions)
saveRDS(MFPCA_2PCs, paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_all_1803_100y_2PCs", run, ".rds"))

pred_MFPCA = predict(MFPCA_2PCs, scores = predictions)
MFPCA_test = MFPCA_2PCs[-index]
funData_test = funData_all[-index]

########################### Plot example functional fits #######################
# Create a data frame from the provided data

j=0
dat = c()
for (iCurve in c(1,88,166,258)){
  j = j+1
  funData_iCurve = funData_test[iCurve]
  pred_iCurve = pred_MFPCA[iCurve]
  for (i in 1:5){
    pft = rep(long_names_pfts(tolower(pfts[i])),200)
    age = rep(1:100,2)
    scenario = rep(long_names_scenarios(scenarios[j]),200)
    values_true = funData_iCurve@.Data[[i]]@X[1,]
    values_pred = pred_iCurve@.Data[[i]]@X[1,]
    values = c(values_true, values_pred)
    
    if (run == "_log") values = exp(values) - 1

    model = c(rep("True",100), rep("Predicted", 100))
    
    dat = rbind(dat,cbind(pft,age,scenario,values,model))
  }
}

dat = as.data.frame(dat) %>%
  mutate(age = as.numeric(age),
         values = as.numeric(values),
         pft = as.factor(pft),
         scenario = as.factor(scenario),
         model = factor(model, levels = c("True", "Predicted")))


# Create the plot
g = ggplot(dat) + 
      geom_line(data = dat, aes(x = age, y = values, col = pft, linetype = model), lwd = 1) +
      facet_grid(rows = vars(scenario)) +
      labs(
        x = "Year after disturbance",
        y = "Share of aboveground carbon",
        title = "True functional fit and predicted one for example grid cells",
        color = "Dominant vegetation",
        linetype = "Model"
      ) +
      theme_bw() +
      theme(
        text = element_text(size = 12),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.box = "vertical"
      ) +
      scale_color_manual(values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                                    "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73"))

if(run == "_cmass") g = g+ labs(y = bquote("Aboveground carbon in kg/m"^2))

ggsave(paste0("Scripts/Plots/Model/Evaluation/pdf/True_vs_fitted_examples", run, ".pdf"), plot = g, width = 10, height = 8)
ggsave(paste0("Scripts/Plots/Model/Evaluation/png/True_vs_fitted_examples", run, ".png"), plot = g, width = 10, height = 8)



######## For cmass data: transform the data to shares
if (run == "_cmass"){
  
  for (mod in c("pred_MFPCA", "funData_test")){
    dat_mod = get(mod)
    matrices <- lapply(dat_mod@.Data, function(obj) {
      x <- obj@X
      x[x < 0] <- 0
      return(x)
    })
    
    total_sum <- Reduce(`+`, matrices)
    
    relative_shares <- lapply(matrices, function(mat) {
      rel_share <- mat / total_sum
      rel_share[is.nan(rel_share)] <- 0  # Set NaN values to 0
      return(rel_share)
    })
    
    dat_mod_shares <- lapply(seq_along(relative_shares), function(i) {
      new_obj <- dat_mod@.Data[[i]]  # Assuming each element in pred_MFPCA@.Data is an object with properties
      new_obj@X <- relative_shares[[i]]
      return(new_obj)
    })
    
    assign(paste0(mod, "_shares"), multiFunData(dat_mod_shares))
  }
 
  j=0
  dat = c()
  for (iCurve in c(1,88,166,258)){
    j = j+1
    funData_iCurve = funData_test_shares[iCurve]
    pred_iCurve = pred_MFPCA_shares[iCurve]
    for (i in 1:5){
      pft = rep(long_names_pfts(tolower(pfts[i])),200)
      age = rep(1:100,2)
      scenario = rep(long_names_scenarios(scenarios[j]),200)
      values_true = funData_iCurve@.Data[[i]]@X[1,]
      values_pred = pred_iCurve@.Data[[i]]@X[1,]
      values = c(values_true, values_pred)
      
      if (run == "_log") values = exp(values)
      
      model = c(rep("True",100), rep("Predicted", 100))
      
      dat = rbind(dat,cbind(pft,age,scenario,values,model))
    }
  }
  
  dat = as.data.frame(dat) %>%
    mutate(age = as.numeric(age),
           values = as.numeric(values),
           pft = as.factor(pft),
           scenario = as.factor(scenario),
           model = factor(model, levels = c("True", "Predicted")))
  
  
  # Create the plot
  ggplot(dat) + 
    geom_line(data = dat, aes(x = age, y = values, col = pft, linetype = model), lwd = 1) +
    facet_grid(rows = vars(scenario)) +
    labs(
      x = "Year after disturbance",
      y = "Share of aboveground carbon",
      title = "True functional fit and predicted one for example grid cells",
      color = "Dominant vegetation",
      linetype = "Model"
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.box = "vertical"
    ) +
    scale_color_manual(values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                                  "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73"))
  
  ggsave(paste0("Scripts/Plots/Model/Evaluation/pdf/True_vs_fitted_examples_cmass_shares.pdf"), width = 10, height = 8)
  ggsave(paste0("Scripts/Plots/Model/Evaluation/png/True_vs_fitted_examples_cmass_shares.png"), width = 10, height = 8)
}



################### Plot predicted scores vs. true scores ######################
pdf(paste0("Scripts/Plots/Model/Evaluation/pdf/True_vs_fitted_test", run, ".pdf"), width = 10, height = 5)
par(mfrow=c(1,2))
plot(MFPCA_2PCs$scores[-index,1], predictions[,1],
     pch = 16,
     xlab = "Predicted PC 1 scores",
     ylab = "True PC 1 scores")
abline(a = 0, b = 1, col = "red", lwd = 2)
plot(predictions[,2], MFPCA_2PCs$scores[-index,2],
     pch = 16,
     xlab = "Predicted PC 2 scores",
     ylab = "True PC 2 scores")
abline(a = 0, b = 1, col = "red", lwd = 2)
mtext("True vs. fitted PC scores for unseen data", side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
dev.off()

# Save as png as well
png(paste0("Scripts/Plots/Model/Evaluation/png/True_vs_fitted_test", run, ".png"), width = 1000, height = 500)
par(mfrow=c(1,2))
plot(MFPCA_2PCs$scores[-index,1], predictions[,1],
     pch = 16,
     xlab = "Predicted PC 1 scores",
     ylab = "True PC 1 scores")
abline(a = 0, b = 1, col = "red", lwd = 2)
plot(predictions[,2], MFPCA_2PCs$scores[-index,2],
     pch = 16,
     xlab = "Predicted PC 2 scores",
     ylab = "True PC 2 scores")
abline(a = 0, b = 1, col = "red", lwd = 2)
mtext("True vs. fitted PC scores for unseen data", side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
dev.off()

######################## Check for non-linear effects ##########################

# Location
png(paste0("Scripts/Plots/Model/Evaluation/png/Resid_location", run, ".png"), width = 800, height = 800)
par(mfrow = c(2,2))
plot(combined_d$Lon[index], resid(lm_1)[,1], xlab = "Longitude", ylab = "Model's residuals", main = "PC 1", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$Lon[index], resid(lm_1)[,2], xlab = "Longitude", ylab = "Model's residuals", main = "PC 2", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$Lat[index], resid(lm_1)[,1], xlab = "Latitude", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$Lat[index], resid(lm_1)[,2], xlab = "Latitude", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
dev.off()

# Of course there are no values between -50 and 0 Longitude, since this area is in the sea --> no grid cells at all. Problem?

# Soil composition
png(paste0("Scripts/Plots/Model/Evaluation/png/Resid_soil_composition", run, ".png"), width = 800, height = 800)
par(mfrow = c(2,2))
plot(combined_d$sand_fraction[index], resid(lm_1)[,1], xlab = "Sand fraction", ylab = "Model's residuals", main = "PC 1", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$sand_fraction[index], resid(lm_1)[,2], xlab = "Sand fraction", ylab = "Model's residuals", main = "PC 2", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$silt_fraction[index], resid(lm_1)[,1], xlab = "Silt fraction", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$silt_fraction[index], resid(lm_1)[,2], xlab = "Silt fraction", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
dev.off()

# Soil attributes
png(paste0("Scripts/Plots/Model/Evaluation/png/Resid_soil_attributes", run, ".png"), width = 800, height = 1200)
par(mfrow = c(3,2))
plot(combined_d$bulkdensity_soil[index], resid(lm_1)[,1], xlab = "Bulk density", ylab = "Model's residuals", main = "PC 1", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$bulkdensity_soil[index], resid(lm_1)[,2], xlab = "Bulk density", ylab = "Model's residuals", main = "PC 2", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$ph_soil[index], resid(lm_1)[,1], xlab = "pH in water", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$ph_soil[index], resid(lm_1)[,2], xlab = "pH in water", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$soilcarbon[index], resid(lm_1)[,1], xlab = "Organic carbon content", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$soilcarbon[index], resid(lm_1)[,2], xlab = "Organic carbon content", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
dev.off()

# Initial recruitment
png(paste0("Scripts/Plots/Model/Evaluation/png/Resid_initial_recruitment", run, ".png"), width = 800, height = 2000)
par(mfrow = c(5,2))
plot(combined_d$initial_recruitment_BNE[index], resid(lm_1)[,1], xlab = "Initial recruitment - BNE", ylab = "Model's residuals", main = "PC 1", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$initial_recruitment_BNE[index], resid(lm_1)[,2], xlab = "Initial recruitment - BNE", ylab = "Model's residuals", main = "PC 2", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$initial_recruitment_IBS[index], resid(lm_1)[,1], xlab = "Initial recruitment - IBS", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$initial_recruitment_IBS[index], resid(lm_1)[,2], xlab = "Initial recruitment - IBS", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$initial_recruitment_otherC[index], resid(lm_1)[,1], xlab = "Initial recruitment - otherC", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$initial_recruitment_otherC[index], resid(lm_1)[,2], xlab = "Initial recruitment - otherC", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$initial_recruitment_TeBS[index], resid(lm_1)[,1], xlab = "Initial recruitment - TeBS", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$initial_recruitment_TeBS[index], resid(lm_1)[,2], xlab = "Initial recruitment - TeBS", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$initial_recruitment_Tundra[index], resid(lm_1)[,1], xlab = "Initial recruitment - Tundra", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$initial_recruitment_Tundra[index], resid(lm_1)[,2], xlab = "Initial recruitment - Tundra", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
dev.off()

# Recruitment 10 years
png(paste0("Scripts/Plots/Model/Evaluation/png/Resid_recruitment_ten_years", run, ".png"), width = 800, height = 2000)
par(mfrow = c(5,2))
plot(combined_d$recruitment_ten_years_BNE[index], resid(lm_1)[,1], xlab = "Recruitment ten years - BNE", ylab = "Model's residuals", main = "PC 1", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$recruitment_ten_years_BNE[index], resid(lm_1)[,2], xlab = "Recruitment ten years - BNE", ylab = "Model's residuals", main = "PC 2", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$recruitment_ten_years_IBS[index], resid(lm_1)[,1], xlab = "Recruitment ten years - IBS", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$recruitment_ten_years_IBS[index], resid(lm_1)[,2], xlab = "Recruitment ten years - IBS", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$recruitment_ten_years_otherC[index], resid(lm_1)[,1], xlab = "Recruitment ten years - otherC", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$recruitment_ten_years_otherC[index], resid(lm_1)[,2], xlab = "Recruitment ten years - otherC", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$recruitment_ten_years_TeBS[index], resid(lm_1)[,1], xlab = "Recruitment ten years - TeBS", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$recruitment_ten_years_TeBS[index], resid(lm_1)[,2], xlab = "Recruitment ten years - TeBS", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$recruitment_ten_years_Tundra[index], resid(lm_1)[,1], xlab = "Recruitment ten years - Tundra", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$recruitment_ten_years_Tundra[index], resid(lm_1)[,2], xlab = "Recruitment ten years - Tundra", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
dev.off()

# Previous state
png(paste0("Scripts/Plots/Model/Evaluation/png/Resid_previous_state", run, ".png"), width = 800, height = 2000)
par(mfrow = c(5,2))
plot(combined_d$previous_state_BNE[index], resid(lm_1)[,1], xlab = "Previous state - BNE", ylab = "Model's residuals", main = "PC 1", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$previous_state_BNE[index], resid(lm_1)[,2], xlab = "Previous state - BNE", ylab = "Model's residuals", main = "PC 2", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$previous_state_IBS[index], resid(lm_1)[,1], xlab = "Previous state - IBS", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$previous_state_IBS[index], resid(lm_1)[,2], xlab = "Previous state - IBS", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$previous_state_otherC[index], resid(lm_1)[,1], xlab = "Previous state - otherC", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$previous_state_otherC[index], resid(lm_1)[,2], xlab = "Previous state - otherC", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$previous_state_TeBS[index], resid(lm_1)[,1], xlab = "Previous state - TeBS", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$previous_state_TeBS[index], resid(lm_1)[,2], xlab = "Previous state - TeBS", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$previous_state_Tundra[index], resid(lm_1)[,1], xlab = "Previous state - Tundra", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$previous_state_Tundra[index], resid(lm_1)[,2], xlab = "Previous state - Tundra", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
dev.off()

# Climate PCs
png(paste0("Scripts/Plots/Model/Evaluation/png/Resid_climate", run, ".png"), width = 800, height = 800)
par(mfrow = c(2,2))
plot(combined_d$PC1_temp[index], resid(lm_1)[,1], xlab = "PC 1 - Average annual temperature", ylab = "Model's residuals", main = "PC 1", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$PC1_temp[index], resid(lm_1)[,2], xlab = "PC 1 - Average annual temperature", ylab = "Model's residuals", main = "PC 2", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$PC1_precip[index], resid(lm_1)[,1], xlab = "PC 1 - Annual precipitation", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
plot(combined_d$PC1_precip[index], resid(lm_1)[,2], xlab = "PC 1 - Annual precipitation", ylab = "Model's residuals", pch = 16)
abline(h = 0, col = "red")
dev.off()


# Plot spacial residuals

dat_locs = d_train[,c(3,4)]
dat_locs$resid_PC1 = residuals(lm_1)[,1]
dat_locs$resid_PC2 = residuals(lm_1)[,2]

library(sf)
df_sf <- st_as_sf(dat_locs, coords = c("Lon", "Lat"), crs = 4326)
ggplot(data = df_sf) +
       geom_sf(aes(color = resid_PC1)) +
       scale_color_viridis_c() + # Optional: for a nice color scale
       theme_minimal() +
       labs(title = "Map of Residuals",
                      color = "Residuals") +
       theme(legend.position = "right")

custom_color_scale <- scale_color_gradient2(
  low = "red", mid = "white", high = "blue", 
  midpoint = 0, space = "Lab", 
  na.value = "grey50", guide = "colourbar", 
  aesthetics = "color"
)

# Plotting
ggplot(data = df_sf) +
  geom_sf(aes(color = resid_PC1)) +
  custom_color_scale +
  theme_minimal() +
  labs(title = "Map of Residuals - PC 1",
       color = "Residuals") +
  theme(legend.position = "right")

ggsave(paste0("Scripts/Plots/Model/Evaluation/pdf/Map_residuals_PC1.pdf"), width = 10, height = 8)
ggsave(paste0("Scripts/Plots/Model/Evaluation/png/Map_residuals_PC1.png"), width = 10, height = 8)

ggplot(data = df_sf) +
  geom_sf(aes(color = resid_PC2)) +
  custom_color_scale +
  theme_minimal() +
  labs(title = "Map of Residuals - PC 2",
       color = "Residuals") +
  theme(legend.position = "right")


ggsave(paste0("Scripts/Plots/Model/Evaluation/pdf/Map_residuals_PC2.pdf"), width = 10, height = 8)
ggsave(paste0("Scripts/Plots/Model/Evaluation/png/Map_residuals_PC2.png"), width = 10, height = 8)

