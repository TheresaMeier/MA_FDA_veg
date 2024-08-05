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

run = "_gam"
######
# The data provided can be classified into these four groups:
# - PFT-dependent and time-varying: Nuptake
# - PFT-dependent only: initial_recruitment, recruitment_ten_years, previous_state, time_since_dist
# - time-varying only: tas_yearlymax, tas_yearlymin, tas_yearlsmeam, pr_yearlysum, Nuptake_total
# - Location dependent only: sand_fraction, silt_fraction, clay_fraction, bulkdensity_soil, ph_soil, soilcarbon

############################ Get the prepared data #############################

# Get target variable: MFPCA scores
plot_data_all = read.table("Scripts/Plots/MFPCA/plot_data/plot_data_all.txt")
# Get soil and ecological variables
d_soil = read.table("Scripts/MA_FDA_veg/04_Model/Data/d_soil.txt")
d_eco_wide = read.table("Scripts/MA_FDA_veg/04_Model/Data/d_eco_wide.txt")

# Get FPCA results for temp and precip
scores_temp = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_temp.txt")
scores_precip = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_precip.txt")

combined_d = cbind(plot_data_all[,c(1:2)], d_soil, d_eco_wide[,c(4:18)], scores_temp[,1], scores_precip[,1])

colnames(combined_d)[c(1:2,28:29)] = c("PC1", "PC2", "PC1_temp", "PC1_precip")

################################### Model ######################################

# Produce train-/test- split
set.seed(1)
index = sample.int(nrow(combined_d), size = ceiling(0.8 * nrow(combined_d)))
d_train = combined_d[index,]
d_test = combined_d[-index,]

## Using non-linear effects
gam_1 = gam(list(PC1 ~ s(Lon) + Lat +
                   sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                   Scenario + time_since_dist +
                   initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                   recruitment_ten_years_BNE + recruitment_ten_years_IBS + recruitment_ten_years_otherC + recruitment_ten_years_TeBS + recruitment_ten_years_Tundra +
                   previous_state_BNE + previous_state_IBS + previous_state_otherC + previous_state_TeBS + previous_state_Tundra +
                   s(PC1_temp) + PC1_precip +
                   PC1_temp:PC1_precip
                 , 
                 PC2 ~ s(Lon) + Lat +
                   sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                   Scenario + time_since_dist +
                   initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                   recruitment_ten_years_BNE + recruitment_ten_years_IBS + recruitment_ten_years_otherC + recruitment_ten_years_TeBS + recruitment_ten_years_Tundra +
                   previous_state_BNE + previous_state_IBS + previous_state_otherC + previous_state_TeBS + previous_state_Tundra +
                   s(PC1_temp) + PC1_precip +
                   PC1_temp:PC1_precip
                  ), 
                  data = d_train %>%
                    select(-c(clay_fraction)),
                  family = mvn(d = 2))

summary(gam_1)
saveRDS(gam_1, paste0("Scripts/MA_FDA_veg/04_Model/ModelRDS/gam_1.rds"))

### Model evaluation
multivariate_AIC(gam_1)
# Create residual plots for PC1
pdf(paste0("Scripts/Plots/Model/Evaluation/pdf/residuals", run, ".pdf"), width = 8, height = 8)
par(mfrow = c(2, 2)) 
plot(fitted(gam_1)[, 1], residuals(gam_1)[, 1], main = "Residuals vs. Fitted PC 1 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals_pc1, main = "Normal Q-Q (PC 1)", pch = 16)
qqline(residuals_pc1, col = "red")

plot(fitted(gam_1)[, 2], residuals(gam_1)[, 2], main = "Residuals vs. Fitted PC 2 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals_pc2, main = "Normal Q-Q (PC 2)", pch = 16)
qqline(residuals_pc2, col = "red")
dev.off()

# Save as png as well
png(paste0("Scripts/Plots/Model/Evaluation/png/residuals", run, ".png"), width = 800, height = 800)
par(mfrow = c(2, 2)) 
plot(fitted(gam_1)[, 1],  residuals(gam_1)[, 1], main = "Residuals vs. Fitted PC 1 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals(gam_1)[, 1], main = "Normal Q-Q (PC 1)", pch = 16)
qqline(residuals(gam_1)[, 1], col = "red")

plot(fitted(gam_1)[, 2],  residuals(gam_1)[, 2], main = "Residuals vs. Fitted PC 2 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals(gam_1)[, 2], main = "Normal Q-Q (PC 2)", pch = 16)
qqline(residuals(gam_1)[, 2], col = "red")
dev.off()

pdf(paste0("Scripts/Plots/Model/Evaluation/pdf/gam_1_smooth_effect.pdf"), width = 5, height = 5)
par(mfrow = c(2,2))
plot(gam_1)
dev.off()

png(paste0("Scripts/Plots/Model/Evaluation/png/gam_1_smooth_effect.png"), width = 800, height = 800)
par(mfrow = c(2,2))
plot(gam_1)
dev.off()

## Predict on the test set
predictions = predict(gam_1, newdata = d_test %>% select(-clay_fraction))

MFPCA_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_all_1803_100y.rds")
funData_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_1803.rds")

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

ggsave(paste0("Scripts/Plots/Model/Evaluation/pdf/True_vs_fitted_examples", run, ".pdf"), plot = g, width = 10, height = 8)
ggsave(paste0("Scripts/Plots/Model/Evaluation/png/True_vs_fitted_examples", run, ".png"), plot = g, width = 10, height = 8)

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
