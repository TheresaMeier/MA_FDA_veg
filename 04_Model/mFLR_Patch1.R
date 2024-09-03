################################################################################
############################ Master's Thesis ###################################
################################################################################

########################### Model: Final Version ###############################

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
run = "_gam_final"

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

# Get MFPCA results for temp (min, max, mean) and precip
scores_climate = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_climate.txt")

# Get FPCA results for temp, precip and nuptake_total
scores_temp = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_temp.txt")
scores_temp_min = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_temp_min.txt")
scores_temp_max = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_temp_max.txt")
scores_precip = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_precip.txt")
scores_nuptake_total = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_nuptake_total.txt")

# Nuptake PFT
scores_nuptake_BNE = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_nuptake_BNE.txt")
scores_nuptake_IBS = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_nuptake_IBS.txt")
scores_nuptake_otherC = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_nuptake_otherC.txt")
scores_nuptake_TeBS = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_nuptake_TeBS.txt")
scores_nuptake_Tundra = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_nuptake_Tundra.txt")


# Combine data set
combined_d = cbind(plot_data_all[,c(1:2)], d_soil, d_eco_wide[,c(4:18)], 
                   scores_climate[,1],
                   scores_temp[,1], scores_temp_min[,1], scores_temp_max[,1],
                   scores_precip[,1], scores_nuptake_total[,1], scores_nuptake_BNE[,1],
                   scores_nuptake_IBS[,1], scores_nuptake_otherC[,1], scores_nuptake_TeBS[,1], scores_nuptake_Tundra[,1])  %>% 
  mutate(Scenario = as.factor(Scenario))

colnames(combined_d)[c(1:2,28:38)] = c("PC1", "PC2", "PC1_climate", "PC1_temp", "PC1_temp_min", "PC1_temp_max",
                                       "PC1_precip", "PC1_nuptake_total",
                                       "PC1_nuptake_BNE", "PC1_nuptake_IBS", "PC1_nuptake_otherC",
                                       "PC1_nuptake_TeBS", "PC1_nuptake_Tundra")

################################### Model ######################################

# Produce train-/test- split
set.seed(1)
index = sample.int(nrow(combined_d), size = ceiling(0.8 * nrow(combined_d)))
d_train = combined_d[index,]
d_test = combined_d[-index,]

# Fit the model using non-linear effects
gam_final = gam(list(PC1 ~ s(Lon,Lat) +
                   sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                   Scenario + time_since_dist +
                   initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                   PC1_climate, 
                 PC2 ~ s(Lon,Lat) +
                   sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                   Scenario + time_since_dist +
                   initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                   PC1_climate), 
      data = d_train,
      family = mvn(d = 2))

summary(gam_final)
saveRDS(gam_final, paste0("Scripts/MA_FDA_veg/04_Model/ModelRDS/gam_final.rds"))

################################################################################
############################# Model evaluation #################################

############################# Residual plots ###################################
pdf(paste0("Scripts/Plots/Model/Evaluation/pdf/residuals", run, ".pdf"), width = 8, height = 8)
par(mfrow = c(2, 2)) 
plot(fitted(gam_final)[, 1], residuals(gam_final)[, 1], main = "Residuals vs. Fitted PC 1 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16, cex.main = 1.5, cex.lab = 1.3)
abline(h = 0, col = "red", lwd = 2)

qqnorm(residuals(gam_final)[, 1], main = "Normal Q-Q (PC 1)", pch = 16, cex.main = 1.5, cex.lab = 1.3)
qqline(residuals(gam_final)[, 1], col = "red", lwd = 2)

plot(fitted(gam_final)[, 2], residuals(gam_final)[, 2], main = "Residuals vs. Fitted PC 2 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16, cex.main = 1.5, cex.lab = 1.3)
abline(h = 0, col = "red", lwd = 2)

qqnorm(residuals(gam_final)[, 2], main = "Normal Q-Q (PC 2)", pch = 16, cex.main = 1.5, cex.lab = 1.3)
qqline(residuals(gam_final)[, 2], col = "red", lwd = 2)
dev.off()

# Save as png as well
png(paste0("Scripts/Plots/Model/Evaluation/png/residuals", run, ".png"), width = 800, height = 800)
par(mfrow = c(2, 2)) 
plot(fitted(gam_final)[, 1],  residuals(gam_final)[, 1], main = "Residuals vs. Fitted PC 1 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals(gam_final)[, 1], main = "Normal Q-Q (PC 1)", pch = 16)
qqline(residuals(gam_final)[, 1], col = "red")

plot(fitted(gam_final)[, 2],  residuals(gam_final)[, 2], main = "Residuals vs. Fitted PC 2 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals(gam_final)[, 2], main = "Normal Q-Q (PC 2)", pch = 16)
qqline(residuals(gam_final)[, 2], col = "red")
dev.off()

#################### True vs. fitted on training data set ######################
pdf(paste0("Scripts/Plots/Model/Evaluation/pdf/True_vs_fitted_train", run, ".pdf"), width = 7, height = 3.5)
par(mfrow=c(1,2))
plot(fitted(gam_final)[,1], d_train$PC1, 
     pch = 16,
     xlab = "Predicted PC 1 scores",
     ylab = "True PC 1 scores")
abline(a = 0, b = 1, col = "red", lwd = 2)
plot(fitted(gam_final)[,2], d_train$PC2, 
     pch = 16,
     xlab = "Predicted PC 2 scores",
     ylab = "True PC 2 scores")
abline(a = 0, b = 1, col = "red", lwd = 2)
mtext("True vs. fitted PC scores for the training data", side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
dev.off()

# Save as png as well
png(paste0("Scripts/Plots/Model/Evaluation/png/True_vs_fitted_train", run, ".png"), width = 1000, height = 500)
par(mfrow=c(1,2))
plot(fitted(gam_final)[,1], d_train$PC1, 
     pch = 16,
     xlab = "Predicted PC 1 scores",
     ylab = "True PC 1 scores")
abline(a = 0, b = 1, col = "red", lwd = 2)
plot(fitted(gam_final)[,2], d_train$PC2, 
     pch = 16,
     xlab = "Predicted PC 2 scores",
     ylab = "True PC 2 scores")
abline(a = 0, b = 1, col = "red", lwd = 2)
mtext("True vs. fitted PC scores for the training data", side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
dev.off()

######################## Predict on the test set ###############################
predictions = predict(gam_final, newdata = d_test)

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
for (iCurve in c(2,89,167,259)){
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

ggsave(paste0("Scripts/Plots/Model/Evaluation/pdf/True_vs_fitted_examples", run, ".pdf"), width = 10, height = 8)
ggsave(paste0("Scripts/Plots/Model/Evaluation/png/True_vs_fitted_examples", run, ".png"), width = 10, height = 8)

################### Plot predicted scores vs. true scores ######################
pdf(paste0("Scripts/Plots/Model/Evaluation/pdf/True_vs_fitted_test", run, ".pdf"), width = 7, height = 3.5)
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

################################################################################
############################# Data Transformation ##############################
run = "_gam_final_trafo"

# Define the generalized logistic function
generalized_logistic <- function(x, A, k, x0, C) {
  A / (1 + exp(-k * (x - x0))) + C
}

# Inverse of the generalized logistic function
inverse_generalized_logistic <- function(y, A, k, x0, C) {
  x0 - (1 / k) * log((A / (y - C)) - 1)
}

### Estimate the parameters of the generalized logistic function on the whole data set

# Fit the nonlinear model

fit = nls(fitted(gam_final)[,1] ~ x0 - (1/k) * log((A/(d_train$PC1 - C))-1), 
          data = d_train,
          start = list(A = 12, k = 1, x0 = 0, C = -6),
          trace = TRUE)
params = coef(fit)

# Plot the fit
pdf(paste0("Scripts/Plots/Model/Evaluation/pdf/Parameters", run, "_orig.pdf"), width = 4, height = 4)
par(mfrow = c(1, 1))

plot(fitted(gam_final)[, 1], d_train$PC1, 
     pch = 16,
     xlab = "Predicted PC 1 scores",
     ylab = "True PC 1 scores",
     main = "Original Scores", 
     ylim = c(-6, 6.1), cex.main = 1.3, cex.lab = 1.1)

dev.off()

pdf(paste0("Scripts/Plots/Model/Evaluation/pdf/Parameters", run, "_params.pdf"), width = 4, height = 4)
# Construct the legend text using expression
leg.text = c(
  expression(A == 12.21),
  expression(k == 0.72),
  expression(x[0] == 0.30),
  expression(C == -6.03)
)

plot(fitted(gam_final)[, 1], d_train$PC1, 
     pch = 16,
     xlab = "Predicted PC 1 scores",
     ylab = "True PC 1 scores",
     main = "Estimated Parameter Setting", 
     ylim = c(-6, 6.1), cex.main = 1.3, cex.lab = 1.1)

curve(generalized_logistic(x, A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]), 
      from = -10, to = 10, add = TRUE, col = "red", lwd = 3)

legend(4, -2, leg.text, cex = 0.7)

dev.off()

# Save as png as well
png(paste0("Scripts/Plots/Model/Evaluation/png/Parameters", run, ".png"), width = 1000, height = 500)
leg.text = c(
  expression(A == 12),
  expression(k == 1 ),
  expression(x[0] == 0 ),
  expression(C == -6 )
)
par(mfrow = c(1, 2))

plot(fitted(gam_final)[, 1], d_train$PC1, 
     pch = 16,
     xlab = "Predicted PC 1 scores",
     ylab = "True PC 1 scores",
     main = "Manual Parameter Setting", 
     ylim = c(-6, 6.1))

curve(generalized_logistic(x, A = 11.8, k = 0.8, x0 = 0, C = -6), 
      from = -10, to = 10, add = TRUE, col = "red", lwd = 3)
legend(5, -2, leg.text, cex = 0.8)

# Construct the legend text using expression
leg.text = c(
  expression(A == 12.21),
  expression(k == 0.72),
  expression(x[0] == 0.30),
  expression(C == -6.03)
)

plot(fitted(gam_final)[, 1], d_train$PC1, 
     pch = 16,
     xlab = "Predicted PC 1 scores",
     ylab = "True PC 1 scores",
     main = "Estimated Parameter Setting", 
     ylim = c(-6, 6.1))

curve(generalized_logistic(x,  A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]), 
      from = -10, to = 10, add = TRUE, col = "red", lwd = 3)

legend(5, -2, leg.text, cex = 0.8)
dev.off()

### Use this data transformation on the model:
d_train_trafo = d_train %>%
  mutate(PC1_trafo = inverse_generalized_logistic(PC1, A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]))

gam_trafo = gam(list(PC1_trafo ~ s(Lon,Lat) +
                                   sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                                   Scenario + time_since_dist +
                                   initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                                   PC1_climate, 
                      PC2 ~ s(Lon,Lat) +
                                   sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                                   Scenario + time_since_dist +
                                   initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                                   PC1_climate), 
                  data = d_train_trafo,
                  family = mvn(d = 2))

summary(gam_trafo)
saveRDS(gam_trafo, paste0("Scripts/MA_FDA_veg/04_Model/ModelRDS/gam_final_trafo.rds"))

################################## Evaluation ##################################

### Plot true vs. fitted for training data set
pdf(paste0("Scripts/Plots/Model/Evaluation/pdf/True_vs_fitted_train", run, ".pdf"), width = 8, height = 4)
par(mfrow = c(1, 2)) 
plot(generalized_logistic(fitted(gam_trafo)[,1], A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]),
     d_train$PC1, 
     pch = 16,
     xlab = "Predicted PC 1 scores",
     ylab = "True PC 1 scores",
     cex = 1.1)
abline(a = 0, b = 1, col = "red", lwd = 2)

plot(fitted(gam_trafo)[,2], d_train$PC2, 
     pch = 16,
     xlab = "Predicted PC 2 scores",
     ylab = "True PC 2 scores",
     cex = 1.1)
abline(a = 0, b = 1, col = "red", lwd = 2)
mtext("True vs. fitted PC scores for training data", side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.3)
dev.off()

# Plot residuals vs. fitted for training data set
png(paste0("Scripts/Plots/Model/Evaluation/png/True_vs_fitted_train", run, ".png"), width = 1000, height = 500)
par(mfrow = c(1, 2)) 
plot(generalized_logistic(fitted(gam_trafo)[,1], A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]),
     d_train$PC1, 
     pch = 16,
     xlab = "Predicted PC 1 scores",
     ylab = "True PC 1 scores")
abline(a = 0, b = 1, col = "red", lwd = 2)

plot(fitted(gam_trafo)[,2], d_train$PC2, 
     pch = 16,
     xlab = "Predicted PC 2 scores",
     ylab = "True PC 2 scores")
abline(a = 0, b = 1, col = "red", lwd = 2)
mtext("True vs. fitted PC scores for training data", side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
dev.off()

pdf(paste0("Scripts/Plots/Model/Evaluation/pdf/residuals", run, ".pdf"), width = 8, height = 8)
par(mfrow = c(2, 2)) 
plot(generalized_logistic(fitted(gam_trafo)[,1], A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]),
     residuals(gam_trafo)[, 1], main = "Residuals vs. Fitted PC 1 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16, cex.main = 1.5, cex.lab = 1.3)
abline(h = 0, col = "red", lwd = 2)

qqnorm(generalized_logistic(fitted(gam_trafo)[,1], A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]), main = "Normal Q-Q (PC 1)", pch = 16, cex.main = 1.5, cex.lab = 1.3)
qqline(generalized_logistic(fitted(gam_trafo)[,1], A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]), col = "red", lwd = 2)

plot(fitted(gam_trafo)[, 2], residuals(gam_trafo)[, 2], main = "Residuals vs. Fitted PC 2 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16, cex.main = 1.5, cex.lab = 1.3)
abline(h = 0, col = "red", lwd = 2)

qqnorm(residuals(gam_trafo)[, 2], main = "Normal Q-Q (PC 2)", pch = 16, cex.main = 1.5, cex.lab = 1.3)
qqline(residuals(gam_trafo)[, 2], col = "red", lwd = 2)
dev.off()

# Save as png as well
png(paste0("Scripts/Plots/Model/Evaluation/png/residuals", run, ".png"), width = 1000, height = 1000)
par(mfrow = c(2, 2)) 
plot(generalized_logistic(fitted(gam_trafo)[,1], A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]),
     residuals(gam_trafo)[, 1], main = "Residuals vs. Fitted PC 1 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(generalized_logistic(fitted(gam_trafo)[,1], A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]), main = "Normal Q-Q (PC 1)", pch = 16)
qqline(generalized_logistic(fitted(gam_trafo)[,1], A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]), col = "red")

plot(fitted(gam_trafo)[, 2], residuals(gam_trafo)[, 2], main = "Residuals vs. Fitted PC 2 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

qqnorm(residuals(gam_trafo)[, 2], main = "Normal Q-Q (PC 2)", pch = 16)
qqline(residuals(gam_trafo)[, 2], col = "red")
dev.off()

## Predict on the test set
d_test_trafo = d_test %>%
  mutate(PC1_trafo = inverse_generalized_logistic(PC1, A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]))

predictions = predict(gam_trafo, newdata = d_test_trafo)

MFPCA_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_all_1803_100y.rds")
funData_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_1803.rds")

uniExpansions <- list(list(type = "uFPCA", npc = 10),
                      list(type = "uFPCA", npc = 10),
                      list(type = "uFPCA", npc = 10),
                      list(type = "uFPCA", npc = 10),
                      list(type = "uFPCA", npc = 10))


MFPCA_2PCs <- MFPCA(funData_all, M = 2, fit = TRUE, uniExpansions = uniExpansions)
#saveRDS(MFPCA_2PCs, paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_all_1803_100y_2PCs", run, ".rds"))

pred_MFPCA = predict(MFPCA_2PCs, 
                     scores = cbind(generalized_logistic(predictions[,1], A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]), predictions[,2]))

MFPCA_test = MFPCA_2PCs[-index]
funData_test = funData_all[-index]

########################### Plot example functional fits #######################
# Create a data frame from the provided data

j=0
dat = c()
for (iCurve in c(2,89,167,259)){
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
    x = "Year after Disturbance",
    y = "Share of aboveground carbon",
    title = "True functional fit and predicted one for example grid cells",
    color = "Dominant vegetation",
    linetype = "Model"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.box = "vertical"
  ) +
  scale_color_manual(values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                                "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73"))

ggsave(paste0("Scripts/Plots/Model/Evaluation/pdf/True_vs_fitted_examples", run, ".pdf"), width = 10, height = 6)
ggsave(paste0("Scripts/Plots/Model/Evaluation/png/True_vs_fitted_examples", run, ".png"), width = 10, height = 6)

################### Plot predicted scores vs. true scores ######################
pdf(paste0("Scripts/Plots/Model/Evaluation/pdf/True_vs_fitted_test", run, ".pdf"), width = 8, height = 4)
par(mfrow=c(1,2))
plot(generalized_logistic(predictions[,1], A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]), MFPCA_2PCs$scores[-index,1],
     pch = 16,
     xlab = "Predicted PC 1 scores",
     ylab = "True PC 1 scores",
     cex = 1.1)
abline(a = 0, b = 1, col = "red", lwd = 2)
plot(predictions[,2], MFPCA_2PCs$scores[-index,2], 
     pch = 16,
     xlab = "Predicted PC 2 scores",
     ylab = "True PC 2 scores",
     cex = 1.1)
abline(a = 0, b = 1, col = "red", lwd = 2)
mtext("True vs. fitted PC scores for unseen data", side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.3)
dev.off()

# Save as png as well
png(paste0("Scripts/Plots/Model/Evaluation/png/True_vs_fitted_test", run, ".png"), width = 1000, height = 500)
par(mfrow=c(1,2))
plot(generalized_logistic(predictions[,1], A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]), MFPCA_2PCs$scores[-index,1],
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

################################ Correlations ##################################
library(ggcorrplot)
ggcorrplot(cor(combined_d %>% select(-Scenario)), 
           type = "lower", 
           method = "square", 
           lab = TRUE,
           lab_size = 1.9,  
           ggtheme = theme_bw() + theme(text = element_text(size = 10), 
                                        plot.title = element_text(size = 10, face = "bold", hjust = 0.5)),
           legend.title = "Correlation")


ggsave(paste0("Scripts/Plots/Model/Evaluation/pdf/Correlations.pdf"), width = 10, height = 10)
ggsave(paste0("Scripts/Plots/Model/Evaluation/png/Correlations.png"), width = 15, height = 12)


################################## Plot T values #############################
summary_gam_trafo = summary(gam_trafo)
t_values_PC1 = data.frame(var = attr(summary_gam_trafo$p.t, "names")[1:16], p.t = summary_gam_trafo$p.t[1:16])
t_values_PC1$var <- factor(t_values_PC1$var, levels = t_values_PC1$var)

## PC 1
ggplot(t_values_PC1[c(1:16),], aes(x = p.t, y = var)) +
  geom_bar(stat = "identity", aes(fill = p.t > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(x = "t-value", y = "Predictors") + 
  ggtitle("t-values of PC 1") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black")  # Add dotted vertical lines

ggsave(paste0("Scripts/Plots/Model/Evaluation/pdf/t-values_PC1.pdf"), width = 8, height = 5)
ggsave(paste0("Scripts/Plots/Model/Evaluation/png/t-values_PC1.png"), width = 10, height = 8)

## PC 2

t_values_PC2 = data.frame(var = attr(summary_gam_trafo$p.t, "names")[1:16], p.t = summary_gam_trafo$p.t[17:32])
t_values_PC2$var <- factor(t_values_PC2$var, levels = t_values_PC2$var)

ggplot(t_values_PC2[c(1:16),], aes(x = p.t, y = var)) +
  geom_bar(stat = "identity", aes(fill = p.t > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(x = "t-value", y = "Predictors") + 
  ggtitle("t-values of PC 2") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black")  # Add dotted vertical lines

ggsave(paste0("Scripts/Plots/Model/Evaluation/pdf/t-values_PC2.pdf"), width = 8, height = 5)
ggsave(paste0("Scripts/Plots/Model/Evaluation/png/t-values_PC2.png"), width = 10, height = 8)

######################## Plot smooth effect manually ###########################

df_smooth_PC1 = data.frame(Lon = d_train_trafo$Lon, Lat = d_train_trafo$Lat, values = predict(gam_trafo, type = "terms")[,27])
df_smooth_PC2 = data.frame(Lon = d_train_trafo$Lon, Lat = d_train_trafo$Lat, values = predict(gam_trafo, type = "terms")[,28])

## As maps
library(ggmap)
register_google(key = "AIzaSyA_eUpOhj7hoPyzyynWvyMqcGEA1Z_SZVY")
worldmap <- get_map(location = c(lon =-4.068561, lat = 58.87355), zoom = 1)

ggmap(worldmap) +
  geom_point(data = df_smooth_PC1, aes(x = Lon, y = Lat, color = values), size = 1) +
  #ggtitle(paste(long_names_scenarios(iScen), " - ", iYear)) +
  #geom_ribbon(data = NULL, aes(ymin = -Inf, ymax = 60), fill = "grey", alpha = 0.05) +
  xlab("Longitude") + ylim(40,75) + xlim(-200,200) +
  ylab("Latitude") +
  scale_color_gradient2(name = "Values", low = "cornflowerblue", mid = "white", high = "darkred", midpoint = 0, 
                        limits = c(min(df_smooth_PC1$values), max(df_smooth_PC1$values)),
                        breaks = c(-1, 0, 1),
                        labels = c("-", "", "+")) +
  theme_bw() + ggtitle("Bivariate additive effect for transformed PC 1 scores") +
  theme(text = element_text(size = 15), plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

ggsave(paste0("Scripts/Plots/Model/Evaluation/pdf/map_smooth_PC1.pdf"), width = 10, height = 8)

ggmap(worldmap) +
  geom_point(data = df_smooth_PC2, aes(x = Lon, y = Lat, color = values), size = 1) +
  #ggtitle(paste(long_names_scenarios(iScen), " - ", iYear)) +
  #geom_ribbon(data = NULL, aes(ymin = -Inf, ymax = 60), fill = "grey", alpha = 0.05) +
  xlab("Longitude") + ylim(40,75) + xlim(-200,200) +
  ylab("Latitude") +
  scale_color_gradient2(name = "Values", low = "cornflowerblue", mid = "white", high = "darkred", midpoint = 0, 
                        limits = c(min(df_smooth_PC2$values), max(df_smooth_PC2$values)),
                        breaks = c(-1, 0, 1),
                        labels = c("-", "", "+")) +
  theme_bw() + ggtitle("Bivariate additive effect for PC 2 scores") +
  theme(text = element_text(size = 15), plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

ggsave(paste0("Scripts/Plots/Model/Evaluation/pdf/map_smooth_PC2.pdf"), width = 10, height = 8)

