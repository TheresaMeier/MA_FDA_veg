################################################################################
############################ Master's Thesis ###################################
################################################################################

########################## Model: Separate Models ##############################

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

# Get target variable: MFPCA scores
plot_data_all = read.table("Scripts/Plots/MFPCA/plot_data/plot_data_all.txt")

# Get soil and ecological variables
d_soil = read.table("Scripts/MA_FDA_veg/04_Model/Data/d_soil.txt")
d_eco_wide = read.table("Scripts/MA_FDA_veg/04_Model/Data/d_eco_wide.txt")

# Get FPCA results for temp, precip and nuptake_total
scores_climate = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_climate.txt")

# Combine data set
combined_d = cbind(plot_data_all[,c(1:2)], d_soil, d_eco_wide[,c(4:18)], 
                   scores_climate[,1])  %>% 
  mutate(Scenario = as.factor(Scenario))

colnames(combined_d)[c(1:2,27)] = c("PC1", "PC2", "PC1_climate")


################################### Model ######################################

# Produce train-/test- split
set.seed(1)
index = sample.int(nrow(combined_d), size = ceiling(0.8 * nrow(combined_d)))
d_train = combined_d[index,]
d_test = combined_d[-index,]


# Define the generalized logistic function
generalized_logistic <- function(x, A, k, x0, C) {
  A / (1 + exp(-k * (x - x0))) + C
}

# Inverse of the generalized logistic function
inverse_generalized_logistic <- function(y, A, k, x0, C) {
  x0 - (1 / k) * log((A / (y - C)) - 1)
}

### Estimate the parameters of the generalized logistic function on the whole data set
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

# Fit the nonlinear model

fit = nls(fitted(gam_final)[,1] ~ x0 - (1/k) * log((A/(d_train$PC1 - C))-1), 
          data = d_train,
          start = list(A = 12, k = 1, x0 = 0, C = -6),
          trace = TRUE)
params = coef(fit)

d_train_trafo = d_train %>%
  mutate(PC1_trafo = inverse_generalized_logistic(PC1, A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]))

gam_trafo = gam(list(PC1_trafo ~ s(Lon,Lat) +
                       sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                       #Scenario + 
                       time_since_dist +
                       initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                       PC1_climate, 
                     PC2 ~ s(Lon,Lat) +
                       sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                       #Scenario +
                       time_since_dist +
                       initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                       PC1_climate), 
                data = d_train_trafo,
                family = mvn(d = 2))

sum_combined = summary(gam_trafo)
# saveRDS(gam_trafo, paste0("Scripts/MA_FDA_veg/04_Model/ModelRDS/gam_final_trafo_all.rds"))

##################### Derive separate models for each scenario #################

# Control
gam_picontrol = gam(list(PC1_trafo ~ s(Lon,Lat) +
                           sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                           #Scenario + 
                           time_since_dist +
                           initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                           PC1_climate, 
                          PC2 ~ s(Lon,Lat) +
                           sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                           #Scenario + 
                           time_since_dist +
                           initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                           PC1_climate), 
                  data = d_train_trafo[d_train_trafo$Scenario == "Control",],
                  family = mvn(d = 2))

sum_picontrol = summary(gam_picontrol)

# SSP1-RCP2.6
gam_ssp126 = gam(list(PC1_trafo ~ s(Lon,Lat) +
                        sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                        #Scenario + 
                        time_since_dist +
                        initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                        PC1_climate, 
                      PC2 ~ s(Lon,Lat) +
                        sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                        #Scenario + 
                        time_since_dist +
                        initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                        PC1_climate), 
                 data = d_train_trafo[d_train_trafo$Scenario == "SSP1-RCP2.6",],
                 family = mvn(d = 2))

sum_ssp126 = summary(gam_ssp126)

# SSP3-RCP7.0
gam_ssp370 = gam(list(PC1_trafo ~ s(Lon,Lat) +
                        sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                        #Scenario + 
                        time_since_dist +
                        initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                        PC1_climate, 
                      PC2 ~ s(Lon,Lat) +
                        sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                        #Scenario + 
                        time_since_dist +
                        initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                        PC1_climate), 
                 data = d_train_trafo[d_train_trafo$Scenario == "SSP3-RCP7.0",],
                 family = mvn(d = 2))

sum_ssp370 = summary(gam_ssp370)

# SSP5-RCP8.5
gam_ssp585 = gam(list(PC1_trafo ~ s(Lon,Lat) +
                        sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                        #Scenario + 
                        time_since_dist +
                        initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                        PC1_climate, 
                      PC2 ~ s(Lon,Lat) +
                        sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                        #Scenario + 
                        time_since_dist +
                        initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                        PC1_climate), 
                 data = d_train_trafo[d_train_trafo$Scenario == "SSP5-RCP8.5",],
                 family = mvn(d = 2))

sum_ssp585 = summary(gam_ssp585)

########################## Compare results #####################################

# Load required libraries
library(ggplot2)
library(reshape2)
library(patchwork)  # For arranging the plots

### t-values
## PC 1
t_scen_PC1 = as.data.frame(cbind(attr(sum_combined$p.t, "names"), sum_combined$p.t, sum_picontrol$p.t, sum_ssp126$p.t,
                                 sum_ssp370$p.t, sum_ssp585$p.t)[1:13,])
colnames(t_scen_PC1) = c("Variable", "Combined", "Control", "SSP1-RCP2.6", "SSP3-RCP7.0", "SSP5-RCP8.5")
t_scen_PC1$Variable = factor(t_scen_PC1$Variable, levels = t_scen_PC1$Variable)
t_scen_PC1[,2:6] = as.numeric(unlist(t_scen_PC1[,2:6]))

## PC 2
t_scen_PC2 = as.data.frame(cbind(attr(sum_combined$p.t, "names")[1:13], cbind(sum_combined$p.t, sum_picontrol$p.t, sum_ssp126$p.t,
                                 sum_ssp370$p.t, sum_ssp585$p.t)[14:26,]))
colnames(t_scen_PC2) = c("Variable", "Combined", "Control", "SSP1-RCP2.6", "SSP3-RCP7.0", "SSP5-RCP8.5")
t_scen_PC2$Variable = factor(t_scen_PC2$Variable, levels = t_scen_PC2$Variable)
t_scen_PC2[,2:6] = as.numeric(unlist(t_scen_PC2[,2:6]))

## Plots for PC 1
plot_combined = ggplot(t_scen_PC1[c(1:13),], aes(x = Combined, y = Variable)) +
  geom_bar(stat = "identity", aes(fill = Combined > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(x = "t-value", y = "Predictors", title = "All scenarios combined")  + xlim(-5.5,16) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = "black")  # Add dotted vertical lines

ggsave(paste0("Scripts/Plots/Model/SeparateModels/pdf/t-values_PC1_combined.pdf"), plot_combined, width = 6, height = 3)
ggsave(paste0("Scripts/Plots/Model/SeparateModels/png/t-values_PC1_combined.png"), plot_combined, width = 6, height = 3)

plot_picontrol = ggplot(t_scen_PC1[c(1:13),], aes(x = Control, y = Variable)) +
  geom_bar(stat = "identity", aes(fill = Control > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(x = "t-value", y = "Predictors", title = "Control") + xlim(-5.5,16) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = "black")  # Add dotted vertical lines

ggsave(paste0("Scripts/Plots/Model/SeparateModels/pdf/t-values_PC1_control.pdf"), plot_picontrol, width = 6, height = 3)
ggsave(paste0("Scripts/Plots/Model/SeparateModels/png/t-values_PC1_control.png"), plot_picontrol, width = 6, height = 3)

plot_ssp126 = ggplot(t_scen_PC1[c(1:13),], aes(x = `SSP1-RCP2.6`, y = Variable)) +
  geom_bar(stat = "identity", aes(fill = `SSP1-RCP2.6` > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(x = "t-value", y = "Predictors", title = "SSP1-RCP2.6") + xlim(-5.5,16) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = "black")  # Add dotted vertical lines

ggsave(paste0("Scripts/Plots/Model/SeparateModels/pdf/t-values_PC1_ssp126.pdf"), plot_ssp126, width = 6, height = 3)
ggsave(paste0("Scripts/Plots/Model/SeparateModels/png/t-values_PC1_ssp126.png"), plot_ssp126, width = 6, height = 3)

plot_ssp370 = ggplot(t_scen_PC1[c(1:13),], aes(x = `SSP3-RCP7.0`, y = Variable)) +
  geom_bar(stat = "identity", aes(fill = `SSP3-RCP7.0` > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(x = "t-value", y = "Predictors", title = "SSP3-RCP7.0") + xlim(-5.5,16) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = "black")  # Add dotted vertical lines

ggsave(paste0("Scripts/Plots/Model/SeparateModels/pdf/t-values_PC1_ssp370.pdf"), plot_ssp370, width = 6, height = 3)
ggsave(paste0("Scripts/Plots/Model/SeparateModels/png/t-values_PC1_ssp370.png"), plot_ssp370, width = 6, height = 3)


plot_ssp585 = ggplot(t_scen_PC1[c(1:13),], aes(x = `SSP5-RCP8.5`, y = Variable)) +
  geom_bar(stat = "identity", aes(fill = `SSP5-RCP8.5` > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(x = "t-value", y = "Predictors", title = "SSP5-RCP8.5") + xlim(-5.5,16) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = "black")  # Add dotted vertical lines

ggsave(paste0("Scripts/Plots/Model/SeparateModels/pdf/t-values_PC1_ssp585.pdf"), plot_ssp585, width = 6, height = 3)
ggsave(paste0("Scripts/Plots/Model/SeparateModels/png/t-values_PC1_ssp585.png"), width = 6, height = 3)

# Arrange the plots
layout <- (
  (plot_spacer() / plot_combined / plot_spacer()) | # Center the Combined plot using spacers
    ((plot_picontrol/ plot_ssp370) | (plot_ssp126 / plot_ssp585))
) + 
  plot_layout(widths = c(1,2), heights = c(1, 1, 1))  # Adjust the layout widths and heights

# Show the layout
layout

ggsave(paste0("Scripts/Plots/Model/SeparateModels/pdf/t-values_PC1_all.pdf"), layout, width = 18, height = 6)
ggsave(paste0("Scripts/Plots/Model/SeparateModels/png/t-values_PC1_all.png"), layout, width = 18, height = 6)

## Plots for PC 2
plot_combined = ggplot(t_scen_PC2[c(1:13),], aes(x = Combined, y = Variable)) +
  geom_bar(stat = "identity", aes(fill = Combined > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(x = "t-value", y = "Predictors", title = "All scenarios combined")  + xlim(-11, 5) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = "black")  # Add dotted vertical lines

ggsave(paste0("Scripts/Plots/Model/SeparateModels/pdf/t-values_PC2_combined.pdf"), plot_combined, width = 6, height = 3)
ggsave(paste0("Scripts/Plots/Model/SeparateModels/png/t-values_PC2_combined.png"), plot_combined, width = 6, height = 3)

plot_picontrol = ggplot(t_scen_PC2[c(1:13),], aes(x = Control, y = Variable)) +
  geom_bar(stat = "identity", aes(fill = Control > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(x = "t-value", y = "Predictors", title = "Control") + xlim(-11, 5) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = "black")  # Add dotted vertical lines

ggsave(paste0("Scripts/Plots/Model/SeparateModels/pdf/t-values_PC2_control.pdf"), plot_picontrol, width = 6, height = 3)
ggsave(paste0("Scripts/Plots/Model/SeparateModels/png/t-values_PC2_control.png"), plot_picontrol, width = 6, height = 3)

plot_ssp126 = ggplot(t_scen_PC2[c(1:13),], aes(x = `SSP1-RCP2.6`, y = Variable)) +
  geom_bar(stat = "identity", aes(fill = `SSP1-RCP2.6` > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(x = "t-value", y = "Predictors", title = "SSP1-RCP2.6") + xlim(-11, 5) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = "black")  # Add dotted vertical lines

ggsave(paste0("Scripts/Plots/Model/SeparateModels/pdf/t-values_PC2_ssp126.pdf"), plot_ssp126, width = 6, height = 3)
ggsave(paste0("Scripts/Plots/Model/SeparateModels/png/t-values_PC2_ssp126.png"), plot_ssp126, width = 6, height = 3)

plot_ssp370 = ggplot(t_scen_PC2[c(1:13),], aes(x = `SSP3-RCP7.0`, y = Variable)) +
  geom_bar(stat = "identity", aes(fill = `SSP3-RCP7.0` > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(x = "t-value", y = "Predictors", title = "SSP3-RCP7.0") + xlim(-11, 5) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = "black")  # Add dotted vertical lines

ggsave(paste0("Scripts/Plots/Model/SeparateModels/pdf/t-values_PC2_ssp370.pdf"), plot_ssp370, width = 6, height = 3)
ggsave(paste0("Scripts/Plots/Model/SeparateModels/png/t-values_PC2_ssp370.png"), plot_ssp370, width = 6, height = 3)


plot_ssp585 = ggplot(t_scen_PC2[c(1:13),], aes(x = `SSP5-RCP8.5`, y = Variable)) +
  geom_bar(stat = "identity", aes(fill = `SSP5-RCP8.5` > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(x = "t-value", y = "Predictors", title = "SSP5-RCP8.5") + xlim(-11, 5) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = "black")  # Add dotted vertical lines

ggsave(paste0("Scripts/Plots/Model/SeparateModels/pdf/t-values_PC2_ssp585.pdf"), plot_ssp585, width = 6, height = 3)
ggsave(paste0("Scripts/Plots/Model/SeparateModels/png/t-values_PC2_ssp585.png"), width = 6, height = 3)

# Arrange the plots
layout <- (
  (plot_spacer() / plot_combined / plot_spacer()) | # Center the Combined plot using spacers
    ((plot_picontrol/ plot_ssp370) | (plot_ssp126 / plot_ssp585))
) + 
  plot_layout(widths = c(1,2), heights = c(1, 1, 1))  # Adjust the layout widths and heights

# Show the layout
layout

ggsave(paste0("Scripts/Plots/Model/SeparateModels/pdf/t-values_PC2_all.pdf"), layout, width = 18, height = 6)
ggsave(paste0("Scripts/Plots/Model/SeparateModels/png/t-values_PC2_all.png"), layout, width = 18, height = 6)

