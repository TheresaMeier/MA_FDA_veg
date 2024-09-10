################################################################################
############################ Master's Thesis ###################################
################################################################################

############################ Model: Scenario wise###############################

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

##################### Run model for each scenario and PID ######################
row_nr = 3
df_beta = array(0, dim = c((16 - row_nr) * 2, 25, 4))
df_t = array(0, dim = c((16 - row_nr) * 2, 25, 4))

for (pid in 0:24){
  print(paste0("Start with PID ", pid))
  
  ############################ Get the prepared data ###########################
  # Get target variable: MFPCA scores
  MFPCA_all = readRDS(paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/AllPatches/MFPCA_PID", pid, ".rds"))
  plot_data_all = MFPCA_all$scores
  #plot_data_all = read.table("Scripts/Plots/MFPCA/plot_data/plot_data_all.txt")
  
  # Get soil and ecological variables
  d_soil = read.table(paste0("Scripts/MA_FDA_veg/04_Model/Data/AllPatches/d_soil_PID", pid, ".txt"))
  d_eco_wide = read.table(paste0("Scripts/MA_FDA_veg/04_Model/Data/AllPatches/d_eco_wide_PID", pid, ".txt"))
  
  # Get MFPCA results for temp (min, max, mean) and precip
  scores_climate = read.table(paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/AllPatches/scores_climate_PID", pid, ".txt"))
  scores_climate_picontrol = read.table(paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/AllPatches/scores_climate_picontrol_PID", pid, ".txt"))
  scores_climate_ssp126 = read.table(paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/AllPatches/scores_climate_ssp126_PID", pid, ".txt"))
  scores_climate_ssp370 = read.table(paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/AllPatches/scores_climate_ssp370_PID", pid, ".txt"))
  scores_climate_ssp585 = read.table(paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/AllPatches/scores_climate_ssp585_PID", pid, ".txt"))
  
  # Combine data set
  counts = table(d_soil$Scenario)
  starts <- c(1, cumsum(counts) + 1)  
  ends <- cumsum(counts)  
  # Combine starts and ends into a matrix
  bounds <- cbind(starts[-length(starts)], ends)
  
  i=0
  for (scen in scenarios){
    i=i+1
    scores_climate_scen = get(paste0("scores_climate_", scen))
    combined_d_scen = cbind(plot_data_all[bounds[i,1]:bounds[i,2], c(1:2)], d_soil[bounds[i,1]:bounds[i,2],], d_eco_wide[bounds[i,1]:bounds[i,2],c(4:18)])  %>% 
    mutate(Scenario = as.factor(Scenario))
  
    colnames(combined_d_scen)[c(1:2,27)] = c("PC1", "PC2", "PC1_climate")
    assign(paste0("combined_d_", scen), combined_d_scen)
  }
  
  ################################### Model ######################################
  # Define parameters for glf transformation (taken from PID 1)
  params <- matrix(c(
    17.65680, 0.3900519,  0.26761163, -9.509338,
    12.09339, 0.6703681, -0.19841980, -6.191885,
    12.11618, 0.6728557, -0.22644067, -6.122977,
    13.08880, 0.6087930,  0.09764066, -6.259491
  ), nrow = 4, byrow = TRUE)
  
  rownames(params) <- c("params_picontrol", "params_ssp126", "params_ssp370", "params_ssp585")
  colnames(params) <- c("A", "k", "x0", "C")  
  
  i=0
  for (scen in scenarios){
    combined_d = get(paste0("combined_d_", scen))
    i=i+1
    # Produce train-/test- split
    set.seed(1)
    index = sample.int(nrow(combined_d), size = ceiling(0.8 * nrow(combined_d)))
    d_train = combined_d[index,]
    d_test = combined_d[-index,]
    
    # Fit the model using non-linear effects
    gam_final = gam(list(PC1 ~ s(Lon,Lat) +
                           sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                           time_since_dist +
                           initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                           PC1_climate, 
                         PC2 ~ s(Lon,Lat) +
                           sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                           time_since_dist +
                           initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                           PC1_climate), 
                    data = d_train,
                    family = mvn(d = 2))
    
    #summary(gam_final)
    
    # Define the generalized logistic function
    generalized_logistic <- function(x, A, k, x0, C) {A / (1 + exp(-k * (x - x0))) + C}
    
    # Inverse of the generalized logistic function
    inverse_generalized_logistic <- function(y, A, k, x0, C) {x0 - (1 / k) * log((A / (y - C)) - 1)}
    
  
    ### Use this data transformation on the model:
    d_train_trafo = d_train %>%
      mutate(PC1_trafo = inverse_generalized_logistic(PC1, A = params[paste0("params_",scen), "A"], 
                                                      k = params[paste0("params_",scen), "k"], 
                                                      x0 = params[paste0("params_",scen), "x0"],
                                                      C = params[paste0("params_",scen), "C"]))
    
    gam_trafo = gam(list(PC1_trafo ~ s(Lon,Lat) +
                           sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                           time_since_dist +
                           initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                           PC1_climate, 
                         PC2 ~ s(Lon,Lat) +
                           sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                           time_since_dist +
                           initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                           PC1_climate), 
                    data = d_train_trafo,
                    family = mvn(d = 2))

    df_beta[,pid+1,i] = gam_trafo$coefficients[c(1:(16-row_nr),(46-row_nr):(61 - 2*row_nr))]
    summary_gam_trafo = summary(gam_trafo)
    df_t[,pid+1,i] = summary_gam_trafo$p.t
  }
}

dimnames(df_beta) = list(names(gam_trafo$coefficients[c(1:(16-row_nr),(46-row_nr):(61 - 2*row_nr))]),
                         NULL, NULL)

################################### Plot t-values ##############################
## Calculate CIs
calculate_ci <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  std_error <- 1/sqrt(n) * sd(x) 
  error_margin <- 1.96 * std_error
  lower_bound <- mean_x - error_margin
  upper_bound <- mean_x + error_margin
  return(c(mean = mean_x, lower = lower_bound, upper = upper_bound))
}

for (i in 1:4){
  ci_data_scen = as.data.frame(t(apply(df_beta[,,i], 1, calculate_ci)))
  ci_data_scen = ci_data_scen %>%
    mutate(Effects = if_else(lower > 0 & upper > 0, "Positive",
                             if_else((lower < 0 & upper < 0), "Negative","Non-significant")),
           Effects = factor(Effects, levels = c("Positive", "Negative", "Non-significant")),
           parameter = factor(rownames(df_beta), levels = rownames(df_beta)),
           t = apply(df_beta[,,i], 1, function(x)  sqrt(length(x)) * mean(as.numeric(x)) / sd(as.numeric(x))))
  
  assign(paste0("ci_data_", scenarios[i]), ci_data_scen)
  
  # PC1
  g_PC1_scen = ggplot(ci_data_scen[c(1:(16-row_nr)),], aes(x = t, y = parameter)) +
    geom_bar(stat = "identity", aes(fill = t > 0)) +
    scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
    theme_bw() + theme(text = element_text(size = 15),
                       legend.position = "none",
                       plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
    labs(x = "t-Value", y = "Predictors") + 
    ggtitle("t-Values of PC 1 using all patches") +
    geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", color = "black") +
    ggtitle(paste0("t-Values of PC 1 - ", long_names_scenarios(scenarios[i]))) + xlim(-28,75)
  
  ggsave(paste0("Scripts/Plots/Model/SeparateModels/pdf/t-values_patches_PC1_", scenarios[i],  ".pdf"), g_PC1_scen, width = 7, height = 4.3)
  ggsave(paste0("Scripts/Plots/Model/SeparateModels/png/t-values_patches_PC1_", scenarios[i],  ".png"), g_PC1_scen, width = 7, height = 4.3)
  
  assign(paste0("g_PC1_", scenarios[i]), g_PC1_scen)
  
  # PC2
  ci_PC2_scen = ci_data_scen[(17-row_nr):(32-2*row_nr),]
  ci_PC2_scen$parameter <- factor(rownames(df_beta)[1:(16-row_nr)], levels = rownames(df_beta)[1:(16-row_nr)])
  
  g_PC2_scen = ggplot(ci_PC2_scen, aes(x = t, y = parameter)) +
    geom_bar(stat = "identity", aes(fill = t > 0)) +
    scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
    theme_bw() + theme(text = element_text(size = 15),
                       legend.position = "none",
                       plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
    labs(x = "t-Value", y = "Predictors") + 
    ggtitle("t-Values of PC 2 using all patches") +
    geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", color = "black") +
    ggtitle(paste0("t-Values of PC 2 - ", long_names_scenarios(scenarios[i]))) + xlim(-16,15)
  
  ggsave(paste0("Scripts/Plots/Model/SeparateModels/pdf/t-values_patches_PC2_", scenarios[i],  ".pdf"), g_PC2_scen, width = 7, height = 4.3)
  ggsave(paste0("Scripts/Plots/Model/SeparateModels/png/t-values_patches_PC2_", scenarios[i],  ".png"), g_PC2_scen, width = 7, height = 4.3)
  
  assign(paste0("g_PC2_", scenarios[i]), g_PC2_scen)
}

# Plot them altogether
plot_grid(g_PC1_picontrol, g_PC1_ssp126, g_PC1_ssp370, g_PC1_ssp585)
plot_grid(g_PC2_picontrol, g_PC2_ssp126, g_PC2_ssp370, g_PC2_ssp585)

