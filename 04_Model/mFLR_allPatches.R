################################################################################
############################ Master's Thesis ###################################
################################################################################

######################### Model: Different Patches #############################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/02_FPCA/functions.R")
source("Scripts/MA_FDA_veg/03_MFPCA/MFPCA/MFPCA_calculation.R")
source("Scripts/MA_FDA_veg/01_Description/utils.R")

## Load libraries
library(duckdb)
library(tidyverse)
library(ggplot2)
library(fda)
library(ggmap)
library(gifski)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
library(grid)
library(abind)
library(MFPCA)
library(foreach)
library(funData)
library(mgcv)
library(data.table)

## Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)
scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")
pfts = c("BNE", "IBS", "otherC", "TeBS", "Tundra")
set.seed(1)

## Choose parameters:
start_year = 2015
end_year = 2040
pid = 1           # Choose pid (int) or 'all'
k = 4             # Number of clusters for kmeans algorithm
M = 10            # Number of PCs
doCreate = FALSE
doPrepare = FALSE
whichScen = "all"   # Choose from "Control", "SSP1-RCP2.6", "SSP3-RCP7.0", "SSP5-RCP8.5" and "all"
covScen = FALSE       # Always equal to FALSE if whichScen != "all"
row_nr = 0

######################### Create funData and run MFPCA ########################
if (doCreate){
  for (pid in c(0:24)){
    print(paste("Start with PID", pid))
    # Create funData object 
    for (pft in pfts){
      print(paste("Start with PFT", long_names_pfts(tolower(pft))))
      
      for (scen in scenarios){
        print(long_names_scenarios(scen))
        d_scen = get_data_fpca(scen, start_year, end_year,pid,pft)
        assign(paste0("d_",scen), d_scen)
        print("...done.")
      }
      d_pft = rbind(t(d_picontrol[[1]][,-1]), t(d_ssp126[[1]][,-1]), t(d_ssp370[[1]][,-1]), t(d_ssp585[[1]][,-1]))
      assign(paste0("d_", pft, "_", pid), d_pft)
      
      locs_disturbed_pid = rbind(d_picontrol[[2]], d_ssp126[[2]], d_ssp370[[2]], d_ssp585[[2]])
      locs_disturbed_pid$scenario = c(rep("Control", length(d_picontrol[[1]])-1), rep("SSP1-RCP2.6", length(d_ssp126[[1]])-1), 
                                      rep("SSP3-RCP7.0", length(d_ssp370[[1]])-1), rep("SSP5-RCP8.5",length(d_ssp585[[1]])-1))
      
      write.table(locs_disturbed_pid, paste0("Scripts/Plots/MFPCA/plot_data/AllPatches/locs_disturbed_PID", pid, ".txt"))
    }

    # Create univariate funData objects
    for (pft in pfts) {
      fun_pft <- funData(argvals = 1:100,
                         X = get(paste0("d_", pft, "_", pid)))
      dimnames(fun_pft@X)[[2]] = 1:100
      assign(paste0("fun_", pft), fun_pft)
      saveRDS(fun_pft, paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/AllPatches/funData_", pft, "_PID", pid, ".rds"))
    }

    # Create multiFunData object
    multiFun_pft = multiFunData(fun_BNE, fun_IBS, fun_otherC, fun_TeBS, fun_Tundra)

    saveRDS(multiFun_pft, paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/AllPatches/funData_PID", pid, ".rds"))

    # Run MFPCA
    uniExpansions <- list(list(type = "uFPCA", npc = 10),
                          list(type = "uFPCA", npc = 10),
                          list(type = "uFPCA", npc = 10),
                          list(type = "uFPCA", npc = 10),
                          list(type = "uFPCA", npc = 10))

    MFPCA_all <- MFPCA_2(multiFun_pft, M = M, fit = TRUE, uniExpansions = uniExpansions)
    saveRDS(MFPCA_all, paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/AllPatches/AllPatches/MFPCA_PID", pid, ".rds"))
  }
}


############################### Examine PCs ####################################

# MFPCA_0 = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/AllPatches/MFPCA_PID24.rds")
# plot(MFPCA_0, plotPC = 1:2, combined = TRUE)

########################### Prepare data for modelling #########################
if (doPrepare){
  
  ##### Ecological and soil covariates
  for (scen in scenarios){
    d_scen_2015_2040 = fread(paste0("Data/data_", scen, "_2015_2040.csv"))
    assign(paste0("d_", scen, "_2015_2040"), d_scen_2015_2040)
  }
  
  for (pid in c(0:24)){
    print(paste0("Start with PID ", pid))
    # Get MFPCA scores
    MFPCA_pid = readRDS(paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/AllPatches/MFPCA_PID", pid, ".rds"))
    scores_pid = as.data.frame(MFPCA_pid$scores)
    colnames(scores_pid) = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
    
    # Get disturbed loctions
    locs_disturbed_pid = read.table(paste0("Scripts/Plots/MFPCA/plot_data/AllPatches/locs_disturbed_PID", pid, ".txt"))
    
    # Get other time-constant soil covariates
    for (scen in scenarios){
      d_scen_2015_2040 = get(paste0("d_", scen, "_2015_2040"))
      
      d_soil_scen = d_scen_2015_2040 %>%
        filter(PID == pid) %>%
        distinct(Lon, Lat, sand_fraction, silt_fraction, clay_fraction, bulkdensity_soil, ph_soil, soilcarbon, time_since_dist) %>%
        mutate(Scenario = long_names_scenarios(scen))
      
      d_eco_scen = d_scen_2015_2040 %>%
        filter(PID == pid) %>%
        distinct(Lon,Lat,PFT,Year,initial_recruitment, recruitment_ten_years, previous_state, age, Nuptake) %>%
        mutate(Scenario = long_names_scenarios(scen))
      
      d_climate_scen = d_scen_2015_2040 %>%
        filter(PID == pid) %>%
        distinct(Lon,Lat,Year, age, tas_yearlymeam, tas_yearlymin, tas_yearlymax, pr_yearlysum) %>%
        mutate(Scenario = long_names_scenarios(scen))
      
      assign(paste0("d_climate_", scen), d_climate_scen)
      assign(paste0("d_soil_", scen), d_soil_scen)
      assign(paste0("d_eco_", scen), d_eco_scen)
    }
    
    ### Soil covariates
    for (var in c("soil", "eco")){
      
      d_var = rbind(get(paste0("d_", var, "_picontrol")),
                    get(paste0("d_", var, "_ssp126")), 
                    get(paste0("d_", var, "_ssp370")), 
                    get(paste0("d_", var, "_ssp585")))
      
      if (var == "soil"){
        # Reorder d_var to match target
        d_var$key = paste(d_var$Lon, d_var$Lat, d_var$Scenario)
        locs_disturbed_pid$key = paste(locs_disturbed_pid$Lon, locs_disturbed_pid$Lat, locs_disturbed_pid$scenario)
        
        index = match(locs_disturbed_pid$key, d_var$key)
        d_var = d_var[index,]
        d_var$key = NULL
      }
      
      assign(paste0("d_", var), d_var)
    }
    
    ### Ecological PFT-dependent variables
    
    # Cut data to 100 years of recovery
    d_eco = d_eco %>%
      filter(age <= 100)
    
    # Get time-independent variables
    d_eco_ind = d_eco %>%
      distinct(Lon,Lat,PFT,initial_recruitment, recruitment_ten_years, previous_state, Scenario) %>%
      group_by(PFT) %>%
      ungroup()
    
    # Spread the data so that each PFT category’s values are in separate columns
    d_eco_wide <- d_eco_ind %>%
      pivot_wider(names_from = PFT, values_from = c(initial_recruitment, recruitment_ten_years, previous_state))
    
    # Sort the resulting data according to right order of locations
    
    d_eco_wide$key = paste(d_eco_wide$Lon, d_eco_wide$Lat, d_eco_wide$Scenario)
    
    index = match(locs_disturbed_pid$key, d_eco_wide$key)
    d_eco_wide = d_eco_wide[index,]
    d_eco_wide$key = NULL
    
    ### Climate covariates
    d_climate = rbind(d_climate_picontrol,d_climate_ssp126, d_climate_ssp370, d_climate_ssp585)
    d_climate = d_climate %>%
      filter(age <= 100 & age > 0)
    d_climate = d_climate[!duplicated(d_climate),]
    
    # store Lon/Lat values
    locs_climate = d_climate %>% distinct(Lon,Lat,Scenario)
    
    ### Get wider format
    names.dfs = c("temp","temp_min", "temp_max", "precip")
    
    names.vars = c("tas_yearlymeam", "tas_yearlymin", "tas_yearlymax", "pr_yearlysum")
    
    for (i in 1:4) {
      var <- names.vars[i]
      
      d_climate_wide <- d_climate %>%
        distinct(Lon, Lat, age, !!sym(var), Scenario) %>%
        pivot_wider(names_from = c(Lon, Lat, Scenario), values_from = !!sym(var)) %>%
        arrange(age) %>%
        rename_with(~ gsub("_", "/", .), everything()) 
      
      assign(paste0("d_climate_", names.dfs[i], "_wide"), d_climate_wide)
      
      # Get funData objects
      funData_obj = funData(argvals = 1:100, X = t(d_climate_wide[,-1]))
      assign(paste0("funData_", names.dfs[i]), funData_obj)
      saveRDS(funData_obj, paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/AllPatches/funData_",names.dfs[i],"_PID", pid, ".rds"))
    }
    ## Run MFPCA
    uniExpansions <- list(list(type = "uFPCA", npc = 5),
                          list(type = "uFPCA", npc = 5),
                          list(type = "uFPCA", npc = 5),
                          list(type = "uFPCA", npc = 5))
    MFPCA_climate = MFPCA_2(multiFunData(funData_temp, funData_temp_min, funData_temp_max, funData_precip),
                            M = 3, uniExpansions = uniExpansions, fit = TRUE)
    
    saveRDS(MFPCA_climate, paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/AllPatches/MFPCA_climate_PID", pid, ".rds"))
    
    # Reorder scores according to right order of locations
    locs_climate$key = paste(locs_climate$Lon, locs_climate$Lat, locs_climate$Scenario)
    locs_disturbed_pid$key = paste(locs_disturbed_pid$Lon, locs_disturbed_pid$Lat, locs_disturbed_pid$scenario)
    
    index = match(locs_disturbed_pid$key, locs_climate$key)
    
    scores_climate = as.data.frame(MFPCA_climate$scores)
    scores_climate = scores_climate[index,]
    locs_climate$key = NULL
    
    ### Get full data set
    write.table(d_soil, paste0("Scripts/MA_FDA_veg/04_Model/Data/AllPatches/d_soil_PID", pid, ".txt"))
    write.table(d_eco_wide, paste0("Scripts/MA_FDA_veg/04_Model/Data/AllPatches/d_eco_wide_PID", pid, ".txt"))
    write.table(scores_climate, paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/AllPatches/scores_climate_PID", pid, ".txt"))
    
  }
}

############################### Models for each PID ############################

# Define the generalized logistic function
generalized_logistic <- function(x, A, k, x0, C) {
  A / (1 + exp(-k * (x - x0))) + C
}

# Inverse of the generalized logistic function
inverse_generalized_logistic <- function(y, A, k, x0, C) {
  x0 - (1 / k) * log((A / (y - C)) - 1)
}

# Create df for storing the coefficients
if (whichScen != "all" || !covScen) row_nr = 3
df_beta = as.data.frame(matrix(0, nrow = (16-row_nr)*2, ncol = 25))
colnames(df_beta) = 0:24

for (pid in c(0:24)){#c(0:6, 9, 11:13, 15, 18:21, 23:24)){
  print(paste0("Start with PID ", pid))
  
  # Get target variable: MFPCA scores
  MFPCA_pid = readRDS(paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/AllPatches/MFPCA_PID", pid, ".rds"))
  scores_pid = as.data.frame(MFPCA_pid$scores)
  colnames(scores_pid) = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
  
  # Get soil and ecological variables
  d_soil = read.table(paste0("Scripts/MA_FDA_veg/04_Model/Data/AllPatches/d_soil_PID", pid, ".txt"))
  d_eco_wide = read.table(paste0("Scripts/MA_FDA_veg/04_Model/Data/AllPatches/d_eco_wide_PID", pid, ".txt"))
  
  # Get MFPCA results for temp (mean, min, max) and precip
  scores_climate = read.table(paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/AllPatches/scores_climate_PID", pid, ".txt"))

  # Combine data set
  combined_d = cbind(scores_pid[,c(1:2)], d_soil, d_eco_wide[,c(4:18)], 
                     scores_climate[,1])  %>% 
    mutate(Scenario = as.factor(Scenario))
  
  colnames(combined_d)[c(1:2,28)] = c("PC1", "PC2", "PC1_climate")
  
  if(whichScen != "all") {
    combined_d = combined_d[combined_d$Scenario == whichScen,]
    covScen = FALSE
  }
  ################################### Model ######################################
  
  # Produce train-/test- split
  set.seed(1)
  index = sample.int(nrow(combined_d), size = ceiling(0.8 * nrow(combined_d)))
  d_train = combined_d[index,]
  d_test = combined_d[-index,]
                    
  ### Estimate the parameters of the generalized logistic function on the whole data set
  
  params = c("A" = 12.21, "k" = 0.72, "x0" = 0.3, "C" = -6.03)

  d_train_trafo = d_train %>%
    mutate(PC1_trafo = inverse_generalized_logistic(PC1, A = params["A"], k = params["k"], x0 = params["x0"], C = params["C"]))
  
  formula_PC1 <- PC1_trafo ~ s(Lon,Lat) + 
    sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
    time_since_dist +
    initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
    PC1_climate
  
  formula_PC2 <- PC2 ~ s(Lon,Lat) + 
    sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
    time_since_dist +
    initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
    PC1_climate
  
  # Conditionally add the Scenario covariate if covScen is TRUE
  if (covScen) {
    formula_PC1 <- update(formula_PC1, . ~ . + Scenario)
    formula_PC2 <- update(formula_PC2, . ~ . + Scenario)
  }
  
  # Fit the model with the dynamically constructed formula
  gam_trafo <- gam(list(formula_PC1, formula_PC2), 
                   data = d_train_trafo, 
                   family = mvn(d = 2))
 
  df_beta[,pid+1] = gam_trafo$coefficients[c(1:(16-row_nr),(46-row_nr):(61 - 2*row_nr))]
}
rownames(df_beta) = names(gam_trafo$coefficients[c(1:(16-row_nr),(46-row_nr):(61 - 2*row_nr))])

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

# Apply the function row-wise
ci_data <- as.data.frame(t(apply(df_beta, 1, calculate_ci)))

# Create a confidence interval plot for both PC scores
ci_data$Effects <- ifelse((ci_data$lower > 0 & ci_data$upper > 0), "Positive",
                               ifelse((ci_data$lower < 0 & ci_data$upper < 0), "Negative","Non-significant"))
ci_data$Effects <- factor(ci_data$Effects, levels = c("Positive", "Negative", "Non-significant"))
ci_data$parameter <- factor(rownames(df_beta), levels = rownames(df_beta))

g = ggplot(ci_data[c(1:(16-row_nr)),], aes(y = parameter, x = mean)) +
  geom_point(aes(color = Effects), size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = Effects), height = 0.2) +  # Horizontal error bars
  scale_color_manual(values = c("Positive" = "darkred", "Negative" = "cornflowerblue", "Non-significant" = "grey")) +  
  labs(title = "95% Confidence intervals - PC 1",
       y = "Predictors",
       x = "Mean and 95% confidence interval") + theme_bw() +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black")

if (!covScen) g = g + ggtitle(paste0("95% Confidence intervals - PC 1 - " , whichScen))
g

ggsave(paste0("Scripts/Plots/Model/AllPatches/pdf/CIs_patches_PC1_", whichScen, covScen, ".pdf"), g, width = 7, height = 4.3)
ggsave(paste0("Scripts/Plots/Model/AllPatches/png/CIs_patches_PC1_", whichScen, covScen, ".png"), g, width = 7, height = 4.3)

ci_data_PC2 = ci_data[(17-row_nr):(32-2*row_nr),]
ci_data_PC2$parameter <- factor(rownames(df_beta)[1:(16-row_nr)], levels = rownames(df_beta)[1:(16-row_nr)])

g = ggplot(ci_data_PC2, aes(y = parameter, x = mean)) +
  geom_point(aes(color = Effects), size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = Effects), height = 0.2) +  # Horizontal error bars
  scale_color_manual(values = c("Positive" = "darkred", "Negative" = "cornflowerblue", "Non-significant" = "grey")) +  
  labs(title = "95% Confidence intervals - PC 2",
       y = "Predictors",
       x = "Mean and 95% confidence interval") + theme_bw() +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black")

if (!covScen) g = g + ggtitle(paste0("95% Confidence intervals - PC 2 - " , whichScen))
g

ggsave(paste0("Scripts/Plots/Model/AllPatches/pdf/CIs_patches_PC2_", whichScen, covScen,  ".pdf"), g, width = 7, height = 4.3)
ggsave(paste0("Scripts/Plots/Model/AllPatches/png/CIs_patches_PC2_", whichScen, covScen,  ".png"), g, width = 7, height = 4.3)

### Calculate estimated t-values from df_beta
ci_data$t = apply(df_beta, 1, function(x)  sqrt(length(x)) * mean(as.numeric(x)) / sd(as.numeric(x)))

g = ggplot(ci_data[c(1:(16-row_nr)),], aes(x = t, y = parameter)) +
  geom_bar(stat = "identity", aes(fill = t > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(text = element_text(size = 15),
                     legend.position = "none",
                     plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(x = "t-Value", y = "Predictors") + 
  ggtitle("t-Values of PC 1 using all patches") +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", color = "black")

if (!covScen) g = g + ggtitle(paste0("t-Values of PC 1 - ", whichScen)) + xlim(-28,75)
g

ggsave(paste0("Scripts/Plots/Model/AllPatches/pdf/t-values_patches_PC1_", whichScen, covScen,  ".pdf"), g, width = 7, height = 4.3)
ggsave(paste0("Scripts/Plots/Model/AllPatches/png/t-values_patches_PC1_", whichScen, covScen,  ".png"), g, width = 7, height = 4.3)

ci_data_PC2 = ci_data[(17-row_nr):(32-2*row_nr),]
ci_data_PC2$parameter <- factor(rownames(df_beta)[1:(16-row_nr)], levels = rownames(df_beta)[1:(16-row_nr)])

g = ggplot(ci_data_PC2, aes(x = t, y = parameter)) +
  geom_bar(stat = "identity", aes(fill = t > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(text = element_text(size = 15),
                     legend.position = "none",
                     plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(x = "t-Value", y = "Predictors") + 
  ggtitle("t-Values of PC 2 using all patches") +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", color = "black") 

if (!covScen) g = g + ggtitle(paste0("t-Values of PC 2 - ", whichScen)) + xlim(-12.5,12.5)
g

ggsave(paste0("Scripts/Plots/Model/AllPatches/pdf/t-values_patches_PC2_", whichScen, covScen, ".pdf"), g, width = 7, height = 4.3)
ggsave(paste0("Scripts/Plots/Model/AllPatches/png/t-values_patches_PC2_", whichScen, covScen, ".png"), g, width = 7, height = 4.3)

