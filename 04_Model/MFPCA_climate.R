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
library(rlang)

scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")

#######################
## Get different colors for plotting

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
             tas_yearlymeam, tas_yearlymin, tas_yearlymax, pr_yearlysum) %>%
    mutate(Scenario = long_names_scenarios(scen))
  
  d_nuptake_PFT = d_scen_2015_2040 %>%
    filter(PID == pid) %>%
    distinct(Lon, Lat, age, PFT, Nuptake) %>%
    mutate(Scenario = long_names_scenarios(scen))
  
  assign(paste0("d_climate_", scen), d_climate_scen)
  assign(paste0("d_nuptake_PFT_", scen), d_nuptake_PFT)
}

d_climate = rbind(d_climate_picontrol,d_climate_ssp126, d_climate_ssp370, d_climate_ssp585)


d_nuptake_PFT = rbind(d_nuptake_PFT_picontrol,
                      d_nuptake_PFT_ssp126, 
                      d_nuptake_PFT_ssp370, 
                      d_nuptake_PFT_ssp585)


############################### Perform FPCA ##################################
## Get functional fit
# Create funData object(s)

# Cut data to 100 years of recovery
d_climate = d_climate %>%
  filter(age <= 100 & age > 0)

d_climate = d_climate[!duplicated(d_climate),]

d_nuptake_PFT = d_nuptake_PFT %>%
  filter(age <= 100 & age > 0)

d_nuptake_PFT = d_nuptake_PFT[!duplicated(d_nuptake_PFT),]

d_nuptake_PFT_wide <- d_nuptake_PFT %>%
  pivot_wider(names_from = PFT, values_from = c(Nuptake))

# store Lon/Lat values
locs_climate = d_climate %>% distinct(Lon,Lat,Scenario)

### Get wider format
names.dfs = c("temp", "temp_min", "temp_max", "precip", "nuptake_total",
              "nuptake_BNE", "nuptake_IBS", "nuptake_otherC", "nuptake_TeBS", "nuptake_Tundra")

names.vars = c("tas_yearlymeam", "tas_yearlymin", "tas_yearlymax", "pr_yearlysum", "Nuptake_total",
               "BNE", "IBS", "otherC", "TeBS", "Tundra")

for (i in 1:5) {
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
  #saveRDS(funData_obj, paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/funData_",names.dfs[i],".rds"))
  
  ## Run MFPCA
  MFPCA_climate = MFPCA_2(multiFunData(funData_obj), M = 6, uniExpansions = list(list(type = "uFPCA", npc = 10)), fit = TRUE)
  if (i == 1 || i == 4){
    MFPCA_climate$functions@.Data[[1]]@X = - MFPCA_climate$functions@.Data[[1]]@X
    MFPCA_climate$scores = - MFPCA_climate$scores
  }
  assign(paste0("MFPCA_", names.dfs[i]), MFPCA_climate)
  
  
  #saveRDS(MFPCA_climate, paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/MFPCA_", names.dfs[i], ".rds"))
  print(var)
  print(round(MFPCA_climate$values[1:10]/sum(MFPCA_climate$values) *100,3))
}

for (i in 6:10) {
  var <- names.vars[i]

  d_nuptake_PFT <- t(d_nuptake_PFT_wide %>%
                                distinct(Lon, Lat, age, !!sym(var), Scenario) %>%
                                pivot_wider(names_from = c(Lon, Lat, Scenario), values_from = !!sym(var)) %>%
                                arrange(age) %>%
                                rename_with(~ gsub("_", "/", .), everything()) %>%
                                select(-age))
  
  
  assign(paste0("d_nuptake_", names.dfs[i], "_wide"), d_nuptake_PFT)
  
  # Get funData objects
  funData_obj = funData(argvals = 1:100, X = d_nuptake_PFT)
  assign(paste0("funData_", names.dfs[i]), funData_obj)
  #saveRDS(funData_obj, paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/funData_",names.dfs[i],".rds"))
  
  ## Run MFPCA
  MFPCA_climate = MFPCA_2(multiFunData(funData_obj), M = 6, uniExpansions = list(list(type = "uFPCA", npc = 10)), fit = TRUE)
  assign(paste0("MFPCA_", names.dfs[i]), MFPCA_climate)
  
  #saveRDS(MFPCA_climate, paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/MFPCA_", names.dfs[i], ".rds"))
  print(var)
  print(round(MFPCA_climate$values[1:10]/sum(MFPCA_climate$values) *100,3))
}

############# Plot the functional fits including mean functions ################

bounds = t(matrix(c(1,434,435,876,877,1338,1339,1803), nrow = 2))

ylabs = c("Annual mean temperature in K", "Annual minimum temperature in K", "Annual maximum temperature in K",  
          bquote("Annual precipitation in kg/m"^2), bquote("Total nitrogen uptake per grid cell in g/m"^2),
          bquote("Nitrogen uptake per grid cell in g/m"^2), bquote("Nitrogen uptake per grid cell in g/m"^2),
          bquote("Nitrogen uptake per grid cell in g/m"^2), bquote("Nitrogen uptake per grid cell in g/m"^2),
          bquote("Nitrogen uptake per grid cell in g/m"^2))

mains = c("Annual average temperature", "Annual minimum temperature", "Annual maximum temperature",
          "Annual precipitation", "Total nitrogen uptake per grid cell", "Nitrogen uptake per grid cell - BNE",
          "Nitrogen uptake per grid cell - IBS", "Nitrogen uptake per grid cell - otherC", "Nitrogen uptake per grid cell - TeBS",
          "Nitrogen uptake per grid cell - Tundra")

ylims = t(matrix(c(263, 284, 210, 275, 284, 308, 0, 4000, 0, 120, 0, 45, 0, 100, 0, 42, 0, 110, 0, 68), nrow = 2))

for (i in seq_along(names.dfs)) {
  var <- names.dfs[i]
  funData_var = get(paste0("funData_", var))
  MFPCA_var = get(paste0("MFPCA_", var))
  
  ## As pdf
  pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/BasisRep/BasisRep_", var, ".pdf"), width = 5, height = 4.5)
  par(mfrow=c(1,1))
  plot(funData_var, col = pal(1803)[1:1803],
       xlim = c(1,100), type = 'l',
       cex.main = 1.5, cex = 1, cex.lab = 1,
       xlab = "Year after Disturbance",
       ylab = ylabs[i], 
       main = paste(mains[i]))
  lines(x=1:100, y = MFPCA_var$meanFunction@.Data[[1]]@X, col = "black", lwd = 3)
  dev.off()

  pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/BasisRep/BasisRep_", var, "_scen.pdf"), width = 10, height = 11)
  par(mfrow =(c(2,2)))
  plot(funData_var[1:434], col = pal(434)[1:434],
       xlim = c(1,100), ylim = c(ylims[i,1],ylims[i,2]),
       type = 'l',
       cex.main = 1.5, cex = 1.3,  cex.lab = 1.5,
       xlab = "Year after Disturbance",
       ylab = ylabs[i],
       main = paste("Control") )
  lines(x=1:100, y = colMeans(MFPCA_var$fit@.Data[[1]]@X[bounds[1,1]:bounds[1,2],]), col = "black", lwd = 3)

  plot(funData_var[bounds[2,1]:bounds[2,2]], col = pal(500)[1:500],
       xlim = c(1,100), ylim = c(ylims[i,1],ylims[i,2]),
       type = 'l',
       cex.main = 1.5, cex = 1.3,  cex.lab = 1.5,
       xlab = "Year after Disturbance",
       ylab = ylabs[i],
       main = paste("SSP1-RCP2.6") )
  lines(x=1:100, y = colMeans(MFPCA_var$fit@.Data[[1]]@X[bounds[2,1]:bounds[2,2],]), col = "black", lwd = 3)

  plot(funData_var[bounds[3,1]:bounds[3,2]], col = pal(500)[1:500],
       xlim = c(1,100), ylim = c(ylims[i,1],ylims[i,2]),
       type = 'l',
       cex.main = 1.5, cex = 1.3,  cex.lab = 1.5,
       xlab = "Year after Disturbance",
       ylab = ylabs[i],
       main = paste("SSP3-RCP7.0") )
  lines(x=1:100, y = colMeans(MFPCA_var$fit@.Data[[1]]@X[bounds[3,1]:bounds[3,2],]), col = "black", lwd = 3)

  plot(funData_var[bounds[4,1]:bounds[4,2]], col = pal(500)[1:500],
       xlim = c(1,100), ylim = c(ylims[i,1],ylims[i,2]),
       type = 'l',
       cex.main = 1.5, cex = 1.3,  cex.lab = 1.5,
       xlab = "Year after Disturbance",
       ylab = ylabs[i],
       main = paste("SSP5-RCP8.5") )
  lines(x=1:100, y = colMeans(MFPCA_var$fit@.Data[[1]]@X[bounds[4,1]:bounds[4,2],]), col = "black", lwd = 3)

  mtext(paste(mains[i]), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.5)
  dev.off()

  ## As png
  png(paste0("Scripts/Plots/Model/(M)FPCA/png/BasisRep/BasisRep_", var, ".png"), width = 800, height = 650)
  par(mfrow=c(1,1))
  plot(funData_var, col = pal(1803)[1:1803],
       xlim = c(1,100), type = 'l',
       cex.main = 1.5, cex = 1.3,
       xlab = "Year after Disturbance",
       ylab = ylabs[i],
       main = paste(mains[i]))
  lines(x=1:100, y = MFPCA_var$meanFunction@.Data[[1]]@X, col = "black", lwd = 3)
  dev.off()

  png(paste0("Scripts/Plots/Model/(M)FPCA/png/BasisRep/BasisRep_", var, "_scen.png"), width = 1000, height = 1100)
  par(mfrow =(c(2,2)))
  plot(funData_var[1:434], col = pal(434)[1:434],
       xlim = c(1,100), ylim = c(ylims[i,1],ylims[i,2]),
       type = 'l',
       cex.main = 1.5, cex = 1.3,
       xlab = "Year after Disturbance",
       ylab = ylabs[i],
       main = paste("Control") )
  lines(x=1:100, y = colMeans(MFPCA_var$fit@.Data[[1]]@X[bounds[1,1]:bounds[1,2],]), col = "black", lwd = 3)

  plot(funData_var[bounds[2,1]:bounds[2,2]], col = pal(500)[1:500],
       xlim = c(1,100), ylim = c(ylims[i,1],ylims[i,2]),
       type = 'l',
       cex.main = 1.5, cex = 1.3,
       xlab = "Year after Disturbance",
       ylab = ylabs[i],
       main = paste("SSP1-RCP2.6") )
  lines(x=1:100, y = colMeans(MFPCA_var$fit@.Data[[1]]@X[bounds[2,1]:bounds[2,2],]), col = "black", lwd = 3)

  plot(funData_var[bounds[3,1]:bounds[3,2]], col = pal(500)[1:500],
       xlim = c(1,100), ylim = c(ylims[i,1],ylims[i,2]),
       type = 'l',
       cex.main = 1.5, cex = 1.3,
       xlab = "Year after Disturbance",
       ylab = ylabs[i],
       main = paste("SSP3-RCP7.0") )
  lines(x=1:100, y = colMeans(MFPCA_var$fit@.Data[[1]]@X[bounds[3,1]:bounds[3,2],]), col = "black", lwd = 3)

  plot(funData_var[bounds[4,1]:bounds[4,2]], col = pal(500)[1:500],
       xlim = c(1,100), ylim = c(ylims[i,1],ylims[i,2]),
       type = 'l',
       cex.main = 1.5, cex = 1.3,
       xlab = "Year after Disturbance",
       ylab = ylabs[i],
       main = paste("SSP5-RCP8.5") )
  lines(x=1:100, y = colMeans(MFPCA_var$fit@.Data[[1]]@X[bounds[4,1]:bounds[4,2],]), col = "black", lwd = 3)

  mtext(paste(mains[i]), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.5)
  dev.off()
}

################################ Plot the PCs ##################################
# mean temperature
# PC 1

for (i in seq_along(names.dfs)) {
  var <- names.dfs[i]
  MFPCA_var = get(paste0("MFPCA_", var))
  
  for (iPC in c(1:5)){
    pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/PCs/PC", iPC, "_", var, ".pdf"), width = 6.5, height = 6.5)
    plot(MFPCA_var, combined = TRUE, plotPCs = iPC, xlim = c(1,100), 
         xlab = "Year after Disturbance",
         ylab = ylabs[i], 
         cex.main = 1.5, cex = 1.4)
    dev.off()
    
    png(paste0("Scripts/Plots/Model/(M)FPCA/png/PCs/PC", iPC, "_", var, ".png"), width = 650, height = 650)
    plot(MFPCA_var, combined = TRUE, plotPCs = iPC, xlim = c(1,100), 
         xlab = "Year after Disturbance",
         ylab = ylabs[i], 
         cex.main = 1.5, cex = 1.4)
    dev.off()
  }
  
}

############################ Plot reconstructions ##############################
for (i in seq_along(names.dfs)) {
  var <- names.dfs[i]
  MFPCA_var = get(paste0("MFPCA_", var))
  
  # As pdf
  pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/Reconstruction/FPCA_reconstruct_", var, ".pdf"), width = 8.5, height = 6.5)
  plot(MFPCA_var$fit, col = pal(1803)[1:1803],
       xlim = c(1,100), type = 'l',
       cex.main = 2, cex = 1.3, cex.lab = 1.8,
       xlab = "Year after Disturbance",
       ylab = ylabs[i], 
       main = paste("Reconstructed fit using 6 PCs"))
  dev.off()
  
  iScen = 0
  for (scen in scenarios){
    pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/Reconstruction/FPCA_reconstruct_", scen, "_", var, ".pdf"), width = 8, height = 6.5)
    
    iScen = iScen + 1
    nOb = bounds[iScen,2] - bounds[iScen,1] + 1
  
    plot(MFPCA_var$fit, obs = c(bounds[iScen,1]: bounds[iScen,2]), col = pal(nOb)[1:nOb],
         xlim = c(1,100), type = 'l',
         cex.main = 1.8, cex = 1.5,
         xlab = "Year after Disturbance",
         ylab = ylabs[i], 
         main = paste("Reconstructed fit using 6 PCs -", long_names_scenarios(scen)) )
    dev.off()
  }
  
  # As png
  png(paste0("Scripts/Plots/Model/(M)FPCA/png/Reconstruction/FPCA_reconstruct_", var, ".png"), width = 800, height = 650)
  plot(MFPCA_var$fit, col = pal(1803)[1:1803],
       xlim = c(1,100), type = 'l',
       cex.main = 1.8, cex = 1.5,
       xlab = "Year after Disturbance",
       ylab = ylabs[i], 
       main = paste("Reconstructed fit using 6 PCs"))
  dev.off()
  
  iScen = 0
  for (scen in scenarios){
    png(paste0("Scripts/Plots/Model/(M)FPCA/png/Reconstruction/FPCA_reconstruct_", scen, "_", var, ".png"), width = 800, height = 650)
    
    iScen = iScen + 1
    nOb = bounds[iScen,2] - bounds[iScen,1] + 1
    
    plot(MFPCA_var$fit, obs = c(bounds[iScen,1]: bounds[iScen,2]), col = pal(nOb)[1:nOb],
         xlim = c(1,100), type = 'l',
         cex.main = 1.8, cex = 1.5,
         xlab = "Year after Disturbance",
         ylab = ylabs[i], 
         main = paste("Reconstructed fit using 6 PCs -", long_names_scenarios(scen)) )
    dev.off()
  }
}

############################## Store scores ####################################
# Reorder scores according to right order of locations
locs_climate$key = paste(locs_climate$Lon, locs_climate$Lat, locs_climate$Scenario)
locs_disturbed$key = paste(locs_disturbed$Lon, locs_disturbed$Lat, locs_disturbed$Scenario)

index = match(locs_disturbed$key, locs_climate$key)

for (i in seq_along(names.dfs)) {
  var <- names.dfs[i]
  MFPCA_var = get(paste0("MFPCA_", var))
  
  scores_var = as.data.frame(MFPCA_var$scores)
  scores_var = scores_var[index,]
  write.table(scores_var, paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_",var, ".txt"))
  
}
locs_climate$key = NULL


################################################################################
############################# Run MFPCA for climate data #######################
################################################################################
# For all scenarios together
multiFun_climate = multiFunData(funData_temp, funData_temp_min, funData_temp_max, funData_precip)
saveRDS(multiFun_climate, "Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/multiFun_climate.rds")

uniExpansions <- list(list(type = "uFPCA", npc = 5),
                      list(type = "uFPCA", npc = 5),
                      list(type = "uFPCA", npc = 5),
                      list(type = "uFPCA", npc = 5))

MFPCA_climate <- MFPCA_2(multiFun_climate, M = 3, fit = TRUE, uniExpansions = uniExpansions)
assign(paste0("MFPCA_climate"), MFPCA_climate)
saveRDS(MFPCA_climate, paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/MFPCA_climate.rds"))

# For each scenario separately
i = 0
for (scen in scenarios){
  set.seed(1)
  multiFun_climate = multiFunData(funData_temp, funData_temp_min, funData_temp_max, funData_precip)
  i = i+1
  for (iClimate in 1:4){
    multiFun_climate@.Data[[iClimate]]@X = multiFun_climate@.Data[[iClimate]]@X[bounds[i,1]:bounds[i,2],]
    dimnames(multiFun_climate@.Data[[iClimate]]@X) = dimnames(multiFun_climate@.Data[[iClimate]]@X[bounds[i,1]:bounds[i,2]])
  }
  uniExpansions <- list(list(type = "uFPCA", npc = 5),
                        list(type = "uFPCA", npc = 5),
                        list(type = "uFPCA", npc = 5),
                        list(type = "uFPCA", npc = 5))
  
  MFPCA_climate <- MFPCA_2(multiFun_climate, M = 3, fit = TRUE, uniExpansions = uniExpansions)
  
  # Flip the sign for control and SSP1-RCP2.6 for consistency
  if (scen == "picontrol" | scen == "ssp126"){
    MFPCA_climate$scores = - MFPCA_climate$scores
    for (iClimate in 1:4){
      MFPCA_climate$functions@.Data[[iClimate]]@X = - MFPCA_climate$functions@.Data[[iClimate]]@X
    }
  }
  assign(paste0("MFPCA_climate_", scen), MFPCA_climate)
  saveRDS(MFPCA_climate,paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/MFPCA_climate_", scen, ".rds"))
  
  # Save the scores
  scores = MFPCA_climate$scores
  write.table(scores, paste0("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_climate_", scen, ".txt"))
}

################################ Plot the PCs ##################################
for (iPC in 1:3){
  pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/PCs/PCs_climate_PC", iPC, ".pdf"), width = 12, height = 5)
  plot.MFPCAfit_2(MFPCA_climate, plotPCs = iPC, combined = TRUE, xlab = paste("Year after Disturbance"), 
                  xlim = c(1,100), ylab = "Values", cex.main = 2.3, cex.lab = 1.8,
                  main = c("Mean temperature", "Minimum temperature", "Maximum temperature", "Precipitation"))
  dev.off()
  
  png(paste0("Scripts/Plots/Model/(M)FPCA/png/PCs/PCs_climate_PC", iPC, ".png"), width = 1200, height = 500)
  plot.MFPCAfit_2(MFPCA_climate, plotPCs = iPC, combined = TRUE, xlab = paste("Year after Disturbance"), 
                  xlim = c(1,100), ylab = "Values", cex.main = 2.3, cex.lab = 1.7,
                  main = c("Mean temperature", "Minimum temperature", "Maximum temperature", "Precipitation"))
  dev.off()
  
}

# Variance that is accounted for
round(100*MFPCA_climate$values/(sum(MFPCA_climate$values)),2)

# Save the scores
scores = as.data.frame(MFPCA_climate$scores)
colnames(scores) = c("PC1", "PC2", "PC3")
write.table(scores, "Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_climate.txt")

### Plot reconstructions using 3 PCs
ylabs = c("Annual mean temperature in K", "Annual minimum temperature in K", "Annual maximum temperature in K",  
          bquote("Annual precipitation in kg/m"^2))

mains = c("Annual average temperature", "Annual minimum temperature", "Annual maximum temperature",
          "Annual precipitation")
for (i in 1:4){
  # As pdf
  pdf(paste0("Scripts/Plots/Model/(M)FPCA/pdf/Reconstruction/FPCA_reconstruct_climate_", names.dfs[i], ".pdf"), width = 5, height = 4.5)
  plot(MFPCA_climate$fit[[i]], col = pal(1803)[1:1803],
       xlim = c(1,100), type = 'l',
       cex.main = 1.5, cex = 1,, cex.lab = 1,
       xlab = "Year after Disturbance",
       ylab = ylabs[i], 
       main = paste("Reconstructed fit using 3 PCs"))
  dev.off()
  
  # As png
  png(paste0("Scripts/Plots/Model/(M)FPCA/png/Reconstruction/FPCA_reconstruct_climate_", names.dfs[i], ".png"), width = 800, height = 650)
  plot(MFPCA_climate$fit[[i]], col = pal(1803)[1:1803],
       xlim = c(1,100), type = 'l',
       cex.main = 1.5, cex = 1,, cex.lab = 1,
       xlab = "Year after Disturbance",
       ylab = ylabs[i], 
       main = paste("Reconstructed fit using 3 PCs"))
  dev.off()
}

