################################################################################
############################ Master's Thesis ###################################
################################################################################

############# Model: Get absolute carbon values and MFPCA scores ###############

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

## Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)
scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")
pfts = c("BNE", "IBS", "otherC", "TeBS","Tundra")
set.seed(1)

## Choose parameters:
start_year = 2015
end_year = 2040
#pft = "TeBS"      # Choose from Tundra, BNE, IBS, otherC, TeBS
pid = 1           # Choose pid (int) or 'all'
k = 4             # Number of clusters for kmeans algorithm
M = 10            # Number of PCs
createFunData = FALSE
runMFPCA = FALSE
# scen = "picontrol"
locs_disturbed = read.table("Scripts/Plots/MFPCA/plot_data/locs_disturbed.txt")

#######################
## Get 450 different colors for plotting

palettes <- c("Set1", "Set2", "Set3", "Dark2", "Paired", "Accent")
all_colors <- character()
for (palette_name in palettes) {
  palette <- brewer.pal(20, palette_name)
  all_colors <- c(all_colors, palette)
}
palette_450 <- unique(all_colors[1:470])
pal <- colorRampPalette(palette_450)

#######################
# Create funData object 
if(createFunData){
  for (pft in pfts){
    print(paste("Start with PFT", long_names_pfts(tolower(pft))))
    
    for (scen in scenarios){
      print(long_names_scenarios(scen))
      d_scen = get_data_fpca_cmass(scen, start_year, end_year,pid,pft)
      assign(paste0("d_",scen), d_scen)
      print("...done.")
    }
    d_pft = rbind(t(d_picontrol[[1]][,-1]), t(d_ssp126[[1]][,-1]), t(d_ssp370[[1]][,-1]), t(d_ssp585[[1]][,-1]))
    assign(paste0("d_", pft), d_pft)
  }
  
  # Create univariate funData objects
  for (pft in pfts) {
    fun_pft <- funData(argvals = 1:100,
                       X = get(paste0("d_", pft)))
    dimnames(fun_pft@X)[[2]] = 1:100
    assign(paste0("fun_", pft), fun_pft)
    saveRDS(fun_pft, paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_", pft, "_cmass.rds"))
  }
  
  # Create multiFunData object
  multiFun_pft = multiFunData(fun_BNE, fun_IBS, fun_otherC, fun_TeBS, fun_Tundra)
  
  saveRDS(multiFun_pft, "Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_cmass.rds")
  
  print("... all done.")
}

# Run MFPCA
if (runMFPCA){
  multiFun_pft = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_cmass.rds")
  
  uniExpansions <- list(list(type = "uFPCA", npc = 10),
                        list(type = "uFPCA", npc = 10),
                        list(type = "uFPCA", npc = 10),
                        list(type = "uFPCA", npc = 10),
                        list(type = "uFPCA", npc = 10))
  
  MFPCA_all <- MFPCA_2(multiFun_pft, M = M, fit = TRUE, uniExpansions = uniExpansions)
  saveRDS(MFPCA_all, paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_all_cmass.rds"))
}

# Get multiFunData objects
multiFun_pft = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_cmass.rds")

# Get MFPCA results
MFPCA_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_all_cmass.rds")

plot_data_all = as.data.frame(MFPCA_all$scores) %>% 
  mutate(scenario = c(rep("Control", 434),rep("SSP1-RCP2.6", 442),rep("SSP3-RCP7.0", 462),rep("SSP5-RCP8.5", 465)))

write.table(plot_data_all, "Scripts/Plots/Model/plot_data/plot_data_cmass.txt", row.names = TRUE, col.names = TRUE)
