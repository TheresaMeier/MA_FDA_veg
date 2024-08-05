################################################################################
############################ Master's Thesis ###################################
################################################################################

##################### Model: Transformation of the target ######################

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

scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")
pfts = c("BNE", "IBS", "otherC", "TeBS","Tundra")

## Choose parameters:
pid = 1           # Choose pid (int) or 'all'
M = 10            # Number of PCs
eps = 1e-08

locs_disturbed = read.table("Scripts/Plots/MFPCA/plot_data/locs_disturbed.txt")

#######################
############################ Create funData object #############################
  
multiFun_pft = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_1803.rds")

# First attempt: suppose that sum(e^(x_i)) = 0 --> x_i = log(p_i)
multiFun_pft_log <- lapply(multiFun_pft@.Data, function(obj) {
  x <- obj@X
  #x[x == 0] <- eps
  obj@X <- x
  return(obj)
})

multiFun_pft_log = multiFunData(multiFun_pft_log)

for (i in 1:5){
  multiFun_pft_log@.Data[[i]]@X = log(1+multiFun_pft_log@.Data[[i]]@X)
}
saveRDS(multiFun_pft_log,"Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_pft_log.rds")

# Second attempt: Choose Tundra as reference category and set it equal to zero
multiFun_pft_ref = multiFun_pft
for (i in 1:4){
  multiFun_pft_ref@.Data[[i]]@X = log((multiFun_pft_ref@.Data[[i]]@X+eps)/(multiFun_pft_ref@.Data[[5]]@X+eps))
}
multiFun_pft_ref@.Data[[5]]@X = matrix(0, nrow = 1803, ncol = 100)
dimnames(multiFun_pft_ref@.Data[[5]]@X) = dimnames(multiFun_pft_ref@.Data[[4]]@X)
saveRDS(multiFun_pft_ref,"Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_pft_ref.rds")


################################## Run MFPCA ###################################

multiFun_pft_log = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_pft_log.rds")
multiFun_pft_ref = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_pft_ref.rds")

uniExpansions <- list(list(type = "uFPCA", npc = 10),
                      list(type = "uFPCA", npc = 10),
                      list(type = "uFPCA", npc = 10),
                      list(type = "uFPCA", npc = 10),
                      list(type = "uFPCA", npc = 10))

MFPCA_log <- MFPCA_2(multiFun_pft_log, M = M, fit = TRUE, uniExpansions = uniExpansions)
saveRDS(MFPCA_log, paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_log.rds"))

MFPCA_ref <- MFPCA_2(multiFun_pft_ref, M = M, fit = TRUE, uniExpansions = uniExpansions)
# Doesn't work because the Tundra values are all 0!

# saveRDS(MFPCA_ref, paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_ref.rds"))

# Save corresponding plot_data_all file

plot_data_log = as.data.frame(MFPCA_log$scores) %>% 
  mutate(scenario = c(rep("Control", 434),rep("SSP1-RCP2.6", 442),rep("SSP3-RCP7.0", 462),rep("SSP5-RCP8.5", 465)))
# plot_data_ref = as.data.frame(MFPCA_ref$scores) %>% 
#   mutate(scenario = c(rep("Control", 434),rep("SSP1-RCP2.6", 442),rep("SSP3-RCP7.0", 462),rep("SSP5-RCP8.5", 465)))

write.table(plot_data_log, "Scripts/Plots/Model/plot_data/plot_data_log.txt", row.names = TRUE, col.names = TRUE)
# write.table(plot_data_ref, "Scripts/Plots/Model/plot_data/plot_data_ref.txt", row.names = TRUE, col.names = TRUE)

  

