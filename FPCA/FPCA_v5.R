################################################################################
############################ Master's Thesis ###################################
################################################################################

######### Exploratory Analysis: FPCA - Reconstruction of the curves ############

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/FPCA/functions.R")
source("Scripts/Description/utils.R")

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

## Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)
scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")
pfts = c("Tundra", "BNE", "IBS", "otherC", "TeBS")

## Choose parameters:
start_year = 2015
end_year = 2040
#pft = "Tundra"      # Choose from Tundra, BNE, IBS, otherC, TeBS
pid = 1           # Choose pid (int) or 'all'

#######################
## Get 450 different colors for plotting

palettes <- c("Set1", "Set2", "Set3", "Dark2", "Paired", "Accent")
all_colors <- character()
for (palette_name in palettes) {
  palette <- brewer.pal(20, palette_name)
  all_colors <- c(all_colors, palette)
}
palette_450 <- unique(all_colors[1:450])
pal <- colorRampPalette(palette_450)

#######################
## Import fd objects

for (scen in scenarios){
  
  for (pft in pfts){
    
    ## Import fd object
    fit.scen.pft = readRDS(paste0("Scripts/FPCA/FdObjects/Wfdobj_", scen, "_", pft, ".rds"))
    fit.scen.pft_exp = fit.scen.pft
    fit.scen.pft_exp$Wfdobj$coefs = exp(fit.scen.pft_exp$Wfdobj$coefs)
    
    ## Run FPCA
    
    fit.pca = pca.fd(fit.scen.pft_exp$Wfdobj,3)
    
    ## Plot reconstructed and original fits for using 1, 2 and 3 PCs
    for (nPC in 1:3){
      
      # Save as pdf
      pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/unrotated/",pft,"/PCA_reconstruct_",scen, "_",pft, "_",pid,"_",nPC,"PCs.pdf"), width = 10, height = 4.5)
      par(mfrow = c(1,2))
      plot(x = c(1:100),y = rep(0,100), xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',lty = 2, xlab = "Year after Disturbance", ylab = "Share of above ground carbon", main = paste0("Original fit - ", long_names_scenarios(scen), " - ", long_names_pfts(tolower(pft))))
      for (icurve in 1:nrow(fit.pca$scores)){
        lines(fit.scen.pft_exp$Wfdobj[icurve], col = pal(450)[icurve])
      }
      plot(x = c(1:100),y = rep(0,100), xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',lty = 2, xlab = "Year after Disturbance", ylab = "Share of above ground carbon", main = paste("Reconstructed fit using", nPC, "PCs"))
      for (icurve in 1:nrow(fit.pca$scores)){
        if (nPC == 1) lines(fit.pca$harmonics[1,] *fit.pca$scores[icurve,1] + fit.pca$meanfd, col = pal(450)[icurve])
        if (nPC == 2) lines(fit.pca$harmonics[1,] *fit.pca$scores[icurve,1] + fit.pca$harmonics[2,]*fit.pca$scores[icurve,2] + fit.pca$meanfd, col = pal(450)[icurve])
        if (nPC == 3) lines(fit.pca$harmonics[1,] *fit.pca$scores[icurve,1] + fit.pca$harmonics[2,]*fit.pca$scores[icurve,2] + fit.pca$harmonics[3,]*fit.pca$scores[icurve,3] + fit.pca$meanfd, col = pal(450)[icurve])
      }
      dev.off()
      
      # Save as png
      png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/unrotated/",pft,"/PCA_reconstruct_",scen, "_",pft, "_",pid,"_",nPC,"PCs.png"), width = 1000, height = 450)
      par(mfrow = c(1,2))
      plot(x = c(1:100),y = rep(0,100), xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',lty = 2, xlab = "Year after Disturbance", ylab = "Share of above ground carbon", main = paste0("Original fit - ", long_names_scenarios(scen), " - ", long_names_pfts(tolower(pft))))
      for (icurve in 1:nrow(fit.pca$scores)){
        lines(fit.scen.pft_exp$Wfdobj[icurve], col = pal(450)[icurve])
      }
      plot(x = c(1:100),y = rep(0,100), xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',lty = 2, xlab = "Year after Disturbance", ylab = "Share of above ground carbon", main = paste("Reconstructed fit using", nPC, "PCs"))
      for (icurve in 1:nrow(fit.pca$scores)){
        if (nPC == 1) lines(fit.pca$harmonics[1,] *fit.pca$scores[icurve,1] + fit.pca$meanfd, col = pal(450)[icurve])
        if (nPC == 2) lines(fit.pca$harmonics[1,] *fit.pca$scores[icurve,1] + fit.pca$harmonics[2,]*fit.pca$scores[icurve,2] + fit.pca$meanfd, col = pal(450)[icurve])
        if (nPC == 3) lines(fit.pca$harmonics[1,] *fit.pca$scores[icurve,1] + fit.pca$harmonics[2,]*fit.pca$scores[icurve,2] + fit.pca$harmonics[3,]*fit.pca$scores[icurve,3] + fit.pca$meanfd, col = pal(450)[icurve])
      }
      dev.off()
    }
  }
}
