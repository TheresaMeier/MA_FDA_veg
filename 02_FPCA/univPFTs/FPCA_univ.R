################################################################################
############################ Master's Thesis ###################################
################################################################################

###################### Exploratory Analysis: FPCA ##############################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/02_FPCA/functions.R")
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

## Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)
scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")
pfts = c("Tundra", "BNE", "IBS", "otherC", "TeBS")

## Choose parameters:
start_year = 2015
end_year = 2040
#pft = "otherC"      # Choose from Tundra, BNE, IBS, otherC, TeBS
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

for (pft in pfts){
  print(paste("Start with PFT", long_names_pfts(tolower(pft))))
  for (scen in scenarios){
    ## Get data for all four scenarios in the appropriate shape
    d_scen = get_data_fpca(scen, start_year, end_year,pid,pft)
    
    ## Get basis representation
    fit.scen = get_basis_rep(start_year,end_year,d_scen[[1]][,-1])
    saveRDS(fit.scen, paste0("Scripts/MA_FDA_veg/02_FPCA/FdObjects/Wfdobj_", scen, "_", pft, ".rds"))
    
    # Transform the values to exp-scale for plotting
    fit.scen_2 = fit.scen
    fit.scen_2$Wfdobj$coefs = exp(fit.scen_2$Wfdobj$coefs)
    
    ## Run FPCA
    scen.pca = pca.fd(fit.scen_2$Wfdobj,2)
    
    # Save results
    assign(paste0("d_", scen, "_", pft), d_scen)
    assign(paste0("fit.", scen, "_", pft), fit.scen)
    assign(paste0("fit.", scen, "_", pft, "_2"), fit.scen_2)
    assign(paste0(scen, ".", pft, ".pca"), scen.pca)
    print(paste("PFT", long_names_pfts(tolower(pft)), ", scenario", long_names_scenarios(scen), "done."))
  }
}

############################ Plotting the results ##############################

for (pft in pfts){
  
  ## Plot fitted curves
  # pdf
  pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/BasisRep/pdf/BasisRep_",pft, "_",pid,".pdf"),width = 9, height = 9)
  par(mfrow=c(2,2))
  for (scen in scenarios){
    plot(x = c(1:100), y = rep(0,100), xlim = c(1,100), ylim = c(-0.05,1.05), type = 'l',lty = 2, cex.main = 1.1, xlab = "Year after Disturbance", ylab = "Share of aboveground carbon", main = paste0("Smoothed fit - ", long_names_scenarios(scen), " - ", long_names_pfts(tolower(pft))))
    fit.scen.pft = get(paste0("fit.", scen, "_", pft, "_2"))
    for (icurve in 1:nrow(fit.scen.pft$Wfdobj$coefs)){
      lines(fit.scen.pft$Wfdobj[icurve], col = pal(450)[icurve])
    }
  }
  dev.off()
  
  # png
  png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/BasisRep/png/BasisRep_",pft, "_",pid,".png"),width = 600, height = 600)
  par(mfrow=c(2,2))
  for (scen in scenarios){
    plot(x = c(1:100), y = rep(0,100), xlim = c(1,100), ylim = c(-0.05,1.05), type = 'l',lty = 2, cex.main = 1.1, xlab = "Year after Disturbance", ylab = "Share of aboveground carbon", main = paste0("Smoothed fit - ", long_names_scenarios(scen), " - ", long_names_pfts(tolower(pft))))
    fit.scen.pft = get(paste0("fit.", scen, "_", pft, "_2"))
    for (icurve in 1:nrow(fit.scen.pft$Wfdobj$coefs)){
      lines(fit.scen.pft$Wfdobj[icurve], col = pal(450)[icurve])
    }
  }
  dev.off()

  ## Plot PCA results
  # pdf
  pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/FPCA/unrotated/",pft,"/pdf/PCs_",pft, "_",pid,"_unrotated.pdf"), width = 10, height = 10)
  par(mfrow=c(2,2))
  for (scen in scenarios){
    plot.pca.fd(get(paste0(scen, ".", pft, ".pca")), xlab = long_names_scenarios(scen), xlim = c(1,100), cex = 1.1)
  }
  dev.off()

  # png
  png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/FPCA/unrotated/",pft,"/png/PCs_control_126_",pft, "_",pid,"_unrotated.png"),width = 1000, height = 1000)
  par(mfrow=c(2,2))
  plot.pca.fd(get(paste0("picontrol.", pft, ".pca")), xlab = "Control", xlim = c(1,100))
  plot.pca.fd(get(paste0("ssp126.", pft, ".pca")), xlab = "SSP1-RSP2.6", xlim = c(1,100))
  dev.off()

  png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/FPCA/unrotated/",pft,"/png/PCs_370_585_",pft, "_",pid,"_unrotated.png"),width = 1000, height = 1000)
  par(mfrow = c(2,2))
  plot.pca.fd(get(paste0("ssp370.", pft, ".pca")), xlab = "SSP3-RSP7.0", xlim = c(1,100))
  plot.pca.fd(get(paste0("ssp585.", pft, ".pca")), xlab = "SSP5-RSP8.5", xlim = c(1,100))
  dev.off()

  print(paste(pft, "done."))
}
