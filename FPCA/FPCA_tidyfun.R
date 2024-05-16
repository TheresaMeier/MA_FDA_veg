################################################################################
############################ Master's Thesis ###################################
################################################################################

############### Exploratory Analysis: visualization with tidyfun ###############

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
library("tidyfun")

## Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)
scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")

## Choose parameters:
start_year = 2015
end_year = 2040
pft = "Tundra"      # Choose from Tundra, BNE, IBS, otherC, TeBS
pid = 1           # Choose pid (int) or 'all'

#############################


## Get data for all four scenarios in the appropriate shape
for (scen in scenarios){
  print(paste("Start with scenario", long_names_scenarios(scen)))
  d_Tundra = get_data_fpca(scen, start_year, end_year,pid,"Tundra")
  d_BNE = get_data_fpca(scen, start_year, end_year,pid,"BNE")
  d_IBS = get_data_fpca(scen, start_year, end_year,pid,"IBS")
  d_otherC = get_data_fpca(scen, start_year, end_year,pid,"otherC")
  d_TeBS = get_data_fpca(scen, start_year, end_year,pid,"TeBS")
  
  print("Data successfully loaded.")
  
  data_loc = d_Tundra[[2]] %>% 
    ungroup() %>%
    mutate(
      Tundra = t(d_Tundra[[1]][,-1]) %>%
        tfb(arg = 0:125, basis = "spline", k = 10, m = 4, bounds = c(0,1)),
      BNE = t(d_BNE[[1]][,-1]) %>%
        tfb(arg = 0:125, basis = "spline", k = 10, m = 4),
      IBS = t(d_IBS[[1]][,-1]) %>%
        tfb(arg = 0:125, basis = "spline", k = 10, m = 4),
      otherC = t(d_otherC[[1]][,-1]) %>%
        tfb(arg = 0:125, basis = "spline", k = 10, m = 4),
      TeBS = t(d_TeBS[[1]][,-1]) %>%
        tfb(arg = 0:125, basis = "spline", k = 10, m = 4),
      Region = classify_region(Lat,Lon)
    )
  
  Tundra_panel <- data_loc |>
    ggplot(aes(y = Tundra, color = Region)) +
    geom_spaghetti() + theme_bw()
  
  BNE_panel <- data_loc |>
    ggplot(aes(y = BNE, color = Region)) +
    geom_spaghetti() + theme_bw()
  
  IBS_panel <- data_loc |>
    ggplot(aes(y = IBS, color = Region)) +
    geom_spaghetti() + theme_bw() #+ scale_y_manual(labels = long_names_scenarios("bne")) 
  
  otherC_panel <- data_loc |>
    ggplot(aes(y = otherC, color = Region)) +
    geom_spaghetti() + theme_bw() #+ scale_y_manual(labels = long_names_scenarios("bne")) 
  
  TeBS_panel <- data_loc |>
    ggplot(aes(y = TeBS, color = Region)) +
    geom_spaghetti() + theme_bw() #+ scale_y_manual(labels = long_names_scenarios("bne")) 
  
  pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/tidyfun/fit_", scen, ".pdf"), width = 25, height = 5)
  gridExtra::grid.arrange(Tundra_panel, BNE_panel, IBS_panel, otherC_panel, TeBS_panel, nrow = 1, top = long_names_scenarios(scen))
  dev.off()
  
  png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/tidyfun/fit_", scen, ".png"), width = 2500, height = 500)
  gridExtra::grid.arrange(Tundra_panel, BNE_panel, IBS_panel, otherC_panel, TeBS_panel, nrow = 1, top = long_names_scenarios(scen))
  dev.off()
  print("Plots saved.")
}

