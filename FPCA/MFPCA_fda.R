################################################################################
############################ Master's Thesis ###################################
################################################################################

#################### Exploratory Analysis: MFPCA #############################

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
library(abind)

## Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)
scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")

## Choose parameters:
start_year = 2015
end_year = 2040
#pft = "TeBS"      # Choose from Tundra, BNE, IBS, otherC, TeBS
pid = 1           # Choose pid (int) or 'all'

run_fit = FALSE   # If true, the smoothing is done, otherwise, the stored fits are loaded
#######################

## Get data for all four scenarios in the appropriate shape
if (run_fit){
  for (scen in scenarios){
    d_Tundra = get_data_fpca(scen, start_year, end_year,pid,"Tundra")
    d_BNE = get_data_fpca(scen, start_year, end_year,pid,"BNE")
    d_IBS = get_data_fpca(scen, start_year, end_year,pid,"IBS")
    d_otherC = get_data_fpca(scen, start_year, end_year,pid,"otherC")
    d_TeBS = get_data_fpca(scen, start_year, end_year,pid,"TeBS")
    
    d_scen = abind(d_Tundra[[1]][,-1], d_BNE[[1]][,-1], d_IBS[[1]][,-1], d_otherC[[1]][,-1], d_TeBS[[1]][,-1], along = 3)
    
    # Get basis representation
    dif_years = end_year - start_year
    yearrange = seq(0,100+dif_years,by=1)
    yearbasis = create.bspline.basis(c(0,100+dif_years),ncol(d_scen), 6)
    
    WfdPar = fdPar(yearbasis, 3, 1) # 3 means it penalizes the third derivative
    
    yearrange.5 = c(0,yearrange[-1]-0.5)
    
    fit.scen = smooth.pos_TM(yearrange.5, d_scen, WfdPar)
    fit.scen$Wfdobj$fdnames = list('Year after Disturbance' = yearrange, 'Location/PID' = colnames(d_scen), 'Share of aboveground carbon' = c("Tundra", "BNE", "IBS", "otherC", "TeBS"))
    
    saveRDS(fit.scen, paste0("Scripts/FPCA/FdObjects/Wfdobj_", scen, ".rds"))
    assign(paste0("fit.",scen), fit.scen)
    fit.scen$Wfdobj$coef = exp(fit.scen$Wfdobj$coef)
    assign(paste0("fit.",scen,"_exp"), fit.scen)
  }
}

if(!run_fit) {
  fit.picontrol = readRDS("Scripts/FPCA/FdObjects/Wfdobj_picontrol.rds")
  fit.picontrol_exp = fit.picontrol
  fit.picontrol_exp$Wfdobj$coefs = exp(fit.picontrol_exp$Wfdobj$coefs)
  
  fit.ssp126 = readRDS("Scripts/FPCA/FdObjects/Wfdobj_ssp126.rds")
  fit.ssp126_exp = fit.picontrol
  fit.ssp126_exp$Wfdobj$coefs = exp(fit.ssp126_exp$Wfdobj$coefs)
  
  fit.ssp370 = readRDS("Scripts/FPCA/FdObjects/Wfdobj_ssp370.rds")
  fit.ssp370_exp = fit.picontrol
  fit.ssp370_exp$Wfdobj$coefs = exp(fit.ssp370_exp$Wfdobj$coefs)
  
  fit.ssp585 = readRDS("Scripts/FPCA/FdObjects/Wfdobj_ssp585.rds")
  fit.ssp585_exp = fit.picontrol
  fit.ssp585_exp$Wfdobj$coefs = exp(fit.ssp585_exp$Wfdobj$coefs)
  
}

for (scen in scenarios){
  
  fit.scen_exp = get(paste0("fit.", scen,"_exp"))
  
  # # Plot the fits 
  # pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/pdf/MFPCA_",scen,"_",pid,".pdf"),width = 18, height = 8)
  # par(mfrow = c(2,3))
  # plot.fd(fit.scen_exp$Wfdobj, xlab = "Year after Disturbance - Control")
  # dev.off()
  # 
  # png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/png/MFPCA_",scen,"_",pid,".png"),width = 1800, height = 800)
  # par(mfrow = c(2,3))
  # plot.fd(fit.scen_exp$Wfdobj, xlab = "Year after Disturbance - Control")
  # dev.off()
  # 
  # Run MFPCA
  scen.pca = pca.fd(fit.scen_exp$Wfdobj,2)
  assign(paste0(scen,".pca"), scen.pca)
  # pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/pdf/MFPCA_PC12_",scen,"_",pid,".pdf"),width = 18, height = 8)
  # par(mfrow = c(2,5))
  # plot.pca.fd(scen.pca, xlim = c(0,100), ylim = c(-0.5,1), xlab = "Year after Disturbance - Control")
  # dev.off()
  # 
  # png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/png/MFPCA_PC12_",scen,"_",pid,".png"),width = 1800, height = 800)
  # par(mfrow = c(2,5))
  # plot.pca.fd(scen.pca, xlim = c(0,100), ylim = c(-0.5,1), xlab = "Year after Disturbance - Control")
  # dev.off()
}


# Plot principal scores against each other

for (scen in scenarios){
  data_pca = get(paste0(scen,".pca"))
  i=1
  
  for (pft in c("Tundra", "BNE", "IBS", "otherC", "TeBS")) {
    
    plot_data_tmp = as.data.frame(data_pca$scores[,,i]) %>%
      rename(PC1 = V1,
             PC2 = V2) %>%
      mutate(PFT = long_names_pfts(tolower(pft)))
    i=i+1
    assign(paste0("plot_data_tmp_",pft), plot_data_tmp)
  }
  
  plot_data = purrr::reduce(list(plot_data_tmp_Tundra, plot_data_tmp_BNE, plot_data_tmp_IBS, plot_data_tmp_otherC, plot_data_tmp_TeBS), bind_rows)
  
  plot_data = plot_data %>%
    mutate(name = long_names_scenarios(scen))
  
  assign(paste0("plot_data_", scen), plot_data)
  
  
}

plot_data = purrr::reduce(list(plot_data_picontrol, plot_data_ssp126, plot_data_ssp370, plot_data_ssp585), bind_rows)

#write.table(plot_data, file = paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/plot_data/plot_data_",pft,"_",pid))

# Plot the scores

ggplot() + 
  geom_point(data = plot_data, 
             aes(x = PC1, y = PC2, color = PFT)) +
  facet_grid(rows = vars(name)) +
  scale_x_continuous(name = "PC 1") +
  scale_y_continuous(name = "PC 2") +
  ggtitle(paste("Principal Component Scores for one Patch")) +
  theme_bw() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))

ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/pdf/PC1_vs_PC2_",pid,".pdf"), width = 7, height = 4.5)
ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/png/PC1_vs_PC2_",pid,".png"), width = 7, height = 4.5)

