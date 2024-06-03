################################################################################
############################ Master's Thesis ###################################
################################################################################

############## Exploratory Analysis: MFPCA - for each scenario #################

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
library(abind)
library(MFPCA)

## Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)
scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")
pfts = c("Tundra", "BNE", "IBS", "otherC", "TeBS")
set.seed(1)

## Choose parameters:
start_year = 2015
end_year = 2040
#pft = "TeBS"      # Choose from Tundra, BNE, IBS, otherC, TeBS
pid = 1           # Choose pid (int) or 'all'
k = 4             # Number of clusters for kmeans algorithm
M = 10            # Number of PCs

# scen = "picontrol"
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
for (scen in scenarios){
  
  print(paste("Start with scenario", long_names_scenarios(scen)))
  print("Create funData object ...")
  # Create funData objects
  for (pft in pfts){
    print(pft)
    d_pft = get_data_fpca(scen, start_year, end_year,pid,pft)
    pft_funData = funData(argvals = 1:nrow(d_pft[[1]]), X = t(data.matrix(d_pft[[1]][,-1])))
    assign(paste0(pft, "_funData"), pft_funData)
    print("...done.")
  }
  print("... done creating funData object.")
  
  scen_funData = multiFunData(Tundra_funData, BNE_funData, IBS_funData, otherC_funData, TeBS_funData)
 
  # Create list of parameters for univariate FPCA
  uniExpansions <- list(list(type = "uFPCA", npc = 5), list(type = "uFPCA", npc = 5), list(type = "uFPCA", npc = 5), list(type = "uFPCA", npc = 5), list(type = "uFPCA", npc = 5))
  
  # Run MFPCA
  print(paste("Run MFPCA for scenario", long_names_scenarios(scen), "..."))
  MFPCA_scen = MFPCA(scen_funData, M=M, uniExpansions = uniExpansions, fit = TRUE)
  assign(paste0("MFPCA_", scen), MFPCA_scen)
  print("... done.")
  
  # Clustering of the coefficients
  kmeans_scen = kmeans(MFPCA_scen$scores, k)
  scores = cbind(MFPCA_scen$scores,kmeans_scen$cluster)
  
  names = c("PC1", "PC2")
  for (i in (3:M)) names = c(names, paste0("PC", i))
  
  colnames(scores) = c(names, "Cluster")
  assign(paste0("scores_", scen), scores)
  
  print(paste("Scenario", long_names_scenarios(scen), "done."))
}

######################### Plot the scores: PC1 vs. PC2 #########################
for (scen in scenarios){
  plot_data = as.data.frame(get(paste0("scores_", scen)))
  plot_data = plot_data %>% mutate(Cluster = as.factor(Cluster),
                                   name = long_names_scenarios(scen))
  assign(paste0("plot_data_", scen), plot_data)
}

plot_data = purrr::reduce(list(plot_data_picontrol, plot_data_ssp126, plot_data_ssp370, plot_data_ssp585), bind_rows)

ggplot(plot_data) + 
  geom_point(data = plot_data, 
             aes(x = PC1, y = PC2, col = Cluster)) +
  facet_grid(rows = vars(name)) +
  scale_x_continuous(name = "PC 1") + 
  scale_y_continuous(name = "PC 2") +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")) +
  ggtitle(paste("Principal Component Scores for one Patch - MFPCA")) +
  theme_bw() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))

ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca/pdf/PC1vsPC2/PC1_vs_PC2_",pid,"_clustered_", M, "PCs.pdf"), width = 7, height = 4.5)
ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca/png/PC1vsPC2/PC1_vs_PC2_",pid,"_clustered_", M, "PCs.png"), width = 7, height = 4.5)


######################### Plot principal components ############################
for (scen in scenarios){
  MFPCA_scen = get(paste0("MFPCA_", scen))
  
  for (iPC in (1:M)){
    pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca/pdf/PCs/", M, "PCs/PCs_",scen, "_",pid,"_PC", iPC, "_", M, "PCs.pdf"), width = 20, height = 5)
    par(mfrow = c(1,5))
    plot(MFPCA_scen, combined = TRUE, plotPCs = iPC, xlab = paste("Year after Disturbance -", long_names_scenarios(scen)), xlim = c(0,100), ylab = "Share of above ground carbon")
    #plot(MFPCA_scen, combined = TRUE, plotPCs = iPC, xlab = paste("Year after Disturbance -", long_names_scenarios(scen)), xlim = c(0,100), ylab = "Share of above ground carbon")
    dev.off()
    
    png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca/png/PCs/",M, "PCs/PCs_",scen, "_",pid,"_PC", iPC, "_", M, "PCs.png"), width = 2000, height = 500)
    par(mfrow = c(2,5))
    plot(MFPCA_scen, combined = TRUE, plotPCs = iPC, xlab = paste("Year after Disturbance -", long_names_scenarios(scen)), xlim = c(0,100), ylab = "Share of above ground carbon")
    plot(MFPCA_scen, combined = TRUE, plotPCs = iPC, xlab = paste("Year after Disturbance -", long_names_scenarios(scen)), xlim = c(0,100), ylab = "Share of above ground carbon")
    dev.off()
  }
}

######################## Plot reconstructed curves #############################
for (scen in scenarios){
  
  MFPCA_scen = get(paste0("MFPCA_", scen))
  i=1
  
  for (pft in pfts){

    p_rec = autoplot(MFPCA_scen$fit)[[i]] + 
      geom_line(aes(colour = obs), show.legend = FALSE) + ggtitle(paste("Reconstructed fit -", long_names_scenarios(scen), "-", long_names_pfts(tolower(pft)))) +
      theme_bw()  +
      scale_x_continuous(name = "Year after Disturbance", limits = c(0,100)) + 
      scale_y_continuous(name = "Share of above ground carbon", limits = c(-0.01,1.3)) +
      scale_color_manual(values = pal(470)) +
      theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))
    
    
    # Get original fit
    fit.scen.pft = readRDS(paste0("Scripts/FPCA/FdObjects/Wfdobj_", scen, "_", pft, ".rds"))
    fit.scen.pft_exp = fit.scen.pft
    fit.scen.pft_exp$Wfdobj$coefs = exp(fit.scen.pft_exp$Wfdobj$coefs)
    
    fit.funData = fd2funData(fit.scen.pft_exp$Wfdobj, argvals = fit.scen.pft_exp$argvals)
    
    p_orig = autoplot(fit.funData) + 
      geom_line(aes(colour = obs), show.legend = FALSE) + ggtitle(paste("Original fit -", long_names_scenarios(scen), "-", long_names_pfts(tolower(pft)))) +
      theme_bw()  +
      scale_x_continuous(name = "Year after Disturbance", limits = c(0,100)) + 
      scale_y_continuous(name = "Share of above ground carbon", limits = c(-0.01,1.3)) +
      scale_color_manual(values = pal(470)) +
      theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))
    
    
    plot_grid(p_orig, p_rec)
    ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca/pdf/reconstruction/", M, "PCs/PCA_reconstruct_", scen, "_", pft, "_", pid,"_", M, "PCs.pdf"), width = 10, height = 4.5)
    ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca/png/reconstruction/", M, "PCs/PCA_reconstruct_", scen, "_", pft, "_", pid,"_", M, "PCs.png"), width = 10, height = 4.5)
    
    i=i+1
    
  }
}

######################### Plot Cluster Assignments #############################
## Import fd objects
for (scen in scenarios){
  for (pft in pfts){
    fit.scen_pft = readRDS(paste0("Scripts/FPCA/FdObjects/Wfdobj_", scen, "_", pft, ".rds"))
    fit.scen_pft_exp = fit.scen_pft
    fit.scen_pft_exp$Wfdobj$coefs = exp(fit.scen_pft_exp$Wfdobj$coefs)
    
    assign(paste0("fit.", scen, "_", pft, "_exp"), fit.scen_pft_exp)
    
    locs = fit.scen_pft_exp$Wfdobj$fdnames$`Location/PID`
    
    scores = get(paste0("scores_", scen))
    # Resort scores to fit the Wfdobjects
    scores = scores[locs,]
    
    
    cluster1 = scores[, M+1] == 1
    assign(paste0("cluster1_",scen,"_",pft), cluster1)
    cluster2 = scores[, M+1] == 2
    assign(paste0("cluster2_",scen,"_",pft), cluster2)
    cluster3 = scores[, M+1] == 3
    assign(paste0("cluster3_",scen,"_",pft), cluster3)
    cluster4 = scores[, M+1] == 4
    assign(paste0("cluster4_",scen,"_",pft), cluster4)
  }
  

  
  print("Start plotting...")
  for (pft in pfts){
    
    # Save as pdf
    pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca/pdf/clusters/", M, "PCs/Cluster_", k,"-means_",scen, "_",pft, "_",pid,"_", M, "PCs.pdf"),width = 10, height = 10)
    par(mfrow = c(2,2))
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 1 -", sum(cluster1), "curves"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster1_", scen,"_", pft))], col = "#F8766D", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster1_", scen,"_", pft))]), col = "darkred", lwd = 4, lty = 1)
    
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 2 -", sum(cluster2), "curves"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster2_", scen,"_", pft))], col = "#7CAE00", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster2_", scen,"_", pft))]), col = "darkgreen", lwd = 4, lty = 1)
    
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 3 -", sum(cluster3), "curves"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster3_", scen,"_", pft))], col = "#00BFC4", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster3_", scen,"_", pft))]), col = "darkblue", lwd = 4, lty = 1)
    
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey",  main = paste("Cluster 4 -", sum(cluster4), "curves"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster4_", scen,"_", pft))], col = "#C77CFF", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster4_", scen,"_", pft))]), col = "purple4", lwd = 4, lty = 1)
    
    mtext(paste(long_names_scenarios(scen), "-", long_names_pfts(tolower(pft))), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
    dev.off()
    
    # Save as png
    png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca/png/clusters/", M, "PCs/Cluster_", k,"-means_",scen, "_",pft, "_",pid,"_", M, "PCs.png"),width = 1000, height = 1000)
    par(mfrow = c(2,2))
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 1 -", sum(cluster1), "curves"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster1_", scen,"_", pft))], col = "#F8766D", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster1_", scen,"_", pft))]), col = "darkred", lwd = 4, lty = 1)
    
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 2 -", sum(cluster2), "curves"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster2_", scen,"_", pft))], col = "#7CAE00", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster2_", scen,"_", pft))]), col = "darkgreen", lwd = 4, lty = 1)
    
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 3 -", sum(cluster3), "curves"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster3_", scen,"_", pft))], col = "#00BFC4", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster3_", scen,"_", pft))]), col = "darkblue", lwd = 4, lty = 1)
    
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey",  main = paste("Cluster 4 -", sum(cluster4), "curves"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster4_", scen,"_", pft))], col = "#C77CFF", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster4_", scen,"_", pft))]), col = "purple4", lwd = 4, lty = 1)
    
    mtext(paste(long_names_scenarios(scen), "-", long_names_pfts(tolower(pft))), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
    dev.off()
    
    print(paste("Plot", scen, pft, "done."))
  }
}

