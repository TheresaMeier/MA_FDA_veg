################################################################################
############################ Master's Thesis ###################################
################################################################################

############## Exploratory Analysis: MFPCA - for all scenarios #################

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
      d_scen = get_data_fpca(scen, start_year, end_year,pid,pft)
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
    saveRDS(fun_pft, paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_", pft, ".rds"))
  }
  
  # Create multiFunData object
  multiFun_pft = multiFunData(fun_BNE, fun_IBS, fun_otherC, fun_TeBS, fun_Tundra)
  
  saveRDS(multiFun_pft, "Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_1803.rds")
  
  print("... all done.")
}

# Run MFPCA
if (runMFPCA){
  multiFun_pft = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_1803.rds")
  
  uniExpansions <- list(list(type = "uFPCA", npc = 10),
                        list(type = "uFPCA", npc = 10),
                        list(type = "uFPCA", npc = 10),
                        list(type = "uFPCA", npc = 10),
                        list(type = "uFPCA", npc = 10))
  
  MFPCA_all <- MFPCA_2(multiFun_pft, M = M, fit = TRUE, uniExpansions = uniExpansions)
  saveRDS(MFPCA_all, paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_all_1803_100y.rds"))
}

# Get multiFunData objects
multiFun_pft = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_1803.rds")

# Get MFPCA results
MFPCA_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_all_1803_100y.rds")

set.seed(1)
# Compute WCSS for different values of k
wcss <- sapply(1:10, function(k) {
  kmeans(MFPCA_all$scores[,1:M], centers = k)$tot.withinss
})

# Save elbow plot as pdf
pdf(paste0("Scripts/Plots/MFPCA/Clusters/pdf/elbow_plot_", pid, ".pdf"), width = 12, height = 8)
par(mfrow = c(1,1))
# Plot the elbow curve
plot(1:10, wcss, type = "b", pch = 19, frame = TRUE, 
     xlab = "Number of Clusters", ylab = "Within-cluster Sum of Squares",
     main = "Elbow Plot for k-means Clustering",
     cex.axis = 1.3,   # Adjust the font size of axis labels
     cex.lab = 1.3,    # Adjust the font size of axis titles
     cex.main = 1.5)   # Adjust the font size of main title

# Add lines for visual aid
abline(v = 1:10, lty = 3, col = "gray")
dev.off()

# Save it as png
png(paste0("Scripts/Plots/MFPCA/Clusters/png/elbow_plot_", pid, ".png"), width = 1200, height = 800)
par(mfrow = c(1, 1))

# Plot the elbow curve
plot(1:10, wcss, type = "b", pch = 19, frame = TRUE, 
     xlab = "Number of Clusters", ylab = "Within-cluster Sum of Squares",
     main = "Elbow Plot for k-means Clustering",
     cex.axis = 1.3,   # Adjust the font size of axis labels
     cex.lab = 1.3,    # Adjust the font size of axis titles
     cex.main = 1.5)   # Adjust the font size of main title

# Add lines for visual aid
abline(v = 1:10, lty = 3, col = "gray")
dev.off()

# Clustering of the coefficient
set.seed(1)
kmeans_result = kmeans(MFPCA_all$scores[,1:M],4)
scores = cbind(MFPCA_all$scores[,1:M], kmeans_result$cluster)

names = c("PC1", "PC2")
for (i in (3:M)) names = c(names, paste0("PC", i))

colnames(scores) = c(names, "Cluster")
plot_data_all = as.data.frame(scores) %>% 
  mutate(Cluster = as.factor(Cluster),
         # scenario = c(rep("Control", dim(d_picontrol)[1]),rep("SSP1-RCP2.6", dim(d_ssp126)[1]),rep("SSP3-RCP7.0", dim(d_ssp370)[1]),rep("SSP5-RCP8.5", dim(d_ssp585)[1])))
         scenario = c(rep("Control", 434),rep("SSP1-RCP2.6", 442),rep("SSP3-RCP7.0", 462),rep("SSP5-RCP8.5", 465)))

# right_order_all = c(4,2,1,3)
# plot_data_all$Cluster = right_order_all[plot_data_all$Cluster]

write.table(plot_data_all, "Scripts/Plots/MFPCA/plot_data/plot_data_all.txt", row.names = TRUE, col.names = TRUE)

################## Cluster description for the whole time span ##################

list_var = c()
list_pft = c()
list_cl = c()
for (iPFT in c(1:5)){
  d_pft = multiFun_pft@.Data[[iPFT]]@X
    
  for (iCl in c(1:4)) {
    list_var = c(list_var, colMeans(d_pft[plot_data_all$Cluster == iCl,]))
    list_cl = c(list_cl, rep(iCl, 100))
  }
  
  list_pft = c(list_pft, rep(pfts[iPFT],100*4))
}

d_plot = as.data.frame(cbind("year" = rep(c(1:100),20), "value" = list_var, "Cluster" = list_cl, "PFT" = long_names_pfts(tolower(list_pft)))) %>%
  mutate(year = as.numeric(year),
         value = as.numeric(value),
         Cluster = as.factor(Cluster),
         PFT = as.factor(PFT))

# Plot the mean curves per cluster and per PFT

ggplot(d_plot) + 
  geom_line(data = d_plot, aes(x = year, y = value, col = Cluster), lwd = 2) +
  facet_grid(rows = vars(PFT)) + 
  scale_y_continuous(name = "Share of aboveground carbon") + 
  scale_x_continuous(name = "Year after Disturbance", breaks = seq(0,100,by=10), limits = c(0,100)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  theme_bw() + theme(text = element_text(size = 15),plot.title = element_text(size = 20, face = "bold",hjust = 0.5)) +
  ggtitle("Cluster-wise mean share of aboveground carbon over time")

ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/cluster_means_per_cluster_",pid,".pdf"), width = 15, height = 10)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/cluster_means_per_cluster_",pid,".png"), width = 15, height = 10)



d_plot = d_plot %>%
  mutate(Cluster = paste("Cluster", Cluster))

ggplot(d_plot) + 
  geom_line(data = d_plot, aes(x = year, y = value, col = PFT), lwd = 2) +
  facet_grid(rows = vars(Cluster)) + 
  scale_y_continuous(name = "Share of aboveground carbon") + 
  scale_x_continuous(name = "Year after Disturbance", breaks = seq(0,100,by=10), limits = c(0,100)) +
  scale_color_manual(name = "Dominant vegetation", drop = TRUE,
                     values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                                "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  theme_bw() + theme(text = element_text(size = 20), plot.title = element_text(size = 20, face = "bold",hjust = 0.5)) +
  ggtitle("PFT-wise mean share of aboveground carbon over time")


ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/cluster_means_per_PFT_",pid,".pdf"), width = 10, height = 6)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/cluster_means_per_PFT_",pid,".png"), width = 15, height = 10)



######################### Plot the scores: PC1 vs. PC2 #########################
plot_data_all = as.data.frame(scores) %>% 
  mutate(Cluster = as.factor(Cluster),
         # scenario = c(rep("Control", dim(d_picontrol)[1]),rep("SSP1-RCP2.6", dim(d_ssp126)[1]),rep("SSP3-RCP7.0", dim(d_ssp370)[1]),rep("SSP5-RCP8.5", dim(d_ssp585)[1])))
         scenario = c(rep("Control", 434),rep("SSP1-RCP2.6", 442),rep("SSP3-RCP7.0", 462),rep("SSP5-RCP8.5", 465)))

g1 = ggplot(plot_data_all) + 
  geom_point(data = plot_data_all, 
             aes(x = PC1, y = PC2, col = Cluster)) +
  scale_x_continuous(name = "PC 1") + 
  scale_y_continuous(name = "PC 2") +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle(paste("PC Scores derived by MFPCA")) +
  theme_bw() +
  theme(text = element_text(size = 15), plot.title = element_text(size = 20, face = "bold",hjust = 0.5))

ggsave(paste0("Scripts/Plots/MFPCA/PC1vsPC2/pdf/PC1_vs_PC2_",pid,"_clustered_all.pdf"), plot = g1, width = 7, height = 4.5)
ggsave(paste0("Scripts/Plots/MFPCA/PC1vsPC2/png/PC1_vs_PC2_",pid,"_clustered_all.png"), plot = g1, width = 7, height = 4.5)

g2 = ggplot(plot_data_all) + 
  geom_point(data = plot_data_all, 
             aes(x = PC1, y = PC2, col = scenario)) +
  scale_x_continuous(name = "PC 1") + 
  scale_y_continuous(name = "PC 2") +
  scale_color_manual(name = "Scenario", values = c("Control" = "#F8766D", "SSP1-RCP2.6" = "#7CAE00" , "SSP3-RCP7.0" = "#00BFC4", "SSP5-RCP8.5" = "#C77CFF", "5" = "darkgrey")) +
  ggtitle(paste("PC Scores derived by MFPCA")) +
  theme_bw() +
  theme(text = element_text(size = 15), plot.title = element_text(size = 20, face = "bold",hjust = 0.5))

ggsave(paste0("Scripts/Plots/MFPCA/PC1vsPC2/pdf/PC1_vs_PC2_",pid,"_scenarios.pdf"), plot = g2, width = 7, height = 4.5)
ggsave(paste0("Scripts/Plots/MFPCA/PC1vsPC2/png/PC1_vs_PC2_",pid,"_scenarios.png"), plot = g2, width = 7, height = 4.5)


plot_grid(g1,g2)
ggsave(paste0("Scripts/Plots/MFPCA/PC1vsPC2/pdf/PC1_vs_PC2_",pid,"_scenarios_clustered_all.pdf"), width = 14, height = 7)
ggsave(paste0("Scripts/Plots/MFPCA/PC1vsPC2/png/PC1_vs_PC2_",pid,"_scenarios_clustered_all.png"), width = 14, height = 7)

######################### Plot principal components ############################

# Plot 10 first PCs
for (iPC in (1:10)){
  pdf(paste0("Scripts/Plots/MFPCA/PCs/pdf/PCs_PC", iPC, "_Presentation.pdf"), width = 15, height = 5)
  plot.MFPCAfit_2(MFPCA_all, plotPCs = iPC, combined = TRUE, xlab = paste("Year after Disturbance"), 
                  xlim = c(1,100), ylab = "Share of aboveground carbon", cex.main = 2.3, cex.lab = 1.6,
                  cols = c("#0072B2", "#E69F00", "#56B4E9",  "#D55E00", "#009E73"),
                  main = long_names_pfts(tolower(pfts)))
  
  dev.off()
  
  # png(paste0("Scripts/Plots/MFPCA/PCs/png/PCs_PC", iPC, ".png"), width = 1500, height = 500)
  # plot(MFPCA_all, plotPCs = iPC, combined = TRUE, xlab = paste("Year after Disturbance"), xlim = c(1,100), ylab = "PFT", cex = 1.5)
  # #mtext(paste("PC", iPC, "for all four scenarios"), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
  # 
  # dev.off()
}

######################## Plot reconstructed curves #############################
bounds = t(matrix(c(1,434,435,876,877,1338,1339,1803), nrow = 2))
iPFT = 0

for(pft in pfts){
  iPFT = iPFT + 1
  iScen = 0
  for (scen in scenarios){
    iScen = iScen + 1
    nOb = bounds[iScen,2] - bounds[iScen,1] + 1
    # Save as pdf
    pdf(paste0("Scripts/Plots/MFPCA/Reconstruction/pdf/", scen, "/", pft,"/MFPCA_reconstruct_",scen, "_",pft, "_",M,"PCs.pdf"), width = 8, height = 6.5)
    plot(MFPCA_all$fit[[iPFT]], obs = c(bounds[iScen,1]: bounds[iScen,2]), col = pal(nOb)[1:nOb],
         xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',
         cex.main = 1.8, cex = 1.5, 
         xlab = "Year after Disturbance",
         ylab = "Share of aboveground carbon", 
         main = paste("Reconstructed fit using", M, "PCs") )
    lines(x = c(1:100),y = rep(0,100), xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',lty = 2)
    dev.off()
    
    # Save as png
    png(paste0("Scripts/Plots/MFPCA/Reconstruction/png/", scen, "/", pft,"/MFPCA_reconstruct_",scen, "_",pft, "_",M,"PCs.png"), width = 800, height = 650)
    plot(MFPCA_all$fit[[1]], obs = c(1:434), col = pal(450)[1:434],
         xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',
         cex.main = 1.8, cex = 1.5, 
         xlab = "Year after Disturbance",
         ylab = "Share of aboveground carbon", 
         main = paste("Reconstructed fit using", M, "PCs") )
    lines(x = c(1:100),y = rep(0,100), xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',lty = 2)
    dev.off()
  }
}

######################### Plot Cluster Assignments #############################
## Import fd objects
for (scen in scenarios){
  for (pft in pfts){
    fit.scen_pft = readRDS(paste0("Scripts/MA_FDA_veg/02_FPCA/FdObjects/Wfdobj_", scen, "_", pft, ".rds"))
    fit.scen_pft_exp = fit.scen_pft
    fit.scen_pft_exp$Wfdobj$coefs = exp(fit.scen_pft_exp$Wfdobj$coefs)
    
    assign(paste0("fit.", scen, "_", pft, "_exp"), fit.scen_pft_exp)
    
    locs = fit.scen_pft_exp$Wfdobj$fdnames$`Location/PID`
    
    # Resort scores to fit the Wfdobjects
    scores_pft = scores[locs,]
    
    
    cluster1 = scores_pft[, M+1] == 1
    assign(paste0("cluster1_",scen,"_",pft), cluster1)
    cluster2 = scores_pft[, M+1] == 2
    assign(paste0("cluster2_",scen,"_",pft), cluster2)
    cluster3 = scores_pft[, M+1] == 3
    assign(paste0("cluster3_",scen,"_",pft), cluster3)
    cluster4 = scores_pft[, M+1] == 4
    assign(paste0("cluster4_",scen,"_",pft), cluster4)
    cluster5 = scores_pft[, M+1] == 5
    assign(paste0("cluster5_",scen,"_",pft), cluster5)
  }
  
  
  print("Start plotting...")
  for (pft in pfts){
    
    # Save as pdf
    pdf(paste0("Scripts/Plots/MFPCA/Clusters/pdf/Cluster_", k,"-means_",scen, "_",pft, "_",pid, ".pdf"),width = 10, height = 10)
    par(mfrow = c(2,2))
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 1 -", sum(cluster1), "of", sum(scores[,M+1] == 1), "grid cells"), lty = 1, xlim = c(0,100))
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster1_", scen,"_", pft))], col = "#F8766D", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster1_", scen,"_", pft))]), col = "darkred", lwd = 4, lty = 1)
    
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 2 -", sum(cluster2), "of", sum(scores[,M+1] == 2), "grid cells"), lty = 1, xlim = c(0,100))
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster2_", scen,"_", pft))], col = "#7CAE00", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster2_", scen,"_", pft))]), col = "darkgreen", lwd = 4, lty = 1)
    
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 3 -", sum(cluster3), "of", sum(scores[,M+1] == 3), "grid cells"), lty = 1, xlim = c(0,100))
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster3_", scen,"_", pft))], col = "#00BFC4", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster3_", scen,"_", pft))]), col = "darkblue", lwd = 4, lty = 1)
    
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey",  main = paste("Cluster 4 -", sum(cluster4), "of", sum(scores[,M+1] == 4), "grid cells"), lty = 1, xlim = c(0,100))
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster4_", scen,"_", pft))], col = "#C77CFF", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster4_", scen,"_", pft))]), col = "purple4", lwd = 4, lty = 1)
    
    mtext(paste(long_names_scenarios(scen), "-", long_names_pfts(tolower(pft)), "-", length(get(paste0("cluster1_", scen, "_", pft))), "grid cells in total"), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
    dev.off()
    
    # Save as png
    png(paste0("Scripts/Plots/MFPCA/Clusters/png/Cluster_", k,"-means_",scen, "_",pft, "_",pid, ".png"),width = 1000, height = 1000)
    par(mfrow = c(2,2))
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 1 -", sum(cluster1), "of", sum(scores[,M+1] == 1), "grid cells"), lty = 1, xlim = c(0,100))
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster1_", scen,"_", pft))], col = "#F8766D", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster1_", scen,"_", pft))]), col = "darkred", lwd = 4, lty = 1)
    
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 2 -", sum(cluster2), "of", sum(scores[,M+1] == 2), "grid cells"), lty = 1, xlim = c(0,100))
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster2_", scen,"_", pft))], col = "#7CAE00", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster2_", scen,"_", pft))]), col = "darkgreen", lwd = 4, lty = 1)
    
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 3 -", sum(cluster3), "of", sum(scores[,M+1] == 3), "grid cells"), lty = 1, xlim = c(0,100))
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster3_", scen,"_", pft))], col = "#00BFC4", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster3_", scen,"_", pft))]), col = "darkblue", lwd = 4, lty = 1)
    
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey",  main = paste("Cluster 4 -", sum(cluster4), "of", sum(scores[,M+1] == 4), "grid cells"), lty = 1, xlim = c(0,100))
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster4_", scen,"_", pft))], col = "#C77CFF", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster4_", scen,"_", pft))]), col = "purple4", lwd = 4, lty = 1)
    
    mtext(paste(long_names_scenarios(scen), "-", long_names_pfts(tolower(pft)), "-", length(get(paste0("cluster1_", scen, "_", pft))), "grid cells in total"), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
    dev.off()
    
    print(paste("Plot", scen, pft, "done."))
  }
}
############################ Screeplot of MFPCA ################################

pdf(paste0("Scripts/Plots/MFPCA/PCs/pdf/screeplot_",pid, ".pdf"),width = 10, height = 8)
screeplot(MFPCA_all, main = "Screeplot of MFPCA", cex = 1.3)
dev.off()

png(paste0("Scripts/Plots/MFPCA/PCs/png/screeplot_",pid, ".png"),width = 1000, height = 800)
screeplot(MFPCA_all, main = "Screeplot of MFPCA", cex = 1.3)
dev.off()

#################### Means per Cluster with original curves ####################
d_picontrol = get_data_scenario("picontrol", 2015, 2040, "cmass", 1)
d_ssp126 = get_data_scenario("ssp126", 2015, 2040, "cmass", 1)
d_ssp370 = get_data_scenario("ssp370", 2015, 2040, "cmass", 1)
d_ssp585 = get_data_scenario("ssp585", 2015, 2040, "cmass", 1)

df_trajectories = purrr::reduce(list(d_picontrol, d_ssp126, d_ssp370, d_ssp585), bind_rows) %>%
  mutate(PFT = long_names_pfts(tolower(PFT))) # make pft names pretty

df_trajectories = df_trajectories %>%
  filter(age > 0 & age < 101)

locs_disturbed$name = locs_disturbed$Scenario

merged_df <- merge(df_trajectories, locs_disturbed[,c(1,2,5,6)], by = c('Lon', 'Lat', 'name'), all.x = TRUE)

# Step 2: Create the Cluster column in df_large based on the merged data
df_large$Cluster <- merged_df$Cluster

df_mean = merged_df %>%
  group_by(age, PFT, Cluster) %>%
  summarize(relative = mean(relative, na.rm = T)) 

merged_df = merged_df %>%
  mutate(Cluster = if_else(Cluster == 1,"Cluster 1",
                           if_else(Cluster == 2, "Cluster 2", 
                                   if_else(Cluster == 3, "Cluster 3",
                                           if_else(Cluster == 4, "Cluster 4", NA)))))

# Assuming 'd_plot' represents the data frame containing the lines you want to highlight

ggplot(merged_df) + 
  geom_line(data = merged_df, linewidth = .05, alpha = .15,
            aes(x = age, y = relative, color = PFT, group = interaction(Lon, Lat, PID, PFT))) +
  geom_path(data = d_plot, aes(x = year, y = value, color = PFT, group = PFT), size = 1.5, color = "black") +  # Thick black outline
  geom_line(data = d_plot, aes(x = year, y = value, color = PFT, group = PFT), size = 1) +  # Thinner lines on top
  facet_grid(rows = vars(Cluster)) +
  scale_x_continuous(name = "Year after disturbance", expand = c(0, 0), limits = c(0, 100)) +
  scale_y_continuous(name = "Share of aboveground carbon", expand = c(0, 0),
                     breaks = c(0.25, 0.50, 0.75, 1.00)) +
  scale_color_manual(name = "Dominant vegetation", drop = TRUE,
                     values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00", "Needleleaf evergreen" = "#0072B2",   
                                "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  theme_bw() + 
  theme(
    text = element_text(size = 15), 
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  ggtitle("PFT-wise mean share of aboveground carbon over time")

ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/cluster_means_per_PFT_all.pdf"), width = 15, height = 10)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/cluster_means_per_PFT_all.png"), width = 15, height = 10)

