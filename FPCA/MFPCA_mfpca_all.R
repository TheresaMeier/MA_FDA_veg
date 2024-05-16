################################################################################
############################ Master's Thesis ###################################
################################################################################

########## Exploratory Analysis: MFPCA - for all scenarios together ############

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/FPCA/functions.R")
source("Scripts/Description/utils.R")
# source("Scripts/FPCA/MFPCA_calculation.R")

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
library(irlba)

## Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)
scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")
pfts = c("Tundra", "BNE", "IBS", "otherC", "TeBS")
set.seed(1) #4,8

## Choose parameters:
start_year = 2015
end_year = 2040
#pft = "TeBS"      # Choose from Tundra, BNE, IBS, otherC, TeBS
pid = 1           # Choose pid (int) or 'all'
k = 4             # Number of clusters for kmeans algorithm
M = 10            # Number of PCs
createFunData = FALSE
runMFPCA = FALSE
makeZlim = TRUE

# scen = "pTRUE# scen = "picontrol"
######################
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
if (createFunData)
{
  for (scen in scenarios){
    print(paste("Start with scenario", long_names_scenarios(scen)))
    print("Create funData object ...")
    # Create funData objects
    for (pft in pfts){
        print(pft)
        d_pft = get_data_fpca(scen, start_year, end_year,pid,pft)
        assign(paste0("d_",pft), d_pft)
        print("...done.")
    }
    d_scen = abind(t(d_Tundra[[1]][,-1]), t(d_BNE[[1]][,-1]), t(d_IBS[[1]][,-1]), t(d_otherC[[1]][,-1]), t(d_TeBS[[1]][,-1]), along = 3)
    assign(paste0("d_", scen), d_scen)
  }
  
  all_grid_cells <- unique(c(rownames(d_picontrol), rownames(d_ssp126), rownames(d_ssp370), rownames(d_ssp585)))
  
  # Pad the smaller matrices for function multiFunData
  
  for (scen in scenarios){
    d_scen = get(paste0("d_", scen))
  
    padded_matrix = array(0, dim = c(length(all_grid_cells), dim(d_scen)[2],5))
    padded_matrix[1:dim(d_scen)[1],,] = d_scen
    rownames(padded_matrix) = c(rownames(d_scen), setdiff(all_grid_cells, rownames(d_scen)))
  
    # Sort equally
    padded_matrix = padded_matrix[all_grid_cells,,]
    assign(paste0("d_", scen, "_padded"), padded_matrix)
  
    scen_funData = funData(argvals = list(1:dim(padded_matrix)[2], 1:5), X = padded_matrix)
    assign(paste0(scen, "_funData"), scen_funData)
    print("... done creating funData object.")
    saveRDS(scen_funData, paste0("Scripts/FPCA/FdObjects/funData_", scen, ".rds"))
    
  }
  
  funData_all = multiFunData(picontrol_funData, ssp126_funData, ssp370_funData, ssp585_funData)
  saveRDS(funData_all, "Scripts/FPCA/FdObjects/funData_all.rds")
}

# Get funData objects
picontrol_funData = readRDS("Scripts/FPCA/FdObjects/funData_picontrol.rds")
ssp126_funData = readRDS("Scripts/FPCA/FdObjects/funData_ssp126.rds")
ssp370_funData = readRDS("Scripts/FPCA/FdObjects/funData_ssp370.rds")
ssp585_funData = readRDS("Scripts/FPCA/FdObjects/funData_ssp585.rds")
funData_all = readRDS("Scripts/FPCA/FdObjects/funData_all.rds")

# Create list of parameters for univariate FPCA
uniExpansions <- list(list(type = "given", picontrol_funData, ortho = TRUE), list(type = "given", ssp126_funData, ortho = TRUE), list(type = "given", ssp370_funData, ortho = TRUE),list(type = "given", ssp585_funData, ortho = TRUE))

# Run MFPCA
if (runMFPCA){
  print("Run MFPCA...")
  MFPCA_all = MFPCA(funData_all, M=M, uniExpansions = uniExpansions, fit = FALSE, approx.eigen = TRUE)
  saveRDS(MFPCA_all, "Scripts/FPCA/FdObjects/MFPCA_all.rds")
  
  print("... done.")
}

MFPCA_all = readRDS("Scripts/FPCA/FdObjects/MFPCA_all.rds")

# Compute WCSS for different values of k
wcss <- sapply(1:10, function(k) {
  kmeans(MFPCA_all$scores, centers = k)$tot.withinss
})

# Save elbow plot as pdf
pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca_all/pdf/clusters/elbow_plot_", pid, ".pdf"), width = 12, height = 8)
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
png(paste0("Scripts/Plots/FPCA/PCs_", start_year, "_", end_year, "/MFPCA/mfpca_all/png/clusters/elbow_plot_", pid, ".png"), width = 1200, height = 800)
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

kmeans_all <- kmeans(MFPCA_all$scores, centers = 4)
scores = as.matrix(cbind(MFPCA_all$scores,kmeans_all$cluster))

names = c("PC1", "PC2")
for (i in (3:M)) names = c(names, paste0("PC", i))

colnames(scores) = c(names, "Cluster")

######################### Plot the scores: PC1 vs. PC2 #########################
plot_data_all = as.data.frame(scores) %>% mutate(Cluster = as.factor(Cluster))

g1 = ggplot(plot_data_all) + 
  geom_point(data = plot_data_all, 
             aes(x = PC1, y = PC2, col = Cluster)) +
  scale_x_continuous(name = "PC 1") + 
  scale_y_continuous(name = "PC 2") +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle(paste("Principal Component Scores for one Patch - MFPCA")) +
  theme_bw() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))

ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca_all/pdf/PC1vsPC2/PC1_vs_PC2_",pid,"_clustered_all.pdf"), plot = g1, width = 7, height = 4.5)
ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca_all/png/PC1vsPC2/PC1_vs_PC2_",pid,"_clustered_all.png"), plot = g1, width = 7, height = 4.5)

duplicated_grid_cells = all_grid_cells_all[duplicated(all_grid_cells_all)]

plot_data = plot_data_all %>%
  filter(!(rownames(plot_data_all) %in% duplicated_grid_cells)) 
plot_data = plot_data %>%
  mutate(scenario = ifelse(rownames(plot_data) %in% dimnames(picontrol_funData@X[!apply(picontrol_funData@X,1,function(x)all(x==0)),,])[[1]], "Control",
                           ifelse(rownames(plot_data) %in% dimnames(ssp126_funData@X[!apply(ssp126_funData@X,1,function(x)all(x==0)),,])[[1]], "SSP1-RCP2.6",
                                  ifelse(rownames(plot_data) %in% dimnames(ssp370_funData@X[!apply(ssp370_funData@X,1,function(x)all(x==0)),,])[[1]], "SSP3-RCP7.0",
                                         ifelse(rownames(plot_data) %in% dimnames(ssp585_funData@X[!apply(ssp585_funData@X,1,function(x)all(x==0)),,])[[1]], "SSP5-RCP8.5",NA)))))

g2 = ggplot(plot_data) + 
  geom_point(data = plot_data, 
             aes(x = PC1, y = PC2, col = Cluster)) +
  scale_x_continuous(name = "PC 1") + 
  scale_y_continuous(name = "PC 2") +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle(paste("Principal Component Scores for one Patch - MFPCA")) +
  theme_bw() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))

ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca_all/pdf/PC1vsPC2/PC1_vs_PC2_",pid,"_clustered_unique.pdf"), plot = g2, width = 7, height = 4.5)
ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca_all/png/PC1vsPC2/PC1_vs_PC2_",pid,"_clustered_unique.png"), plot = g2, width = 7, height = 4.5)

g3 = ggplot(plot_data) + 
  geom_point(data = plot_data, 
             aes(x = PC1, y = PC2, col = scenario)) +
  scale_x_continuous(name = "PC 1") + 
  scale_y_continuous(name = "PC 2") +
  scale_color_manual(name = "Scenario", values = c("Control" = "#F8766D", "SSP1-RCP2.6" = "#7CAE00" , "SSP3-RCP7.0" = "#00BFC4", "SSP5-RCP8.5" = "#C77CFF", "5" = "darkgrey")) +
  ggtitle(paste("Principal Component Scores for one Patch - MFPCA")) +
  theme_bw() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))

ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca_all/pdf/PC1vsPC2/PC1_vs_PC2_",pid,"_scenarios.pdf"), plot = g3, width = 7, height = 4.5)
ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca_all/png/PC1vsPC2/PC1_vs_PC2_",pid,"_scenarios.png"), plot = g3, width = 7, height = 4.5)


plot_grid(g2,g3)
ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca_all/pdf/PC1vsPC2/PC1_vs_PC2_",pid,"_scenarios_clustered_unique.pdf"), width = 14, height = 7)
ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca_all/png/PC1vsPC2/PC1_vs_PC2_",pid,"_scenarios_clustered_unique.png"), width = 14, height = 7)

plot_grid(g1,g3)
ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca_all/pdf/PC1vsPC2/PC1_vs_PC2_",pid,"_scenarios_clustered_all.pdf"), width = 14, height = 7)
ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca_all/png/PC1vsPC2/PC1_vs_PC2_",pid,"_scenarios_clustered_all.png"), width = 14, height = 7)

######################### Plot principal components ############################

for (iPC in (1:M)){
  pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca_all/pdf/PCs/PCs_PC", iPC, "_commonZlim_", makeZlim, ".pdf"), width = 25, height = 8)
  plot.MFPCAfit2D(MFPCA_all, combined = FALSE, plotPCs = iPC, xlab = paste("Year after Disturbance"), xlim = c(0,100), ylab = "PFT", ylim = c(0.5,5.5), makeZlim = makeZlim, cex = 1.5)
  mtext(paste("PC", iPC, "for all four scenarios"), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
  
  dev.off()
  
  png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca_all/png/PCs/PCs_PC", iPC, "_commonZlim_", makeZlim, ".png"), width = 2500, height = 800)
  plot.MFPCAfit2D(MFPCA_all, combined = FALSE, plotPCs = iPC, xlab = paste("Year after Disturbance"), xlim = c(0,100), ylab = "PFT", ylim = c(0.5,5.5), makeZlim = makeZlim, cex = 1.5)
  mtext(paste("PC", iPC, "for all four scenarios"), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
  
  dev.off()
}


######################## Plot reconstructed curves #############################
# not implemented in the package

######################### Plot Cluster Assignments #############################
## Import fd objects
for (scen in scenarios){
  for (pft in pfts){
    fit.scen_pft = readRDS(paste0("Scripts/FPCA/FdObjects/Wfdobj_", scen, "_", pft, ".rds"))
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
    pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca_all/pdf/clusters/Cluster_", k,"-means_",scen, "_",pft, "_",pid, ".pdf"),width = 10, height = 10)
    par(mfrow = c(2,2))
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 1 -", sum(cluster1), "of", sum(scores[,M+1] == 1), "grid cells"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster1_", scen,"_", pft))], col = "#F8766D", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster1_", scen,"_", pft))]), col = "darkred", lwd = 4, lty = 1)

    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 2 -", sum(cluster2), "of", sum(scores[,M+1] == 2), "grid cells"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster2_", scen,"_", pft))], col = "#7CAE00", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster2_", scen,"_", pft))]), col = "darkgreen", lwd = 4, lty = 1)

    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 3 -", sum(cluster3), "of", sum(scores[,M+1] == 3), "grid cells"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster3_", scen,"_", pft))], col = "#00BFC4", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster3_", scen,"_", pft))]), col = "darkblue", lwd = 4, lty = 1)

    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey",  main = paste("Cluster 4 -", sum(cluster4), "of", sum(scores[,M+1] == 4), "grid cells"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster4_", scen,"_", pft))], col = "#C77CFF", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster4_", scen,"_", pft))]), col = "purple4", lwd = 4, lty = 1)

    mtext(paste(long_names_scenarios(scen), "-", long_names_pfts(tolower(pft)), "-", length(get(paste0("cluster1_", scen, "_", pft))), "grid cells in total"), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
    dev.off()

    # Save as png
    png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/MFPCA/mfpca_all/png/clusters/Cluster_", k,"-means_",scen, "_",pft, "_",pid, ".png"),width = 1000, height = 1000)
    par(mfrow = c(2,2))
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 1 -", sum(cluster1), "of", sum(scores[,M+1] == 1), "grid cells"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster1_", scen,"_", pft))], col = "#F8766D", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster1_", scen,"_", pft))]), col = "darkred", lwd = 4, lty = 1)
    
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 2 -", sum(cluster2), "of", sum(scores[,M+1] == 2), "grid cells"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster2_", scen,"_", pft))], col = "#7CAE00", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster2_", scen,"_", pft))]), col = "darkgreen", lwd = 4, lty = 1)
    
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 3 -", sum(cluster3), "of", sum(scores[,M+1] == 3), "grid cells"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster3_", scen,"_", pft))], col = "#00BFC4", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster3_", scen,"_", pft))]), col = "darkblue", lwd = 4, lty = 1)
    
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey",  main = paste("Cluster 4 -", sum(cluster4), "of", sum(scores[,M+1] == 4), "grid cells"), lty = 1)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster4_", scen,"_", pft))], col = "#C77CFF", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[get(paste0("cluster4_", scen,"_", pft))]), col = "purple4", lwd = 4, lty = 1)
    
    mtext(paste(long_names_scenarios(scen), "-", long_names_pfts(tolower(pft)), "-", length(get(paste0("cluster1_", scen, "_", pft))), "grid cells in total"), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
    dev.off()
    
    print(paste("Plot", scen, pft, "done."))
  }
}

