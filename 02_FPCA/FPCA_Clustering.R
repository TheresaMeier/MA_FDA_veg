################################################################################
############################ Master's Thesis ###################################
################################################################################

############# Exploratory Analysis: FPCA - Clustering of Scores ################

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
pid = 1           # Choose pid (int) or 'all'

#######################

## Import fd objects
for (scen in scenarios){
  for (pft in pfts){
    fit.scen_pft = readRDS(paste0("Scripts/MA_FDA_veg/02_FPCA/FdObjects/Wfdobj_", scen, "_", pft, ".rds"))
    fit.scen_pft_exp = fit.scen_pft
    fit.scen_pft_exp$Wfdobj$coefs = exp(fit.scen_pft_exp$Wfdobj$coefs)
    
    assign(paste0("fit.", scen, "_", pft, "_exp"), fit.scen_pft_exp)
    
    scen.pca_pft = pca.fd(fit.scen_pft_exp$Wfdobj,2)
    
    assign(paste0(scen, ".pca_", pft), scen.pca_pft)
  }
  
  print("Fd objects loaded.")
  # Get scores
  scores.pca.scen <- cbind(get(paste0(scen,".pca_Tundra"))$scores,  
                           get(paste0(scen,".pca_BNE"))$scores,
                           get(paste0(scen,".pca_IBS"))$scores,
                           get(paste0(scen,".pca_otherC"))$scores,
                           get(paste0(scen,".pca_TeBS"))$scores)
  
  colnames(scores.pca.scen) <- c("Tundra_PC1", "Tundra_PC2", "BNE_PC1", "BNE_PC2",
                                 "IBS_PC1", "IBS_PC2", "otherC_PC1", "otherC_PC2",
                                 "TeBS_PC1", "TeBS_PC2")
  
  assign(paste0("scores.pca.", scen), scores.pca.scen)
  
  print( "Scores derived.")
}

set.seed(2)
pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/Clustering/Elbow_",pid,".pdf"), width = 7, height = 7)
par(mfrow = c(2,2))
for (scen in scenarios){
  ## Clustering of the coefficients for different k
  
  scores.pca.scen = get(paste0("scores.pca.", scen))
  wcss <- sapply(1:10, function(k) {
    kmeans(scores.pca.scen, centers = k)$tot.withinss
  })
  
  # Plot the elbow curve
  plot(1:10, wcss, type = "b", pch = 19, frame = TRUE, 
       xlab = "Number of Clusters", ylab = "Within-cluster sum of squares",
       main = long_names_scenarios(scen), cex.main = 1.6, cex.lab = 1.3)
  
  # Add lines for visual aid
  abline(v = 1:10, lty = 3, col = "gray")
}
dev.off()

# Save as png as well
set.seed(2)
png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/Clustering/Elbow_",pid,".png"), width = 1000, height = 1000)
par(mfrow = c(2,2))
for (scen in scenarios){
  ## Clustering of the coefficients for different k
  
  scores.pca.scen = get(paste0("scores.pca.", scen))
  wcss <- sapply(1:10, function(k) {
    kmeans(scores.pca.scen, centers = k)$tot.withinss
  })
  
  # Plot the elbow curve
  plot(1:10, wcss, type = "b", pch = 19, frame = TRUE, 
       xlab = "Number of Clusters", ylab = "Within-cluster sum of squares",
       main = long_names_scenarios(scen))
  
  # Add lines for visual aid
  abline(v = 1:10, lty = 3, col = "gray")
}
dev.off()

## Look at correlations between scores 
library(ggcorrplot)

for (scen in scenarios){
  scores.pca.scen = get(paste0("scores.pca.", scen))
  colnames(scores.pca.scen) = gsub("_", " ", colnames(scores.pca.scen))
  g = ggcorrplot(cor(scores.pca.scen), method = "square", type = "lower", lab = TRUE, 
             title = long_names_scenarios(scen), ggtheme = ggplot2::theme_bw() + ggplot2::theme(text = element_text(size = 15),
                                                                                                plot.title = element_text(size = 20, face = "bold",hjust = 0.5)),
             legend.title = "Correlation", lab_size = 2.7)
  assign(paste0("g_", scen), g)
}

plot_grid(g_picontrol, g_ssp126, g_ssp370, g_ssp585)
ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/FPCA/unrotated/correlations_",pid,".pdf"), width = 10, height = 10)
ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/FPCA/unrotated/correlations_",pid,".png"), width = 10, height = 10)


## Run k-means clustering with 4 clusters
set.seed(1)
par(mfrow = c(1,1))

for (scen in scenarios){
  k <- 4
  scores.pca.scen = get(paste0("scores.pca.", scen))
  
  kmeans_result <- kmeans(scores.pca.scen, centers = k)
  assign(paste0("kmeans_result_", scen), kmeans_result)

  cluster1 = kmeans_result$cluster == 1
  cluster2 = kmeans_result$cluster == 2
  cluster3 = kmeans_result$cluster == 3
  cluster4 = kmeans_result$cluster == 4
  
  #table(kmeans_result$cluster)
  print("Start plotting...")
  for (pft in pfts){

    # Save as pdf
    pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/Clustering/", pft, "/pdf/Cluster_", k,"-means_",scen, "_",pft, "_",pid,".pdf"), width = 8, height = 8)
    par(mfrow = c(2,2))
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 1 -", sum(cluster1), "curves"), lty = 1, cex.main = 1.7, cex.lab = 1.5)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster1], col = "#F8766D", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster1]), col = "darkred", lwd = 4, lty = 1)

    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 2 -", sum(cluster2), "curves"), lty = 1, cex.main = 1.7, cex.lab = 1.5)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster2], col = "#7CAE00", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster2]), col = "darkgreen", lwd = 4, lty = 1)

    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 3 -", sum(cluster3), "curves"), lty = 1, cex.main = 1.7, cex.lab = 1.5)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster3], col = "#00BFC4", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster3]), col = "darkblue", lwd = 4, lty = 1)

    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey",  main = paste("Cluster 4 -", sum(cluster4), "curves"), lty = 1, cex.main = 1.7, cex.lab = 1.5)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster4], col = "#C77CFF", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster4]), col = "purple4", lwd = 4, lty = 1)

    mtext(paste(long_names_scenarios(scen), "-", long_names_pfts(tolower(pft))), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.5)
    dev.off()

    # Save as png
    png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/Clustering/", pft, "/png/Cluster_", k,"-means_",scen, "_",pft, "_",pid,".png"), width = 700, height = 700)
    par(mfrow = c(2,2))
    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 1 -", sum(cluster1), "curves"), lty = 1, cex.main = 1.7, cex.lab = 1.5)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster1], col = "#F8766D", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster1]), col = "darkred", lwd = 4, lty = 1)

    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 2 -", sum(cluster2), "curves"), lty = 1, cex.main = 1.7, cex.lab = 1.5)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster2], col = "#7CAE00", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster2]), col = "darkgreen", lwd = 4, lty = 1)

    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey", main = paste("Cluster 3 -", sum(cluster3), "curves"), lty = 1, cex.main = 1.7, cex.lab = 1.5)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster3], col = "#00BFC4", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster3]), col = "darkblue", lwd = 4, lty = 1)

    plot(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj, col = "grey",  main = paste("Cluster 4 -", sum(cluster4), "curves"), lty = 1, cex.main = 1.7, cex.lab = 1.5)
    lines(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster4], col = "#C77CFF", lwd = 3, lty = 1)
    lines(mean.fd(get(paste0("fit.", scen, "_", pft, "_exp"))$Wfdobj[cluster4]), col = "purple4", lwd = 4, lty = 1)

    mtext(paste(long_names_scenarios(scen), "-", long_names_pfts(tolower(pft))), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.5)
    dev.off()

    print(paste("Plot", scen, pft, "done."))
  }
}

print("All done.")

## Display PC1 vs. PC2 for each Cluster and PFT
for (pft in pfts){
  for (scen in scenarios){
    for (cl in (1:4)){
      data_pca = get(paste0(scen,".pca_", pft))
      cluster = get(paste0("kmeans_result_", scen))$cluster == cl
      
      plot_data = as.data.frame(data_pca$scores[cluster,]) %>%
        rename(PC1 = V1,
               PC2 = V2) %>%
        mutate(Cluster = cl)
      
      assign(paste0("plot_data_cluster", cl), plot_data)
    }
    plot_data = purrr::reduce(list(plot_data_cluster1, plot_data_cluster2, plot_data_cluster3, plot_data_cluster4), bind_rows)
    
    plot_data = plot_data %>% mutate(name = long_names_scenarios(scen))
    assign(paste0("plot_data_", scen), plot_data)
  }
  
  plot_data = purrr::reduce(list(plot_data_picontrol, plot_data_ssp126, plot_data_ssp370, plot_data_ssp585), bind_rows) %>%
    mutate(Cluster = as.factor(Cluster))

  ggplot(plot_data) + 
    geom_point(data = plot_data, 
               aes(x = PC1, y = PC2, col = Cluster), size = 2.5) +
    facet_grid(rows = vars(name)) +
    scale_x_continuous(name = "PC 1") +
    scale_y_continuous(name = "PC 2") +
    scale_color_manual(name = "Cluster", drop = TRUE, values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")) +
    ggtitle(paste("First and second PC scores -", long_names_pfts(tolower(pft)))) +
    theme_bw() +
    theme(text = element_text(size = 11), plot.title = element_text(size = 16, face = "bold",hjust = 0.5))
  
  ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/Clustering/",pft,"/pdf/PC1_vs_PC2_",pft,"_",pid,"_clustered.pdf"), width = 7, height = 4.5)
  ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/Clustering/",pft,"/png/PC1_vs_PC2_",pft,"_",pid,"_clustered.png"), width = 7, height = 4.5)
}

