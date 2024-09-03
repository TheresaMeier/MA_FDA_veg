################################################################################
############################ Master's Thesis ###################################
################################################################################

############## Exploratory Analysis: MFPCA - for all scenarios #################

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
library(foreach)
library(funData)
library(ggsankey)

# Get MFPCA and multiFunData object
MFPCA_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_all_1803_100y.rds")
funData_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_1803.rds")
M = 10      # 10 PCs
scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")
pfts = c("Tundra", "BNE", "IBS", "otherC", "TeBS")
pid = 1

# Center the data by subtracting the mean function
funData_centered = funData_all
for (iPFT in 1:5){
  for (iCurve in 1:1803) {
    funData_centered@.Data[[iPFT]]@X[iCurve,] = funData_centered@.Data[[iPFT]]@X[iCurve,] - MFPCA_all$meanFunction@.Data[[iPFT]]@X
  }
}

################# Calculate and cluster temporal scores ########################
# Cluster assignment of previous clustering for consistency
right_order_all = c(4,1,3,2)  # for control only: c(1,4,2,3)
for (iYear in c(seq(10,100,by=10))){
  # Calculate scores manually
  scores_manual = funData_centered@.Data[[1]]@X[,iYear] %*% t(MFPCA_all$functions[[1]]@X[,iYear]) +
    funData_centered@.Data[[2]]@X[,iYear] %*% t(MFPCA_all$functions[[2]]@X[,iYear])+
    funData_centered@.Data[[3]]@X[,iYear] %*% t(MFPCA_all$functions[[3]]@X[,iYear]) +
    funData_centered@.Data[[4]]@X[,iYear] %*% t(MFPCA_all$functions[[4]]@X[,iYear])+
    funData_centered@.Data[[5]]@X[,iYear] %*% t(MFPCA_all$functions[[5]]@X[,iYear])
  
  # Cluster temporal scores
  set.seed(1)
  kmeans_iYear = kmeans(scores_manual, 4)
  scores = cbind(scores_manual, kmeans_iYear$cluster)
  
  names = c("PC1", "PC2")
  for (i in (3:M)) names = c(names, paste0("PC", i))
  
  colnames(scores) = c(names, "Cluster")
  
  assign(paste0("scores_", iYear), scores)
  
  plot_data_iYear = as.data.frame(scores) %>% 
    mutate(Cluster = as.factor(Cluster),
           scenario = c(rep("Control", 434),rep("SSP1-RCP2.6", 442),rep("SSP3-RCP7.0", 462),rep("SSP5-RCP8.5", 465)))
  
  plot_data_iYear$Cluster = right_order_all[plot_data_iYear$Cluster]
  assign(paste0("plot_data_", iYear), plot_data_iYear)  
}

# Resort clusters if needed
combis = cbind(c(rev(seq(20,100,by=10))), c(rev(seq(10,90,by=10))))

for (cb in 1:nrow(combis)){
  y1_data = get(paste0("plot_data_", combis[cb,1]))
  y2_data = get(paste0("plot_data_", combis[cb,2]))
  
  lens = matrix(0, nrow = 4, ncol = 4)
  for (iCl1 in c(1:4)){
    for (iCl2 in c(1:4)){
      lens[iCl1, iCl2] = length(intersect(rownames(y1_data[y1_data$Cluster==iCl1,]), rownames(y2_data[y2_data$Cluster==iCl2,])))
    }
  }
  
  right_order = max.col(t(lens), "first")
  
  if (!any(duplicated(right_order))) y2_data$Cluster = right_order[y2_data$Cluster]
  assign(paste0("plot_data_", combis[cb,2]), y2_data)
  print(cb)
}

# Save resulting plot data
for (iYear in c(seq(10,100,by=10))){
  write.table(plot_data_iYear, paste0("Scripts/Plots/MFPCA/plot_data/plot_data_", iYear, "_single.txt"), row.names = TRUE, col.names = TRUE)
}

########################## Plot cluster development ############################

for (scen in scenarios){
  
  cluster_list = list()
  
  for (iYear in c(seq(10,100,by=10))){
    cluster_list[[iYear]] <- get(paste0("plot_data_", iYear))[
      get(paste0("plot_data_", iYear))$scenario == long_names_scenarios(scen), "Cluster"]
  }
  
  sankey_data_scen = as.data.frame(do.call(cbind,cluster_list))
  colnames(sankey_data_scen) = paste0("y", c(seq(10,100,by=10))) 
  assign(paste0("sankey_data_", scen), sankey_data_scen)
  
  sankey_data_scen_tmp = sankey_data_scen %>% make_long(y10,y20,y30,y40,y50,y60,y70,y80,y90,y100) 
  
  # Create Sankey plot  
  g = ggplot(sankey_data_scen_tmp, aes(x = x, 
                                       next_x = next_x, 
                                       node = node, 
                                       next_node = next_node,
                                       fill = factor(node)))  +
    labs(x = "Year after disturbance",
         y = NULL,
         fill = "Cluster",
         color = NULL) +
    geom_sankey(flow.alpha = .6,
                node.color = "gray30") +
    #geom_sankey_label(size = 3, color = "white", fill = "gray40") +
    scale_fill_discrete(drop=FALSE) + ggtitle(paste(long_names_scenarios(scen))) +
    theme_bw() +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          text = element_text(size = 14),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_x_discrete(labels = c(seq(10,100,by=10)))
  
  assign(paste0("sankey_plot_", scen, "_sorted"),g)
  ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/clusters_with_years_",pid,"_single_", scen, ".pdf"), g, width = 5, height = 3)
  ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/clusters_with_years_",pid,"_single_", scen, ".png"), g, width = 5, height = 3)
  
}

# For all scenarios together
cluster_list = list()

for (iYear in c(seq(10,100,by=10))){
  cluster_list[[iYear]] <- get(paste0("plot_data_", iYear))[, "Cluster"]
}

sankey_data = as.data.frame(do.call(cbind,cluster_list))
colnames(sankey_data) = paste0("y", c(seq(10,100,by=10))) 

sankey_data_tmp = sankey_data %>% make_long(y10,y20,y30,y40,y50,y60,y70,y80,y90,y100) 

# Create Sankey plot  
ggplot(sankey_data_tmp, aes(x = x, 
                                     next_x = next_x, 
                                     node = node, 
                                     next_node = next_node,
                                     fill = factor(node)))  +
  labs(x = "Year after disturbance",
       y = NULL,
       fill = "Cluster",
       color = NULL) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  #geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_discrete(drop=FALSE) + ggtitle("All scenarios combined") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        text = element_text(size = 14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_discrete(labels = c(seq(10,100,by=10)))

ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/clusters_with_years_",pid,"_single_combined.pdf"), width = 5, height = 3)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/clusters_with_years_",pid,"_single_combined.png"), width = 5, height = 3)


# Check Adjusted Rand Index (ARI) for Cluster consistency:
library(mclust)
adjusted_rand_for_consecutive_columns <- function(df) {
  n <- ncol(df)
  results <- numeric(n - 1)  # to store results
  for (i in 1:(n - 1)) {
    results[i] <- adjustedRandIndex(df[[i]], df[[i + 1]])
  }
  return(results)
}

pastel_colors <- c("Control" = "skyblue2",  # Light blue
                   "SSP1-RCP2.6" = "#77DD77",    # Light green
                   "SSP3-RCP7.0" = "#FFB347",    # Light orange
                   "SSP5-RCP8.5" = "plum2",    # Light purple
                   "Combined" = "grey30")  # Light red

ari_values <- data.frame(
  Year = rep(c(seq(10,100,by=10))[-1], 5),  # Exclude the first year as ARI is calculated between consecutive years
  ARI = c(
    round(adjusted_rand_for_consecutive_columns(sankey_data_picontrol), 2),
    round(adjusted_rand_for_consecutive_columns(sankey_data_ssp126), 2),
    round(adjusted_rand_for_consecutive_columns(sankey_data_ssp370), 2),
    round(adjusted_rand_for_consecutive_columns(sankey_data_ssp585), 2),
    round(adjusted_rand_for_consecutive_columns(sankey_data), 2)
  ),
  Scenario = rep(c("picontrol", "ssp126", "ssp370", "ssp585", "Combined"), each = 9)
)

ari_values$Scenario = long_names_scenarios(ari_values$Scenario)
ari_values$Scenario <- factor(ari_values$Scenario, levels = c("Control", "SSP1-RCP2.6", "SSP3-RCP7.0", "SSP5-RCP8.5", "Combined"))


ggplot(ari_values, aes(x = Year, y = ARI, color = Scenario)) +
  geom_line(aes(linewidth = ifelse(Scenario == "Combined", 2, 1))) +  # Adjust linewidth based on Scenario
  geom_point(size = 2) +
  geom_hline(yintercept = 1, color = "darkgrey", size = 0.5) +
  labs(title = "ARI over time for different scenarios",
       x = "Year after disturbance",
       y = "ARI",
       color = "Scenario") +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  scale_x_continuous(breaks = c(seq(10, 100, by = 10))[-1]) +
  scale_linewidth_identity() +
  scale_color_manual(values = pastel_colors)

ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/ARI_clusters_with_years_",pid,"_single.pdf"), width = 10, height = 6)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/ARI_clusters_with_years_",pid,"_single.png"), width = 10, height = 6)
