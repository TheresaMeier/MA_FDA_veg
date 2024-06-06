################################################################################
############################ Master's Thesis ###################################
################################################################################

############## Exploratory Analysis: MFPCA - for all scenarios #################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/02_FPCA/functions.R")
source("Scripts/MA_FDA_veg/03_MFPCA/MFPCA_calculation.R")
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

# Extract features for readibility
mean_function <- MFPCA_all$meanFunction@.Data[[1]]@X
PCs <- MFPCA_all$functions@.Data[[1]]@X
X_orig = funData_all@.Data[[1]]@X

# Center the data by subtracting the mean function
X_orig_centered <- array(NA, dim = c(1803, 100, 5))
for (i in 1:1803) {
  for (j in 1:5){
    X_orig_centered[i,,j] <- X_orig[i,,j] - mean_function[,,j]
  }
}

################# Calculate and cluster temporal scores ########################

for (iYear in seq(10,100,by=10)){
  # Calculate scores manually
  scores_manual <- matrix(0, nrow = 1803, ncol = 10)
  for (i in 1:1803) {
    for (k in 1:10) {
      scores_manual[i, k] <- sum(X_orig_centered[i,1:iYear,] * PCs[k,1:iYear,])
    }
  }
  
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
           # scenario = c(rep("Control", dim(d_picontrol)[1]),rep("SSP1-RCP2.6", dim(d_ssp126)[1]),rep("SSP3-RCP7.0", dim(d_ssp370)[1]),rep("SSP5-RCP8.5", dim(d_ssp585)[1])))
           scenario = c(rep("Control", 434),rep("SSP1-RCP2.6", 442),rep("SSP3-RCP7.0", 462),rep("SSP5-RCP8.5", 465)))
  
  assign(paste0("plot_data_", iYear), plot_data_iYear)  
}

# # Resort clusters if needed
# right_order_126 = c(4,2,1,3)
# plot_data_126$Cluster = right_order_126[plot_data_126$Cluster]

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
for (iYear in seq(10,100,by=10)){
  write.table(plot_data_iYear, paste0("Scripts/Plots/MFPCA/plot_data/plot_data_", iYear, ".txt"), row.names = TRUE, col.names = TRUE)
}

########################## Plot cluster development ############################

for (scen in scenarios){
  
  cluster_list = list()
  
  for (iYear in seq(10,100,by=10)){
    cluster_list[[iYear]] <- get(paste0("plot_data_", iYear))[
      get(paste0("plot_data_", iYear))$scenario == long_names_scenarios(scen), "Cluster"]
  }
  
  sankey_data_scen = as.data.frame(do.call(cbind,cluster_list))
  #rownames(sankey_data_scen) = rownames(d_picontrol)
  colnames(sankey_data_scen) = paste0("y", seq(10,100,by=10)) 
  assign(paste0("sankey_data_", scen), sankey_data_scen)
  
  sankey_data_scen_tmp = sankey_data_scen %>% make_long(y10,y20,y30,y40,y50,y60,y70,y80,y90,y100) 
  
  # Create Sankey plot  
  g = ggplot(sankey_data_scen_tmp, aes(x = x, 
                                       next_x = next_x, 
                                       node = node, 
                                       next_node = next_node,
                                       fill = factor(node)))  +
    labs(x = "Years after disturbance",
         y = NULL,
         fill = "Cluster",
         color = NULL) +
    geom_sankey(flow.alpha = .6,
                node.color = "gray30") +
    #geom_sankey_label(size = 3, color = "white", fill = "gray40") +
    scale_fill_discrete(drop=FALSE) + ggtitle(long_names_scenarios(scen)) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          text = element_text(size = 15),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_x_discrete(labels = seq(10,100,by=10))
  
  assign(paste0("sankey_plot_", scen, "_sorted"),g)
}

plot_grid(sankey_plot_picontrol_sorted, sankey_plot_ssp126_sorted, sankey_plot_ssp370_sorted, sankey_plot_ssp585_sorted, nrow = 2)

ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/clusters_with_years_",pid,".pdf"), width = 15, height = 10)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/clusters_with_years_",pid,".png"), width = 15, height = 10)

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

ari_values <- data.frame(
  Year = rep(seq(10,100,by=10)[-1], 4),  # Exclude the first year as ARI is calculated between consecutive years
  ARI = c(
    round(adjusted_rand_for_consecutive_columns(sankey_data_picontrol), 2),
    round(adjusted_rand_for_consecutive_columns(sankey_data_ssp126), 2),
    round(adjusted_rand_for_consecutive_columns(sankey_data_ssp370), 2),
    round(adjusted_rand_for_consecutive_columns(sankey_data_ssp585), 2)
  ),
  Scenario = rep(c("picontrol", "ssp126", "ssp370", "ssp585"), each = 9)
)

ari_values$Scenario = long_names_scenarios(ari_values$Scenario)

ggplot(ari_values, aes(x = Year, y = ARI, color = Scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, color = "darkgrey", size = 0.5) +
  labs(title = "Adjusted Rand Index (ARI) Over Time for Different Scenarios",
       x = "Years after disturbance",
       y = "Adjusted Rand Index (ARI)",
       color = "Scenario") +
  theme_bw() + theme(text = element_text(size = 15),plot.title = element_text(size = 20, face = "bold",hjust = 0.5)) +
  scale_x_continuous(breaks = seq(10,100,by=10)[-1])

ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/ARI_clusters_with_years_",pid,".pdf"), width = 15, height = 10)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/ARI_clusters_with_years_",pid,".png"), width = 15, height = 10)


####################### Get temporal cluster description #######################

for (pft in pfts) 
  for (cl in c(1:4)) assign(paste0(pft,"_cl", cl), rep(0,10))

for (iPFT in c(1:5)){
  d_pft = funData_all@.Data[[1]]@X[,,iPFT]
  i=1
  for (iYear in seq(10,100,by=10)){
    plot_data_iYear = get(paste0("plot_data_", iYear))
    
    for (iCl in c(1:4)) {
      
      pft_cl = get(paste0(pfts[iPFT], "_cl", iCl))
      
      pft_cl[i] = mean(d_pft[plot_data_iYear$Cluster == iCl,c(1:iYear)])
      assign(paste0(pfts[iPFT],"_cl", iCl), pft_cl)
    }
    i=i+1
  }
}

# Create full data set for plotting
for (pft in pfts){
  d_pft = as.data.frame(cbind(c(get(paste0(pft, "_cl1")), get(paste0(pft, "_cl2")), get(paste0(pft, "_cl3")), get(paste0(pft, "_cl4"))), rep(c(5,seq(10,120,by=10),126),4), c(rep(1,14),rep(2,14),rep(3,14), rep(4,14)), rep(pft,14*4)))
  colnames(d_pft) = c("value", "year", "Cluster", "PFT")
  assign(paste0("d_", pft), d_pft)
}

d_pft = rbind(d_Tundra, d_BNE, d_IBS, d_otherC, d_TeBS) %>%
  mutate(value = as.numeric(value),
         year = as.numeric(year),
         Cluster = as.factor(Cluster),
         PFT = long_names_pfts(tolower(PFT)))


# Plot the results


ggplot(d_pft) + 
  geom_line(data = d_pft,
            aes(x = year, y = value, color = Cluster), lwd = 2) +
  facet_grid(rows = vars(PFT)) + 
  scale_y_continuous(name = "Share of above ground carbon") + 
  scale_x_continuous(name = "Year after Disturbance", breaks = seq(10,100,by=10)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  geom_vline(xintercept = seq(10,100,by=10), linetype = "dashed", color = "black", lwd = 0.5) +
  theme_bw() + theme(plot.title = element_text(hjust = .5)) +
  ggtitle("Cluster-wise mean share of above ground carbon over time")

ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/cluster_means_with_years_",pid,".pdf"), width = 15, height = 10)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/cluster_means_with_years_",pid,".png"), width = 15, height = 10)

