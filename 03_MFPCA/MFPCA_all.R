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
createFunData = FALSE
runMFPCA = FALSE
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
# Create funData object 
if(createFunData){
  for (scen in scenarios){
    print(paste("Start with scenario", long_names_scenarios(scen)))
    
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
  
  d_all = abind(d_picontrol, d_ssp126, d_ssp370, d_ssp585, along = 1)
  funData_tmp = funData(argvals = list(1:126, 1:5), X = d_all)
  
  funData_all = multiFunData(funData_tmp)
  saveRDS(funData_all, "Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_1803.rds")
  saveRDS(funData_tmp, "Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_tmp_1803.rds")
  
  
  print("... all done.")
}

# Run MFPCA
if (runMFPCA){
  funData_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_1803.rds")
  funData_tmp = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_tmp_1803.rds")
  
  for (iYear in c(5,seq(10,120,by=10),126)){
    funData_iYear = funData_all
    funData_iYear@.Data[[1]]@argvals[[1]] = c(1:iYear)
    funData_iYear@.Data[[1]]@X = funData_iYear@.Data[[1]]@X[,c(1:iYear),]
    
    funData_tmp_iYear = funData_tmp
    funData_tmp_iYear@argvals[[1]] = c(1:iYear)
    funData_tmp_iYear@X = funData_tmp_iYear@X[,c(1:iYear),]
    
    MFPCA_iYear = MFPCA2(funData_iYear, M=M, uniExpansions = list(list(type = "given", funData_tmp_iYear, ortho = TRUE)), fit = TRUE, approx.eigen = FALSE)
    saveRDS(MFPCA_iYear, paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_all_1803_noApprox_allEVs_", iYear,".rds"))
    assign(paste0("MFPCA_", iYear), MFPCA_iYear)
    print(paste("... year", iYear, "done."))
  }
  
}

# Get multiFunData objects
funData_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/funData_all_1803.rds")

# Get MFPCA results
# MFPCA_all = readRDS("Scripts/FPCA/FdObjects/MFPCA_all_1803.rds")
MFPCA_all = readRDS("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_all_1803_noApprox_allEVs_126.rds")

for (iYear in c(5,seq(10,120,by=10),126)){
  assign(paste0("MFPCA_", iYear), readRDS(paste0("Scripts/MA_FDA_veg/03_MFPCA/FdObjects/MFPCA_all_1803_noApprox_allEVs_", iYear, ".rds")))
}

# Compute WCSS for different values of k
wcss <- sapply(1:10, function(k) {
  kmeans(MFPCA_all$scores, centers = k)$tot.withinss
})

# Save elbow plot as pdf
pdf(paste0("Scripts/Plots/MFPCA/Clustrers/pdf/elbow_plot_", pid, ".pdf"), width = 12, height = 8)
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
png(paste0("Scripts/Plots/MFPCA/Clustrers/pdf/elbow_plot_", pid, ".png"), width = 1200, height = 800)
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

# Clustering of the coefficients
kmeans_result = kmeans(MFPCA_all$scores,4)
scores = cbind(MFPCA_all$scores, kmeans_result$cluster)

names = c("PC1", "PC2")
for (i in (3:M)) names = c(names, paste0("PC", i))

colnames(scores) = c(names, "Cluster")
plot_data_all = as.data.frame(scores) %>% 
  mutate(Cluster = as.factor(Cluster),
         # scenario = c(rep("Control", dim(d_picontrol)[1]),rep("SSP1-RCP2.6", dim(d_ssp126)[1]),rep("SSP3-RCP7.0", dim(d_ssp370)[1]),rep("SSP5-RCP8.5", dim(d_ssp585)[1])))
         scenario = c(rep("Control", 434),rep("SSP1-RCP2.6", 442),rep("SSP3-RCP7.0", 462),rep("SSP5-RCP8.5", 465)))

right_order_all = c(4,2,1,3)
plot_data_all$Cluster = right_order_all[plot_data_all$Cluster]

write.table(plot_data_all, "Scripts/Plots/MFPCA/plot_data/plot_data_all.txt", row.names = TRUE, col.names = TRUE)

############################## Temporal Clustering #############################

for (iYear in c(5,seq(10,120,by=10),126)){
  set.seed(1)
  MFPCA_tmp = get(paste0("MFPCA_", iYear))
  kmeans_iYear = kmeans(MFPCA_tmp$scores, k)
  scores = cbind(MFPCA_tmp$scores,kmeans_iYear$cluster)
  
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

# Resort clusters if needed
right_order_126 = c(4,2,1,3)
plot_data_126$Cluster = right_order_126[plot_data_126$Cluster]

combis = cbind(c(rev(c(seq(10,120,by=10),126))), rev(c(5,seq(10,120,by=10))))

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
  
  
########################## Plot cluster development ############################
library(ggsankey)

for (scen in scenarios){
  
  cluster_list = list()
  
  for (iYear in c(5,seq(10,120,by=10),126)){
    cluster_list[[iYear]] <- get(paste0("plot_data_", iYear))[
      get(paste0("plot_data_", iYear))$scenario == long_names_scenarios(scen), "Cluster"]
  }
  
  sankey_data_scen = as.data.frame(do.call(cbind,cluster_list))
  #rownames(sankey_data_scen) = rownames(d_picontrol)
  colnames(sankey_data_scen) = paste0("y", c(5,seq(10,120,by=10),126)) 
  assign(paste0("sankey_data_", scen), sankey_data_scen)
  
  sankey_data_scen_tmp = sankey_data_scen %>% make_long(y5,y10,y20,y30,y40,y50,y60,y70,y80,y90,y100,y110,y120,y126) 

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
    theme(plot.title = element_text(hjust = .5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_x_discrete(labels = c(5,seq(10,120,by=10),126))
  
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
  Year = rep(c(5, seq(10, 120, by=10), 126)[-1], 4),  # Exclude the first year as ARI is calculated between consecutive years
  ARI = c(
    round(adjusted_rand_for_consecutive_columns(sankey_data_picontrol), 2),
    round(adjusted_rand_for_consecutive_columns(sankey_data_ssp126), 2),
    round(adjusted_rand_for_consecutive_columns(sankey_data_ssp370), 2),
    round(adjusted_rand_for_consecutive_columns(sankey_data_ssp585), 2)
  ),
  Scenario = rep(c("picontrol", "ssp126", "ssp370", "ssp585"), each = 13)
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
  theme_bw() + theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5)) +
  scale_x_continuous(breaks = c(5, seq(10, 120, by=10), 126)[-1])

ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/ARI_clusters_with_years_",pid,".pdf"), width = 15, height = 10)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/ARI_clusters_with_years_",pid,".png"), width = 15, height = 10)


####################### Get temporal cluster description #######################

for (pft in pfts) 
  for (cl in c(1:4)) assign(paste0(pft,"_cl", cl), rep(0,14))

for (iPFT in c(1:5)){
  d_pft = funData_all@.Data[[1]]@X[,,iPFT]
  i=1
  for (iYear in c(5,seq(10,120,by=10),126)){
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
  scale_x_continuous(name = "Year after Disturbance", breaks = c(5,seq(10,120,by=10),126)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  geom_vline(xintercept = c(5, seq(10, 120, by = 10), 126), linetype = "dashed", color = "black", lwd = 0.5) +
  theme_bw() + theme(plot.title = element_text(hjust = .5)) +
  ggtitle("Cluster-wise mean share of above ground carbon over time")

ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/cluster_means_with_years_",pid,".pdf"), width = 15, height = 10)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/cluster_means_with_years_",pid,".png"), width = 15, height = 10)


################## Cluster description for the whole time span ##################

list_var = c()
list_pft = c()
list_cl = c()
for (iPFT in c(1:5)){
  d_pft = funData_all@.Data[[1]]@X[,,iPFT]
    
  for (iCl in c(1:4)) {
    list_var = c(list_var, colMeans(d_pft[plot_data_126$Cluster == iCl,]))
    list_cl = c(list_cl, rep(iCl, 126))
  }
  
  list_pft = c(list_pft, rep(pfts[iPFT],126*4))
}

d_plot = as.data.frame(cbind("year" = rep(c(1:126),20), "value" = list_var, "Cluster" = list_cl, "PFT" = long_names_pfts(tolower(list_pft)))) %>%
  mutate(year = as.numeric(year),
         value = as.numeric(value),
         Cluster = as.factor(Cluster),
         PFT = as.factor(PFT))

# Plot the mean curves per cluster and per PFT

ggplot(d_plot) + 
  geom_line(data = d_plot, aes(x = year, y = value, col = Cluster), lwd = 2) +
  facet_grid(rows = vars(PFT)) + 
  scale_y_continuous(name = "Share of above ground carbon") + 
  scale_x_continuous(name = "Year after Disturbance", breaks = seq(0,100,by=10), limits = c(0,100)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  theme_bw() + theme(plot.title = element_text(hjust = .5)) +
  ggtitle("Cluster-wise mean share of above ground carbon over time")

ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/cluster_means_per_cluster_",pid,".pdf"), width = 15, height = 10)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/cluster_means_per_cluster_",pid,".png"), width = 15, height = 10)



d_plot = d_plot %>%
  mutate(Cluster = paste("Cluster", Cluster))

ggplot(d_plot) + 
  geom_line(data = d_plot, aes(x = year, y = value, col = PFT), lwd = 2) +
  facet_grid(rows = vars(Cluster)) + 
  scale_y_continuous(name = "Share of above ground carbon") + 
  scale_x_continuous(name = "Year after Disturbance", breaks = seq(0,100,by=10), limits = c(0,100)) +
  scale_color_manual(name = "Dominant vegetation", drop = TRUE,
                     values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                                "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  theme_bw() + theme(plot.title = element_text(hjust = .5)) +
  ggtitle("PFT-wise mean share of above ground carbon over time")


ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/cluster_means_per_PFT_",pid,".pdf"), width = 15, height = 10)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/cluster_means_per_PFT_",pid,".png"), width = 15, height = 10)



######################### Plot the scores: PC1 vs. PC2 #########################
plot_data_all = as.data.frame(scores_126) %>% 
  mutate(Cluster = as.factor(Cluster),
         # scenario = c(rep("Control", dim(d_picontrol)[1]),rep("SSP1-RCP2.6", dim(d_ssp126)[1]),rep("SSP3-RCP7.0", dim(d_ssp370)[1]),rep("SSP5-RCP8.5", dim(d_ssp585)[1])))
         scenario = c(rep("Control", 434),rep("SSP1-RCP2.6", 442),rep("SSP3-RCP7.0", 462),rep("SSP5-RCP8.5", 465)))

g1 = ggplot(plot_data_all) + 
  geom_point(data = plot_data_all, 
             aes(x = PC1, y = PC2, col = Cluster)) +
  scale_x_continuous(name = "PC 1") + 
  scale_y_continuous(name = "PC 2") +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle(paste("Principal Component Scores for one Patch - MFPCA")) +
  theme_bw() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))

ggsave(paste0("Scripts/Plots/MFPCA/PC1vsPC2/pdf/PC1_vs_PC2_",pid,"_clustered_all.pdf"), plot = g1, width = 7, height = 4.5)
ggsave(paste0("Scripts/Plots/MFPCA/PC1vsPC2/png/PC1_vs_PC2_",pid,"_clustered_all.png"), plot = g1, width = 7, height = 4.5)



g2 = ggplot(plot_data_all) + 
  geom_point(data = plot_data_all, 
             aes(x = PC1, y = PC2, col = scenario)) +
  scale_x_continuous(name = "PC 1") + 
  scale_y_continuous(name = "PC 2") +
  scale_color_manual(name = "Scenario", values = c("Control" = "#F8766D", "SSP1-RCP2.6" = "#7CAE00" , "SSP3-RCP7.0" = "#00BFC4", "SSP5-RCP8.5" = "#C77CFF", "5" = "darkgrey")) +
  ggtitle(paste("Principal Component Scores for one Patch - MFPCA")) +
  theme_bw() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))

ggsave(paste0("Scripts/Plots/MFPCA/PC1vsPC2/pdf/PC1_vs_PC2_",pid,"_scenarios.pdf"), plot = g2, width = 7, height = 4.5)
ggsave(paste0("Scripts/Plots/MFPCA/PC1vsPC2/png/PC1_vs_PC2_",pid,"_scenarios.png"), plot = g2, width = 7, height = 4.5)


plot_grid(g1,g2)
ggsave(paste0("Scripts/Plots/MFPCA/PC1vsPC2/pdf/PC1_vs_PC2_",pid,"_scenarios_clustered_all.pdf"), width = 14, height = 7)
ggsave(paste0("Scripts/Plots/MFPCA/PC1vsPC2/png/PC1_vs_PC2_",pid,"_scenarios_clustered_all.png"), width = 14, height = 7)

######################### Plot principal components ############################

# Plot mean function
pdf("Scripts/Plots/MFPCA/PCs/pdf/mean.pdf", width = 15, height = 5)
plot(MFPCA_all$meanFunction, xlab = paste("Year after Disturbance"), xlim = c(1,100), ylab = "PFT", ylim = c(0.5,5.5), cex = 1.5, main = "Mean function - MFPCA")
dev.off()

png("Scripts/Plots/MFPCA/PCs/png/mean.png", width = 1500, height = 500)
plot(MFPCA_all$meanFunction, xlab = paste("Year after Disturbance"), xlim = c(1,100), ylab = "PFT", ylim = c(0.5,5.5), cex = 1.5, main = "Mean function - MFPCA")
dev.off()


# Plot deviations from the mean for each PC
for (iPC in (1:M)){
  pdf(paste0("Scripts/Plots/MFPCA/PCs/pdf/PCs_PC", iPC, ".pdf"), width = 15, height = 10)
  plot.MFPCAfit2D(MFPCA_all, combined = FALSE, plotPCs = iPC, xlab = paste("Year after Disturbance"), xlim = c(1,100), ylab = "PFT", ylim = c(0.5,5.5), cex = 1.5)
  mtext(paste("PC", iPC, "for all four scenarios"), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
  
  dev.off()
  
  png(paste0("Scripts/Plots/MFPCA/PCs/png/PCs_PC", iPC, ".png"), width = 1500, height = 1000)
  plot.MFPCAfit2D(MFPCA_all, combined = FALSE, plotPCs = iPC, xlab = paste("Year after Disturbance"), xlim = c(1,100), ylab = "PFT", ylim = c(0.5,5.5), cex = 1.5)
  mtext(paste("PC", iPC, "for all four scenarios"), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
  
  dev.off()
}

######################## Plot reconstructed curves #############################
# not implemented in the package, but self-implemented:

bounds = t(matrix(c(1,434,434+1,434+1+441,434+1+441+1, 434+1+441+1+461, 434+1+441+1+461+1, 434+1+441+1+461+1+464), ncol = 4))

for (scen in scenarios){
  i = 1
  j = 1
  print(scen)
  for (pft in pfts){
    print(pft)
    # Get original curves
    ## Import fd object
    fit.scen.pft = readRDS(paste0("Scripts/MA_FDA_veg/02_FPCA/FdObjects/Wfdobj_", scen, "_", pft, ".rds"))
    fit.scen.pft_exp = fit.scen.pft
    fit.scen.pft_exp$Wfdobj$coefs = exp(fit.scen.pft_exp$Wfdobj$coefs)
    
    ## Run FPCA
    
    fit.pca = pca.fd(fit.scen.pft_exp$Wfdobj,3)
    
    for (nPC in 2:10){
      print(nPC)
      funs = MFPCA_all$functions[[1]]
      funs@X = funs@X[1:nPC,,]
      univExp <- univExpansion(type = "default", 
                               scores = as.matrix(MFPCA_all$normFactors[1:nPC]/MFPCA_all$scores[,1:nPC]),
                               argvals = MFPCA_all$functions[[1]]@argvals,
                               functions = funs)
      
      multiExp = multiFunData(univExp) + MFPCA_all$meanFunction
      
      # plot(multiExp@.Data[[1]]@X[180,,1], type = "l", ylim = c(0,1), xlim = c(0,100))
      
      ## Plot reconstructed and original fits for using 2,...,10 PCs
  
      # Save as pdf
      pdf(paste0("Scripts/Plots/MFPCA/Reconstruction/pdf/", scen, "/", pft, "/MFPCA_reconstruct_",scen, "_",pft, "_",pid,"_",nPC,"PCs.pdf"), width = 10, height = 4.5)
      par(mfrow = c(1,2))
      plot(x = c(1:100),y = rep(0,100), xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',lty = 2, xlab = "Year after Disturbance", ylab = "Share of above ground carbon", main = paste0("Original fit - ", long_names_scenarios(scen), " - ", long_names_pfts(tolower(pft))))
      for (icurve in 1:nrow(fit.pca$scores)){
        lines(fit.scen.pft_exp$Wfdobj[icurve], col = pal(450)[icurve])
      }
      plot(x = c(1:100),y = rep(0,100), xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',lty = 2, xlab = "Year after Disturbance", ylab = "Share of above ground carbon", main = paste("Reconstructed fit using", nPC, "PCs"))
      
      for (icurve in bounds[j,1]:bounds[j,2]){
        lines(multiExp@.Data[[1]]@X[icurve,,i], col = pal(450)[icurve])
      }
      dev.off()
      
      # Save single plots
      pdf(paste0("Scripts/Plots/MFPCA/Reconstruction/pdf/", scen, "/", pft, "/MFPCA_reconstruct_",scen, "_",pft, "_",pid,"_orig.pdf"), width = 10, height = 8)
      par(mfrow = c(1,1))
      plot(x = c(1:100),y = rep(0,100), xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',lty = 2, xlab = "Year after Disturbance", ylab = "Share of above ground carbon", main = paste0("Original fit - ", long_names_scenarios(scen), " - ", long_names_pfts(tolower(pft))))
      for (icurve in 1:nrow(fit.pca$scores)){
        lines(fit.scen.pft_exp$Wfdobj[icurve], col = pal(450)[icurve])
      }
      dev.off()
      
      pdf(paste0("Scripts/Plots/MFPCA/Reconstruction/pdf/", scen, "/", pft, "/MFPCA_reconstruct_",scen, "_",pft, "_",pid,"_",nPC,"PCs_only.pdf"), width = 10, height = 8)
      plot(x = c(1:100),y = rep(0,100), xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',lty = 2, xlab = "Year after Disturbance", ylab = "Share of above ground carbon", main = paste("Reconstructed fit using", nPC, "PCs"))
      
      for (icurve in bounds[j,1]:bounds[j,2]){
        lines(multiExp@.Data[[1]]@X[icurve,,i], col = pal(450)[icurve])
      }
      dev.off()
      
      # # Save as png
      png(paste0("Scripts/Plots/MFPCA/Reconstruction/png/", scen, "/", pft, "/MFPCA_reconstruct_",scen, "_",pft, "_",pid,"_",nPC,"PCs.png"), width = 1000, height = 450)
      
      par(mfrow = c(1,2))
      plot(x = c(1:100),y = rep(0,100), xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',lty = 2, xlab = "Year after Disturbance", ylab = "Share of above ground carbon", main = paste0("Original fit - ", long_names_scenarios(scen), " - ", long_names_pfts(tolower(pft))))
      for (icurve in 1:nrow(fit.pca$scores)){
        lines(fit.scen.pft_exp$Wfdobj[icurve], col = pal(450)[icurve])
      }
      plot(x = c(1:100),y = rep(0,100), xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',lty = 2, xlab = "Year after Disturbance", ylab = "Share of above ground carbon", main = paste("Reconstructed fit using", nPC, "PCs"))
      
      for (icurve in bounds[j,1]:bounds[j,2]){
        lines(multiExp@.Data[[1]]@X[icurve,,i], col = pal(450)[icurve])
      }
      dev.off()
      
      # Save single plots
      png(paste0("Scripts/Plots/MFPCA/Reconstruction/png/", scen, "/", pft, "/MFPCA_reconstruct_",scen, "_",pft, "_",pid,"_orig.png"), width = 1000, height = 800)
      par(mfrow = c(1,1))
      plot(x = c(1:100),y = rep(0,100), xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',lty = 2, xlab = "Year after Disturbance", ylab = "Share of above ground carbon", main = paste0("Original fit - ", long_names_scenarios(scen), " - ", long_names_pfts(tolower(pft))))
      for (icurve in 1:nrow(fit.pca$scores)){
        lines(fit.scen.pft_exp$Wfdobj[icurve], col = pal(450)[icurve])
      }
      dev.off()
      
      png(paste0("Scripts/Plots/MFPCA/Reconstruction/png/", scen, "/", pft, "/MFPCA_reconstruct_",scen, "_",pft, "_",pid,"_",nPC,"PCs_only.png"), width = 1000, height = 800)
      plot(x = c(1:100),y = rep(0,100), xlim = c(0,100), ylim = c(-0.05,1.2), type = 'l',lty = 2, xlab = "Year after Disturbance", ylab = "Share of above ground carbon", main = paste("Reconstructed fit using", nPC, "PCs"))
      
      for (icurve in bounds[j,1]:bounds[j,2]){
        lines(multiExp@.Data[[1]]@X[icurve,,i], col = pal(450)[icurve])
      }
      dev.off()
      
    }
    
    i = i+1
  }
  j = j+1
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
    png(paste0("Scripts/Plots/MFPCA/Clusters/png/Cluster_", k,"-means_",scen, "_",pft, "_",pid, ".png"),width = 1000, height = 1000)
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
############################ Screeplot of MFPCA ################################

pdf(paste0("Scripts/Plots/MFPCA/PCs/pdf/screeplot_",pid, ".pdf"),width = 10, height = 8)
screeplot(MFPCA_all, main = "Screeplot of MFPCA", cex = 1.3)
dev.off()

png(paste0("Scripts/Plots/MFPCA/PCs/png/screeplot_",pid, ".png"),width = 1000, height = 800)
screeplot(MFPCA_all, main = "Screeplot of MFPCA", cex = 1.3)
dev.off()


