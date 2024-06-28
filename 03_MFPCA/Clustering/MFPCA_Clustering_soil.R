################################################################################
############################ Master's Thesis ###################################
################################################################################

########################### Clustering: Soil data ##############################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/01_Description/utils.R")

## Load libraries
library(dplyr)
library(data.table)
library(stringr)
library(cowplot)
library(ggplot2)
library(ggridges)
library(forcats)

scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")

## Load data
pid = 1

# Get cluster assignments
plot_data_all = read.table("Scripts/Plots/MFPCA/plot_data/plot_data_all.txt")

for (scen in scenarios){
  d_scen_2015_2040 = fread(paste0("Data/data_", scen, "_2015_2040.csv"))
  
  d_soil_scen = d_scen_2015_2040 %>%
    distinct(Lon,Lat,sand_fraction,silt_fraction, clay_fraction, bulkdensity_soil, ph_soil, soilcarbon)
  
  d_scen = as.data.frame(cbind(rownames(plot_data_all[plot_data_all$scenario == long_names_scenarios(scen),]), plot_data_all[plot_data_all$scenario == long_names_scenarios(scen),"Cluster"]))
  colnames(d_scen) = c("Scenario", "Cluster")
  
  # Transform Location to Lon and Lat values
  numeric_values <- as.numeric(unlist(str_extract_all(d_scen$Scenario, "\\d+")))
  
  # Remove the ones from the vector
  numeric_values <- numeric_values[!(numeric_values < 8)]
  
  d_scen$Lon = numeric_values[seq(1, length(numeric_values), by = 4)] + 0.01 * numeric_values[seq(2, length(numeric_values), by = 4)]
  d_scen$Lat = numeric_values[seq(3, length(numeric_values), by = 4)] + 0.01 * numeric_values[seq(4, length(numeric_values), by = 4)]
  
  d_scen[which(!paste(d_scen$Lon, d_scen$Lat) %in% paste(d_soil_scen$Lon, d_soil_scen$Lat)),]$Lon = -d_scen[which(!paste(d_scen$Lon, d_scen$Lat) %in% paste(d_soil_scen$Lon, d_soil_scen$Lat)),]$Lon
  print(nrow(anti_join(d_scen, d_soil_scen, by = c("Lon", "Lat"))))
  
  # Merge data sets
  d_scen_all = inner_join(d_scen, d_soil_scen, by = c("Lon", "Lat")) %>%
    mutate(Scenario = long_names_scenarios(scen))
  
  assign(paste0("d_", scen, "_all"), d_scen_all)
}

# Merge data sets for all scenarios
d_all = rbind(d_picontrol_all, d_ssp126_all, d_ssp370_all, d_ssp585_all) %>%
  mutate(Cluster = factor(Cluster, levels = c("1", "2", "3", "4"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")),
         Cluster = fct_rev(Cluster))

d_all_long_fraction = melt(setDT(d_all[,c(1:7)]), id.vars = c("Scenario", "Cluster", "Lon", "Lat"), variable.name = "property") 
d_all_long_fraction$property = gsub("_", " ", d_all_long_fraction$property)
d_all_long_fraction$property = gsub("sand", "Sand ", d_all_long_fraction$property)
d_all_long_fraction$property = gsub("clay", "Clay ", d_all_long_fraction$property)
d_all_long_fraction$property = gsub("silt", "Silt ", d_all_long_fraction$property)


d_all_long_fraction <- d_all_long_fraction %>%
  mutate(Cluster = factor(Cluster, levels = c("1", "2", "3", "4"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")),
         Cluster = fct_rev(Cluster))

ggplot(d_all_long_fraction, aes(x = value, y = Cluster, fill = Cluster)) + 
  geom_density_ridges(aes(height = after_stat(density)), stat = "density", scale = 0.7) + 
  theme_bw() +  facet_grid(Scenario~property) +
  ylab("Scenario") + xlab("Fraction") + 
  #xlim(0, 100) +
  scale_fill_manual(
    name = "Cluster", 
    values = c("Cluster 1" = "#F8766D", "Cluster 2" = "#7CAE00", "Cluster 3" = "#00BFC4", "Cluster 4" = "#C77CFF"),
    guide = "none"  # Remove legend
  ) +
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  ggtitle("Soil composition per scenario and cluster")

ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/soil_composition.pdf"), width = 10, height = 8)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/soil_composition.png"), width = 10, height = 8)

ggplot(d_all, aes(x = bulkdensity_soil, y = Cluster, fill = Cluster)) + 
  geom_density_ridges(aes(height = after_stat(density)), stat = "density", scale = 1) + 
  theme_bw() + facet_grid(~Scenario) +
  ylab("Scenario") + xlab("Value") + 
  #xlim(0, 100) +
  scale_fill_manual(
    name = "Cluster", 
    values = c("Cluster 1" = "#F8766D", "Cluster 2" = "#7CAE00", "Cluster 3" = "#00BFC4", "Cluster 4" = "#C77CFF"),
    guide = "none"  # Remove legend
  ) +
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  ggtitle("Bulk density per scenario and cluster")

ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/bulk_density.pdf"), width = 10, height = 8)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/bulk_density.png"), width = 10, height = 8)

ggplot(d_all, aes(x = ph_soil, y = Cluster, fill = Cluster)) + 
  geom_density_ridges(aes(height = after_stat(density)), stat = "density", scale = 1) + 
  theme_bw() + facet_grid(~Scenario) +
  ylab("Scenario") + xlab("Value") + 
  #xlim(0, 100) +
  scale_fill_manual(
    name = "Cluster", 
    values = c("Cluster 1" = "#F8766D", "Cluster 2" = "#7CAE00", "Cluster 3" = "#00BFC4", "Cluster 4" = "#C77CFF"),
    guide = "none"  # Remove legend
  ) +
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  ggtitle("pH in water per scenario and cluster")

ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/ph_soil.pdf"), width = 10, height = 8)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/ph_soil.png"), width = 10, height = 8)

ggplot(d_all, aes(x = soilcarbon, y = Cluster, fill = Cluster)) + 
  geom_density_ridges(aes(height = after_stat(density)), stat = "density", scale = 1) + 
  theme_bw() + facet_grid(~Scenario) +
  ylab("Scenario") + xlab("Value") + 
  #xlim(0, 100) +
  scale_fill_manual(
    name = "Cluster", 
    values = c("Cluster 1" = "#F8766D", "Cluster 2" = "#7CAE00", "Cluster 3" = "#00BFC4", "Cluster 4" = "#C77CFF"),
    guide = "none"  # Remove legend
  ) +
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  ggtitle("Organic carbon content per scenario and cluster")

ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/soilcarbon.pdf"), width = 10, height = 5)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/soilcarbon.png"), width = 10, height = 5)
