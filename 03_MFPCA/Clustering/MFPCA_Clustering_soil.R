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
d_all = rbind(d_picontrol_all, d_ssp126_all, d_ssp370_all, d_ssp585_all)

d_all_long_fraction = melt(setDT(d_all[,c(1:7)]), id.vars = c("Scenario", "Cluster", "Lon", "Lat"), variable.name = "property") 
d_all_long_fraction$property = gsub("_", " ", d_all_long_fraction$property)
d_all_long_fraction$property = gsub("sand", "Sand ", d_all_long_fraction$property)
d_all_long_fraction$property = gsub("clay", "Clay ", d_all_long_fraction$property)
d_all_long_fraction$property = gsub("silt", "Silt ", d_all_long_fraction$property)

g_fraction = ggplot(d_all_long_fraction, aes(x=Scenario, y=value, fill=Cluster)) + 
  geom_boxplot() + facet_wrap(~property) + theme_bw() + ylab("Fraction") +
  theme(text = element_text(size = 10), plot.title = element_text(size = 15, face = "bold",hjust = 0.5), axis.text.x = element_text(angle = 30,hjust=1)) +
  ggtitle("Soil properties per scenario and cluster")


g_bulk = ggplot(d_all, aes(x=Scenario, y=bulkdensity_soil, fill=Cluster)) + 
  geom_boxplot() + theme_bw() + ylab("Value") + xlab("") +
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust=1),
    legend.position = "none"  # This line removes the legend
  ) +
  ggtitle("Bulk density")


g_ph = ggplot(d_all, aes(x=Scenario, y=ph_soil, fill=Cluster)) + 
  geom_boxplot() + theme_bw() + ylab("") + 
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust=1),
    legend.position = "none"  # This line removes the legend
  ) +
  ggtitle("PH value")

g_carbon = ggplot(d_all, aes(x=Scenario, y=soilcarbon, fill=Cluster)) + 
  geom_boxplot() + theme_bw() + ylab("") + xlab("") +
  theme(text = element_text(size = 10), plot.title = element_text(size = 15, face = "bold",hjust = 0.5), axis.text.x = element_text(angle = 30,hjust=1)) +
  ggtitle("Soil carbon")

plot_grid(
  g_fraction, 
  plot_grid(g_bulk, g_ph, g_carbon, ncol = 3, rel_widths = c(0.9, 0.9, 1)), 
  ncol = 1, 
  rel_heights = c(1, 1)
)

ggsave(paste0("Scripts/Plots/MFPCA/Clusters/pdf/soil_props_",pid,".pdf"), width = 10, height = 5)
ggsave(paste0("Scripts/Plots/MFPCA/Clusters/png/soil_props_",pid,".png"), width = 10, height = 5)
