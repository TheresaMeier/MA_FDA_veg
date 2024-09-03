################################################################################
############################ Master's Thesis ###################################
################################################################################

########################## Clustering: Climate data ############################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/01_Description/utils.R")

## Load libraries
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(cowplot)
library(ggplot2)

scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")

## Load data
pid = 1

# Get locations and cluster assignments
d_loc = read.table("Scripts/Plots/MFPCA/plot_data/locs_disturbed.txt")

for (scen in scenarios){
  d_scen_2015_2040 = fread(paste0("Data/data_", scen, "_2015_2040.csv"))
  
  d_climate_scen = d_scen_2015_2040 %>%
    filter(PID == pid) %>%
    distinct(Lon,Lat,Year,tas_yearlymin,tas_yearlymeam, tas_yearlymax,pr_yearlysum) %>%
    mutate(tas_yearlymin = tas_yearlymin - 273.15,
           tas_yearlymeam = tas_yearlymeam - 273.15,
           tas_yearlymax = tas_yearlymax - 273.15,
           Scenario = long_names_scenarios(scen))
  
  d_scen = d_loc[d_loc$Scenario == long_names_scenarios(scen),]
  
  # Merge data sets
  counts <- d_climate_scen %>%
    group_by(Lon, Lat) %>%
    summarise(count = n()) %>%
    ungroup()
  
  d_scen =  d_scen %>%
    inner_join(counts, by = c("Lon", "Lat")) %>%
    uncount(count)
  
  d_climate_scen$Cluster = d_scen$Cluster
  assign(paste0("d_", scen, "_all"), d_climate_scen)
}

# Merge data sets for all scenarios
d_all = rbind(d_picontrol_all, d_ssp126_all, d_ssp370_all, d_ssp585_all) %>%
  mutate(Cluster = as.factor(Cluster))

################################ Clusters only  ################################

# calculate mean values per Cluster
d_climate_mean_clusterOnly = d_all %>%
  group_by(Cluster, Year) %>%
  summarize(mean_tas_min = mean(tas_yearlymin, na.rm = T),
            mean_tas_mean = mean(tas_yearlymeam, na.rm = T),
            mean_tas_max = mean(tas_yearlymax, na.rm = T),
            mean_precip = mean(pr_yearlysum, na.rm = T)) %>%
  ungroup %>%
  mutate(Cluster = as.factor(Cluster))

# Yearly mean temperature
ggplot(d_climate_mean_clusterOnly) + 
  geom_line(data = d_climate_mean_clusterOnly, linewidth = 1, 
            aes(x = Year, y = mean_tas_mean, col = Cluster)) +
  geom_hline(yintercept = 0, color = "black", size = 0.25, alpha = .95) +
  scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 20)) +
  scale_y_continuous(name = "Yearly mean temperature in °C", expand = c(0,0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Annual mean temperature") +
  theme_bw() + theme(
    text = element_text(size = 22), 
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust=1))

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/temp_mean.pdf", width = 8, height = 6)
ggsave("Scripts/Plots/MFPCA/Clusters/png/temp_mean.png", width = 8, height = 6)

# Yearly min temperature
ggplot(d_climate_mean_clusterOnly) + 
  geom_line(data = d_climate_mean_clusterOnly, linewidth = 1, 
            aes(x = Year, y = mean_tas_min, col = Cluster)) +
  scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 20)) +
  scale_y_continuous(name = "Yearly minimum temperature in °C", expand = c(0,0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Annual minimum temperature") +
  theme_bw() + theme(
    text = element_text(size = 22), 
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust=1))

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/temp_min.pdf", width = 8, height = 6)
ggsave("Scripts/Plots/MFPCA/Clusters/png/temp_min.png", width = 8, height = 6)

# Yearly max temperature
ggplot(d_climate_mean_clusterOnly) + 
  geom_line(data = d_climate_mean_clusterOnly, linewidth = 1,
            aes(x = Year, y = mean_tas_max, col = Cluster)) +
  scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 20)) +
  scale_y_continuous(name = "Yearly maximum temperature in °C", expand = c(0,0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Annual maximum temperature") +
  theme_bw() + theme(
    text = element_text(size = 22), 
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust=1))

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/temp_max.pdf", width = 8, height = 6)
ggsave("Scripts/Plots/MFPCA/Clusters/png/temp_max.png", width = 8, height = 6)

# Precipitation

ggplot(d_climate_mean_clusterOnly) + 
  geom_line(data = d_climate_mean_clusterOnly, linewidth = 1, 
            aes(x = Year, y = mean_precip, col = Cluster)) +
  scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 20)) +
  scale_y_continuous(name = bquote("Yearly precipitation in kg/m"^2), expand = c(0, 0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Annual precipitation") +
  theme_bw() + theme(
    text = element_text(size = 22), 
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust=1))

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/precip.pdf", width = 8, height = 6)
ggsave("Scripts/Plots/MFPCA/Clusters/png/precip.png", width = 8, height = 6)

