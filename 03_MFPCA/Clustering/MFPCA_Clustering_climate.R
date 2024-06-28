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

############################ Cluster and Scenarios #############################

# calculate mean values per Cluster and Scenario
d_climate_mean_cluster = d_all %>%
  group_by(Cluster, Scenario, Year) %>%
  summarize(mean_tas_min = mean(tas_yearlymin, na.rm = T),
            mean_tas_mean = mean(tas_yearlymeam, na.rm = T),
            mean_tas_max = mean(tas_yearlymax, na.rm = T),
            mean_precip = mean(pr_yearlysum, na.rm = T)) %>%
  ungroup %>%
  mutate(Cluster = as.factor(Cluster))

# Yearly mean temperature
ggplot(d_climate_mean_cluster) + 
  geom_line(data = d_climate_mean_cluster, linewidth = 1,
            aes(x = Year, y = mean_tas_mean, col = Cluster)) +
  # geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_tas_mean, group = Cluster), 
  #             method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
  # geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_tas_mean, color = Cluster),
  #             method = "loess", se = FALSE, linewidth = 1) +  
  facet_wrap(~Scenario) +
  geom_hline(yintercept = 0, color = "black", size = 0.25, alpha = .95) +
  # geom_hline(yintercept = 5, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
  # geom_hline(yintercept = -5, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
  scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 20)) +
  scale_y_continuous(name = "Yearly mean temperature in °C", expand = c(0,0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Yearly mean temperature averaged over disturbed grid cells") +
  theme_bw() + theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust=1)
  )

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/temp_mean_scen.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/temp_mean_scen.png", width = 8, height = 5)

# # Cut y-axis
# ggplot(d_all) + 
#   geom_line(data = d_all, linewidth = 0.05, alpha = .25,
#             aes(x = Year, y = tas_yearlymeam, col = Cluster, group = interaction(Lon, Lat))) +
#   geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_tas_mean, group = Cluster), 
#               method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
#   geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_tas_mean, color = Cluster), 
#               method = "loess", se = FALSE, linewidth = 1) +  facet_wrap(~Scenario) +
#   geom_hline(yintercept = 0, color = "black", size = 0.25, alpha = .95) +
#   geom_hline(yintercept = 5, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
#   geom_hline(yintercept = -5, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
#   scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 10)) +
#   scale_y_continuous(name = "Yearly mean temperature in °C", expand = c(0,0), limits = c(-2,6)) +
#   scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
#   ggtitle("Yearly mean temperature for each disturbed grid cell and scenario") +
#   theme_bw() + theme(
#     text = element_text(size = 10), 
#     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
#   )
# 
# ggsave("Scripts/Plots/MFPCA/Clusters/pdf/temp_mean_scen_cut.pdf", width = 8, height = 5)
# ggsave("Scripts/Plots/MFPCA/Clusters/png/temp_mean_scen_cut.png", width = 8, height = 5)


# Yearly min temperature
ggplot(d_climate_mean_cluster) + 
  geom_line(data = d_climate_mean_cluster, linewidth = 1,
            aes(x = Year, y = mean_tas_min, col = Cluster)) +
  # geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_tas_min, group = Cluster), 
  #             method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
  # geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_tas_min, color = Cluster), 
  #             method = "loess", se = FALSE, linewidth = 1) +  
  facet_wrap(~Scenario) +
  # geom_hline(yintercept = -20, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
  # geom_hline(yintercept = -40, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
  scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 20)) +
  scale_y_continuous(name = "Yearly minimum temperature in °C", expand = c(0,0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Daily minimum temperature per year averaged over disturbed grid cells") +
  theme_bw() + theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust=1)
  )


ggsave("Scripts/Plots/MFPCA/Clusters/pdf/temp_min_scen.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/temp_min_scen.png", width = 8, height = 5)

# # Cut y-axis
# ggplot(d_all) + 
#   geom_line(data = d_all, linewidth = 0.05, alpha = .25,
#             aes(x = Year, y = tas_yearlymin, col = Cluster, group = interaction(Lon, Lat))) +
#   geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_tas_min, group = Cluster), 
#               method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
#   geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_tas_min, color = Cluster), 
#               method = "loess", se = FALSE, linewidth = 1) +  facet_wrap(~Scenario) +
#   geom_hline(yintercept = -20, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
#   geom_hline(yintercept = -40, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
#   scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 10)) +
#   scale_y_continuous(name = "Yearly minimum temperature in °C", expand = c(0,0), limits = c(-40,-20)) +
#   scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
#   ggtitle("Daily minimum temperature per year for each disturbed grid cell and scenario") +
#   theme_bw() + theme(
#     text = element_text(size = 10), 
#     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
#   )
# 
# ggsave("Scripts/Plots/MFPCA/Clusters/pdf/temp_min_scen_cut.pdf", width = 8, height = 5)
# ggsave("Scripts/Plots/MFPCA/Clusters/png/temp_min_scen_cut.png", width = 8, height = 5)

# Yearly max temperature
ggplot(d_climate_mean_cluster) + 
  geom_line(data = d_climate_mean_cluster, linewidth = 1,
            aes(x = Year, y = mean_tas_max, col = Cluster)) +
  # geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_tas_max, group = Cluster), 
  #             method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
  # geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_tas_max, color = Cluster), 
  #             method = "loess", se = FALSE, linewidth = 1) +  
  facet_wrap(~Scenario) +
  # geom_hline(yintercept = 20, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
  # geom_hline(yintercept = 30, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
  scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 20)) +
  scale_y_continuous(name = "Yearly maximum temperature in °C", expand = c(0,0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Daily maximum temperature per year averaged over disturbed grid cells") +
  theme_bw() + theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust=1)
  )

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/temp_max_scen.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/temp_max_scen.png", width = 8, height = 5)

# # Cut y-axis
# ggplot(d_all) + 
#   geom_line(data = d_all, linewidth = 0.05, alpha = .25,
#             aes(x = Year, y = tas_yearlymax, col = Cluster, group = interaction(Lon, Lat))) +
#   geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_tas_max, group = Cluster), 
#               method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
#   geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_tas_max, color = Cluster), 
#               method = "loess", se = FALSE, linewidth = 1) +  facet_wrap(~Scenario) +
#   geom_hline(yintercept = 20, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
#   geom_hline(yintercept = 30, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
#   scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 10)) +
#   scale_y_continuous(name = "Yearly maximum temperature in °C", expand = c(0,0), limits = c(21,27.5)) +
#   scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
#   ggtitle("Daily maximum temperature per year for each disturbed grid cell and scenario") +
#   theme_bw() + theme(
#     text = element_text(size = 10), 
#     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
#   )
# 
# ggsave("Scripts/Plots/MFPCA/Clusters/pdf/temp_max_scen_cut.pdf", width = 8, height = 5)
# ggsave("Scripts/Plots/MFPCA/Clusters/png/temp_max_scen_cut.png", width = 8, height = 5)
# 

# Precipitation

ggplot(d_climate_mean_cluster) + 
  geom_line(data = d_climate_mean_cluster, linewidth = 1,
            aes(x = Year, y = mean_precip, col = Cluster)) +
  # geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_precip, group = Cluster), 
  #             method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
  # geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_precip, color = Cluster), 
  #             method = "loess", se = FALSE, linewidth = 1) +  
  facet_wrap(~Scenario) +
  scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 20)) +
  scale_y_continuous(name = bquote("Yearly precipitation in kg/m"^2), expand = c(0, 0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Yearly precipitation averaged over disturbed grid cells") +
  theme_bw() + theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust=1)
  )

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/precip_scen.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/precip_scen.png", width = 8, height = 5)

# # Cut y-axis
# ggplot(d_all) + 
#   geom_line(data = d_all, linewidth = 0.05, alpha = .25,
#             aes(x = Year, y = pr_yearlysum, col = Cluster, group = interaction(Lon, Lat))) +
#   geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_precip, group = Cluster), 
#               method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
#   geom_smooth(data = d_climate_mean_cluster, aes(x = Year, y = mean_precip, color = Cluster), 
#               method = "loess", se = FALSE, linewidth = 1) +  facet_wrap(~Scenario) +
#   scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 10)) +
#   scale_y_continuous(name = bquote("Yearly precipitation in kg/m"^2), expand = c(0, 0), limits = c(550,900)) +
#   scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
#   ggtitle("Yearly precipitation for each disturbed grid cell and scenario") +
#   theme_bw() + theme(
#     text = element_text(size = 10), 
#     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
#   )
# 
# ggsave("Scripts/Plots/MFPCA/Clusters/pdf/precip_scen_cut.pdf", width = 8, height = 5)
# ggsave("Scripts/Plots/MFPCA/Clusters/png/precip_scen_cut.png", width = 8, height = 5)
# 

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
  # geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_tas_mean, group = Cluster), 
  #             method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
  # geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_tas_mean, color = Cluster), 
  #             method = "loess", se = FALSE, linewidth = 1) +
  geom_hline(yintercept = 0, color = "black", size = 0.25, alpha = .95) +
  # geom_hline(yintercept = 5, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
  # geom_hline(yintercept = -5, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
  scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 20)) +
  scale_y_continuous(name = "Yearly mean temperature in °C", expand = c(0,0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Yearly mean temperature averaged over disturbed grid cells") +
  theme_bw() + theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust=1)
  )

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/temp_mean.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/temp_mean.png", width = 8, height = 5)

# # Cut y-axis
# ggplot(d_all) + 
#   geom_line(data = d_all, linewidth = 0.05, alpha = .25,
#             aes(x = Year, y = tas_yearlymeam, col = Cluster, group = interaction(Lon, Lat))) +
#   geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_tas_mean, group = Cluster), 
#               method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
#   geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_tas_mean, color = Cluster), 
#               method = "loess", se = FALSE, linewidth = 1) +
#   geom_hline(yintercept = 0, color = "black", size = 0.25, alpha = .95) +
#   geom_hline(yintercept = 5, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
#   geom_hline(yintercept = -5, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
#   scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 10)) +
#   scale_y_continuous(name = "Yearly mean temperature in °C", expand = c(0,0), limits = c(-0.5,3.5)) +
#   scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
#   ggtitle("Yearly mean temperature for each disturbed grid cell") +
#   theme_bw() + theme(
#     text = element_text(size = 10), 
#     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
#   )
# 
# ggsave("Scripts/Plots/MFPCA/Clusters/pdf/temp_mean_cut.pdf", width = 8, height = 5)
# ggsave("Scripts/Plots/MFPCA/Clusters/png/temp_mean_cut.png", width = 8, height = 5)


# Yearly min temperature
ggplot(d_climate_mean_clusterOnly) + 
  geom_line(data = d_climate_mean_clusterOnly, linewidth = 1, 
            aes(x = Year, y = mean_tas_min, col = Cluster)) +
  # geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_tas_min, group = Cluster), 
  #             method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
  # geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_tas_min, color = Cluster), 
  #             method = "loess", se = FALSE, linewidth = 1) +
  # geom_hline(yintercept = -20, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
  # geom_hline(yintercept = -40, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
  scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 20)) +
  scale_y_continuous(name = "Yearly minimum temperature in °C", expand = c(0,0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Daily minimum temperature per year averaged over disturbed grid cells") +
  theme_bw() + theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust=1)
  )


ggsave("Scripts/Plots/MFPCA/Clusters/pdf/temp_min.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/temp_min.png", width = 8, height = 5)

# # Cut y-axis
# ggplot(d_all) + 
#   geom_line(data = d_all, linewidth = 0.05, alpha = .25,
#             aes(x = Year, y = tas_yearlymin, col = Cluster, group = interaction(Lon, Lat))) +
#   geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_tas_min, group = Cluster), 
#               method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
#   geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_tas_min, color = Cluster), 
#               method = "loess", se = FALSE, linewidth = 1) +
#   geom_hline(yintercept = -20, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
#   geom_hline(yintercept = -40, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
#   scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 10)) +
#   scale_y_continuous(name = "Yearly minimum temperature in °C", expand = c(0,0), limits = c(-36,-25)) +
#   scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
#   ggtitle("Daily minimum temperature per year for each disturbed grid cell") +
#   theme_bw() + theme(
#     text = element_text(size = 10), 
#     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
#   )
# 
# ggsave("Scripts/Plots/MFPCA/Clusters/pdf/temp_min_cut.pdf", width = 8, height = 5)
# ggsave("Scripts/Plots/MFPCA/Clusters/png/temp_min_cut.png", width = 8, height = 5)

# Yearly max temperature
ggplot(d_climate_mean_clusterOnly) + 
  geom_line(data = d_climate_mean_clusterOnly, linewidth = 1,
            aes(x = Year, y = mean_tas_max, col = Cluster)) +
  # geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_tas_max, group = Cluster), 
  #             method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
  # geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_tas_max, color = Cluster), 
  #             method = "loess", se = FALSE, linewidth = 1) + 
  # geom_hline(yintercept = 20, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
  # geom_hline(yintercept = 30, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
  scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 20)) +
  scale_y_continuous(name = "Yearly maximum temperature in °C", expand = c(0,0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Daily maximum temperature per year averaged over disturbed grid cells") +
  theme_bw() + theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust=1)
  )

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/temp_max.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/temp_max.png", width = 8, height = 5)

# # Cut y-axis
# ggplot(d_all) + 
#   geom_line(data = d_all, linewidth = 0.05, alpha = .25,
#             aes(x = Year, y = tas_yearlymax, col = Cluster, group = interaction(Lon, Lat))) +
#   geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_tas_max, group = Cluster), 
#               method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
#   geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_tas_max, color = Cluster), 
#               method = "loess", se = FALSE, linewidth = 1) + 
#   geom_hline(yintercept = 20, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
#   geom_hline(yintercept = 30, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
#   scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 10)) +
#   scale_y_continuous(name = "Yearly maximum temperature in °C", expand = c(0,0), limits = c(22,26)) +
#   scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
#   ggtitle("Daily maximum temperature per year for each disturbed grid cell") +
#   theme_bw() + theme(
#     text = element_text(size = 10), 
#     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
#   )
# 
# ggsave("Scripts/Plots/MFPCA/Clusters/pdf/temp_max_cut.pdf", width = 8, height = 5)
# ggsave("Scripts/Plots/MFPCA/Clusters/png/temp_max_cut.png", width = 8, height = 5)

# Precipitation

ggplot(d_climate_mean_clusterOnly) + 
  geom_line(data = d_climate_mean_clusterOnly, linewidth = 1, 
            aes(x = Year, y = mean_precip, col = Cluster)) +
  # geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_precip, group = Cluster), 
  #             method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
  # geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_precip, color = Cluster), 
  #             method = "loess", se = FALSE, linewidth = 1) +
  scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 20)) +
  scale_y_continuous(name = bquote("Yearly precipitation in kg/m"^2), expand = c(0, 0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Yearly precipitation averaged over disturbed grid cells") +
  theme_bw() + theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust=1)
  )

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/precip.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/precip.png", width = 8, height = 5)

# # Cut y-axis
# ggplot(d_all) + 
#   geom_line(data = d_all, linewidth = 0.05, alpha = .25,
#             aes(x = Year, y = pr_yearlysum, col = Cluster, group = interaction(Lon, Lat))) +
#   geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_precip, group = Cluster), 
#               method = "loess", se = FALSE, linewidth = 2, color = "black", show.legend = FALSE) +
#   geom_smooth(data = d_climate_mean_clusterOnly, aes(x = Year, y = mean_precip, color = Cluster), 
#               method = "loess", se = FALSE, linewidth = 1) +
#   scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 20)) +
#   scale_y_continuous(name = bquote("Yearly precipitation in kg/m"^2), expand = c(0, 0), limits = c(600, 850)) +
#   scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
#   ggtitle("Yearly precipitation for each disturbed grid cell") +
#   theme_bw() + theme(
#     text = element_text(size = 10), 
#     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
#   )
# 
# ggsave("Scripts/Plots/MFPCA/Clusters/pdf/precip_cut.pdf", width = 8, height = 5)
# ggsave("Scripts/Plots/MFPCA/Clusters/png/precip_cut.png", width = 8, height = 5)


