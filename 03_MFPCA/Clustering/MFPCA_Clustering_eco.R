################################################################################
############################ Master's Thesis ###################################
################################################################################

######################### Clustering: Ecological data ##########################

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
  
  d_eco_scen = d_scen_2015_2040 %>%
    filter(PID == pid) %>%
    distinct(Lon,Lat,Year,age,PFT,relative,initial_recruitment,
             recruitment_ten_years, previous_state, time_since_dist,
             Nuptake_total,Nuptake) %>%
    mutate(Scenario = long_names_scenarios(scen),
           PFT = long_names_pfts(tolower(PFT)))
  
  d_scen = d_loc[d_loc$Scenario == long_names_scenarios(scen),]
  
  # Merge data sets
  counts <- d_eco_scen %>%
    group_by(Lon, Lat) %>%
    summarise(count = n()) %>%
    ungroup()
  
  d_scen =  d_scen %>%
    inner_join(counts, by = c("Lon", "Lat")) %>%
    uncount(count)
  
  d_eco_scen$Cluster = d_scen$Cluster
  assign(paste0("d_", scen, "_all"), d_eco_scen)
}

# Merge data sets for all scenarios
d_all = rbind(d_picontrol_all, d_ssp126_all, d_ssp370_all, d_ssp585_all) %>%
  mutate(Cluster = as.factor(Cluster))


############################ Cluster and Scenarios #############################

# calculate mean values per Cluster and Scenario
d_eco_mean_cluster = d_all %>%
  group_by(Scenario, Cluster, PFT, Year) %>%
  summarize(mean_nuptake = mean(Nuptake, na.rm = T),
            mean_nuptake_total = mean(Nuptake_total, na.rm = T)) %>%
  ungroup %>%
  mutate(Cluster = as.factor(Cluster))

# Nuptake
ggplot(d_all) + 
  geom_line(data = d_all, linewidth = .05, alpha = .25,
            aes(x = Year, y = Nuptake, color = Cluster, group = interaction(Lon, Lat, Cluster,PFT))) +
  geom_smooth(data = d_eco_mean_cluster, aes(x = Year, y = mean_nuptake, color = Cluster, group = Cluster),
              method = "loess", se = FALSE, linewidth = 2.5, color = "black", show.legend = FALSE) +
  geom_smooth(data = d_eco_mean_cluster, aes(x = Year, y = mean_nuptake, color = Cluster, group = Cluster), 
              method = "loess", se = FALSE, linewidth = 1.5) +  facet_grid(Scenario~PFT) +
  scale_x_continuous(name = "Year after disturbance", expand = c(0,0), 
                     breaks = seq(2020, 2120, by = 20)) +
  scale_y_continuous(name = "Nitrogen uptake per PFT on the grid cell level", 
                     expand = c(0,0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Nitrogen uptake per PFT for each disturbed grid cell and scenario") +
  theme_bw() + 
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45)
  )
ggsave("Scripts/Plots/MFPCA/Clusters/pdf/nuptake_scen.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/nuptake_scen.png", width = 8, height = 5)

# Cut y-axis
ggplot(d_all) + 
  geom_line(data = d_all, linewidth = .05, alpha = .25,
            aes(x = Year, y = Nuptake, color = Cluster, group = interaction(Lon, Lat, Cluster,PFT))) +
  geom_smooth(data = d_eco_mean_cluster, aes(x = Year, y = mean_nuptake, color = Cluster, group = Cluster),
              method = "loess", se = FALSE, linewidth = 2.5, color = "black", show.legend = FALSE) +
  geom_smooth(data = d_eco_mean_cluster, aes(x = Year, y = mean_nuptake, color = Cluster, group = Cluster), 
              method = "loess", se = FALSE, linewidth = 1.5) +  facet_grid(Scenario~PFT) +
  scale_x_continuous(name = "Year after disturbance", expand = c(0,0), 
                     breaks = seq(2020, 2120, by = 20)) +
  scale_y_continuous(name = "Nitrogen uptake per PFT on the grid cell level", 
                     expand = c(0,0), limits = c(0,30)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Nitrogen uptake per PFT for each disturbed grid cell and scenario") +
  theme_bw() + 
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45)
    )

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/nuptake_scen_cut.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/nuptake_scen_cut.png", width = 8, height = 5)

# Nuptake total
ggplot(d_all) + 
  geom_line(data = d_all, linewidth = .05, alpha = .25,
            aes(x = Year, y = Nuptake_total, color = Cluster, group = interaction(Lon, Lat, Cluster))) +
  geom_smooth(data = d_eco_mean_cluster, aes(x = Year, y = mean_nuptake_total, color = Cluster, group = Cluster),
              method = "loess", se = FALSE, linewidth = 2.5, color = "black", show.legend = FALSE) +
  geom_smooth(data = d_eco_mean_cluster, aes(x = Year, y = mean_nuptake_total, color = Cluster, group = Cluster), 
              method = "loess", se = FALSE, linewidth = 1.5) +  facet_grid(~Scenario) +
  scale_x_continuous(name = "Year after disturbance", expand = c(0,0), 
                     breaks = seq(2020, 2120, by = 20)) +
  scale_y_continuous(name = "Total nitrogen uptake of the grid cell", 
                     expand = c(0,0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Total nitrogen uptake of the grid cell") +
  theme_bw() + 
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45)
  )

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/nuptake_total_scen.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/nuptake_total_scen.png", width = 8, height = 5)

# Cut y-axis
ggplot(d_all) + 
  geom_line(data = d_all, linewidth = .05, alpha = .25,
            aes(x = Year, y = Nuptake_total, color = Cluster, group = interaction(Lon, Lat, Cluster))) +
  geom_smooth(data = d_eco_mean_cluster, aes(x = Year, y = mean_nuptake_total, color = Cluster, group = Cluster),
              method = "loess", se = FALSE, linewidth = 2.5, color = "black", show.legend = FALSE) +
  geom_smooth(data = d_eco_mean_cluster, aes(x = Year, y = mean_nuptake_total, color = Cluster, group = Cluster), 
              method = "loess", se = FALSE, linewidth = 1.5) +  facet_grid(~Scenario) +
  scale_x_continuous(name = "Year after disturbance", expand = c(0,0), 
                     breaks = seq(2020, 2120, by = 20)) +
  scale_y_continuous(name = "Total nitrogen uptake of the grid cell", 
                     expand = c(0,0), limits = c(25,50)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Total nitrogen uptake of the grid cell") +
  theme_bw() + 
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45)
  )

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/nuptake_total_scen_cut.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/nuptake_total_scen_cut.png", width = 8, height = 5)

################################ Per cluster ###################################
# calculate mean values per Cluster
d_eco_mean = d_all %>%
  group_by(Cluster, Year) %>%
  summarize(mean_nuptake = mean(Nuptake, na.rm = T),
            mean_nuptake_total = mean(Nuptake_total, na.rm = T)) %>%
  ungroup %>%
  mutate(Cluster = as.factor(Cluster))

# Nuptake
ggplot(d_all) + 
  geom_line(data = d_all, linewidth = .05, alpha = .25,
            aes(x = Year, y = Nuptake, color = Cluster, group = interaction(Lon, Lat, Cluster))) +
  geom_smooth(data = d_eco_mean, aes(x = Year, y = mean_nuptake, color = Cluster, group = Cluster),
              method = "loess", se = FALSE, linewidth = 2.5, color = "black", show.legend = FALSE) +
  geom_smooth(data = d_eco_mean, aes(x = Year, y = mean_nuptake, color = Cluster, group = Cluster), 
              method = "loess", se = FALSE, linewidth = 1.5) +
  scale_x_continuous(name = "Year after disturbance", expand = c(0,0), 
                     breaks = seq(2020, 2120, by = 20)) +
  scale_y_continuous(name = "Nitrogen uptake on the grid cell level", 
                     expand = c(0,0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Nitrogen uptake per cluster for each disturbed grid cell") +
  theme_bw() + 
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45)
  )

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/nuptake.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/nuptake.png", width = 8, height = 5)

# Cut y-axis

ggplot(d_all) + 
  geom_line(data = d_all, linewidth = .05, alpha = .25,
            aes(x = Year, y = Nuptake, color = Cluster, group = interaction(Lon, Lat, Cluster))) +
  geom_smooth(data = d_eco_mean, aes(x = Year, y = mean_nuptake, color = Cluster, group = Cluster),
              method = "loess", se = FALSE, linewidth = 2.5, color = "black", show.legend = FALSE) +
  geom_smooth(data = d_eco_mean, aes(x = Year, y = mean_nuptake, color = Cluster, group = Cluster), 
              method = "loess", se = FALSE, linewidth = 1.5) +
  scale_x_continuous(name = "Year after disturbance", expand = c(0,0), 
                     breaks = seq(2020, 2120, by = 20)) +
  scale_y_continuous(name = "Nitrogen uptake on the grid cell level", 
                     expand = c(0,0), limits = c(4,10)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Nitrogen uptake per cluster for each disturbed grid cell") +
  theme_bw() + 
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45)
  )

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/nuptake_cut.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/nuptake_cut.png", width = 8, height = 5)

# Nuptake total
ggplot(d_all) + 
  geom_line(data = d_all, linewidth = .05, alpha = .25,
            aes(x = Year, y = Nuptake_total, color = Cluster, group = interaction(Lon, Lat, Cluster))) +
  geom_smooth(data = d_eco_mean, aes(x = Year, y = mean_nuptake_total, color = Cluster, group = Cluster),
              method = "loess", se = FALSE, linewidth = 2.5, color = "black", show.legend = FALSE) +
  geom_smooth(data = d_eco_mean, aes(x = Year, y = mean_nuptake_total, color = Cluster, group = Cluster), 
              method = "loess", se = FALSE, linewidth = 1.5) +
  scale_x_continuous(name = "Year after disturbance", expand = c(0,0), 
                     breaks = seq(2020, 2120, by = 20)) +
  scale_y_continuous(name = "Total nitrogen uptake of the grid cell", 
                     expand = c(0,0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Total nitrogen uptake per cluster of the grid cell") +
  theme_bw() + 
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45)
  )

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/nuptake_total.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/nuptake_total.png", width = 8, height = 5)

# Cut y-axis
ggplot(d_all) + 
  geom_line(data = d_all, linewidth = .05, alpha = .25,
            aes(x = Year, y = Nuptake_total, color = Cluster, group = interaction(Lon, Lat, Cluster))) +
  geom_smooth(data = d_eco_mean, aes(x = Year, y = mean_nuptake_total, color = Cluster, group = Cluster),
              method = "loess", se = FALSE, linewidth = 2.5, color = "black", show.legend = FALSE) +
  geom_smooth(data = d_eco_mean, aes(x = Year, y = mean_nuptake_total, color = Cluster, group = Cluster), 
              method = "loess", se = FALSE, linewidth = 1.5) + 
  scale_x_continuous(name = "Year after disturbance", expand = c(0,0), 
                     breaks = seq(2020, 2120, by = 20)) +
  scale_y_continuous(name = "Total nitrogen uptake of the grid cell", 
                     expand = c(0,0), limits = c(30,45)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Total nitrogen uptake per cluster of the grid cell") +
  theme_bw() + 
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45)
  )

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/nuptake_total_cut.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/nuptake_total_cut.png", width = 8, height = 5)

###

d_recruit = d_all %>%
  distinct(Lon,Lat,PFT,initial_recruitment, recruitment_ten_years, previous_state, time_since_dist, Scenario, Cluster) %>%
  mutate(PFT = factor(PFT, levels = c("Needleleaf evergreen", "Pioneering broadleaf", "Conifers (other)", "Temperate broadleaf", "Tundra")))


############################# Initial recruitment ##############################
ggplot(d_recruit, aes(x=Scenario, y=initial_recruitment, fill = Cluster)) + 
  geom_boxplot() + theme_bw() + ylab("Number of new seedlings") + facet_grid(~PFT) + #ylim(0,100) +
  scale_fill_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Number of new seedlings per PFT right after disturbance ")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/initial_recruitment_scen_PFT.pdf", width = 10, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/initial_recruitment_scen_PFT.png", width = 10, height = 5)

# Cut y-axis
ggplot(d_recruit, aes(x=Scenario, y=initial_recruitment, fill = Cluster)) + 
  geom_boxplot() + theme_bw() + ylab("Number of new seedlings") + facet_grid(~PFT) + ylim(0,100) +
  scale_fill_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Number of new seedlings per PFT right after disturbance ")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/initial_recruitment_scen_PFT_cut.pdf", width = 10, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/initial_recruitment_scen_PFT_cut.png", width = 10, height = 5)

# For all scenarios and PFTs together
ggplot(d_recruit, aes(x = Cluster, y=initial_recruitment, fill = Cluster)) + 
  geom_boxplot() + theme_bw() + ylab("Number of new seedlings") + #ylim(0,100) +
  scale_fill_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  ggtitle("Number of new seedlings per cluster right after disturbance ")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/initial_recruitment.pdf", width = 10, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/initial_recruitment.png", width = 10, height = 5)

# Cut y-axis
ggplot(d_recruit, aes(x = Cluster, y=initial_recruitment, fill = Cluster)) + 
  geom_boxplot() + theme_bw() + ylab("Number of new seedlings") + ylim(0,100) +
  scale_fill_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  ggtitle("Number of new seedlings per cluster right after disturbance ")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/initial_recruitment_cut.pdf", width = 10, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/initial_recruitment_cut.png", width = 10, height = 5)


############################# Recruitment 10 years #############################
ggplot(d_recruit, aes(x=Scenario, y=recruitment_ten_years, fill = Cluster)) + 
  geom_boxplot() + theme_bw() + ylab("Number of new seedlings") + facet_grid(~PFT) + #ylim(0,500) +
  scale_fill_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Number of new seedling per PFT in the ten years after disturbance (summed up)")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/recruitment_10_scen_PFT.pdf", width = 10, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/recruitment_10_scen_PFT.png", width = 10, height = 5)

# Cut y-axis
ggplot(d_recruit, aes(x=Scenario, y=recruitment_ten_years, fill = Cluster)) + 
  geom_boxplot() + theme_bw() + ylab("Number of new seedlings") + facet_grid(~PFT) + ylim(0,500) +
  scale_fill_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Number of new seedling per PFT in the ten years after disturbance (summed up)")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/recruitment_10_scen_PFT_cut.pdf", width = 10, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/recruitment_10_scen_PFT_cut.png", width = 10, height = 5)

# For all scenarios and PFTs together
ggplot(d_recruit, aes(x = Cluster, y=recruitment_ten_years, fill = Cluster)) + 
  geom_boxplot() + theme_bw() + ylab("Number of new seedlings") + #ylim(0,30) +
  scale_fill_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  ggtitle("Number of new seedling per cluster in the ten years after disturbance (summed up)")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/recruitment_10.pdf", width = 10, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/recruitment_10.png", width = 10, height = 5)

# Cut y-axis
ggplot(d_recruit, aes(x = Cluster, y=recruitment_ten_years, fill = Cluster)) + 
  geom_boxplot() + theme_bw() + ylab("Number of new seedlings") + ylim(0,30) +
  scale_fill_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  ggtitle("Number of new seedling per cluster in the ten years after disturbance (summed up)")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/recruitment_10_cut.pdf", width = 10, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/recruitment_10_cut.png", width = 10, height = 5)


############################## previous state ##################################

ggplot(d_recruit, aes(x=Scenario, y=previous_state, fill = Cluster)) + 
  geom_boxplot() + theme_bw() + facet_grid(~PFT) +
  ylab(bquote("Aboveground carbon in kg/m"^2)) + #ylim(0,10) +
  scale_fill_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Vegetation composition before the disturbance")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/previous_state_scen_PFT.pdf", width = 10, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/previous_state_scen_PFT.png", width = 10, height = 5)

# Cut y-axis
ggplot(d_recruit, aes(x=Scenario, y=previous_state, fill = Cluster)) + 
  geom_boxplot() + theme_bw() + facet_grid(~PFT) +
  ylab(bquote("Aboveground carbon in kg/m"^2)) + ylim(0,10) +
  scale_fill_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Vegetation composition before the disturbance")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/previous_state_scen_PFT_cut.pdf", width = 10, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/previous_state_scen_PFT_cut.png", width = 10, height = 5)


# For all scenarios and PFTs together
ggplot(d_recruit, aes(x = Cluster, y=previous_state, fill = Cluster)) + 
  geom_boxplot() + theme_bw() + ylab(bquote("Aboveground carbon in kg/m"^2)) + #ylim(0,30) +
  scale_fill_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  ggtitle("Vegetation composition per cluster before the disturbance")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/previous_state.pdf", width = 10, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/previous_state.png", width = 10, height = 5)

# Cut y-axis
ggplot(d_recruit, aes(x = Cluster, y=previous_state, fill = Cluster)) + 
  geom_boxplot() + theme_bw() + ylab(bquote("Aboveground carbon in kg/m"^2)) + ylim(0,5) +
  scale_fill_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  ggtitle("Vegetation composition per cluster before the disturbance")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/previous_state_cut.pdf", width = 10, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/previous_state_cut.png", width = 10, height = 5)
