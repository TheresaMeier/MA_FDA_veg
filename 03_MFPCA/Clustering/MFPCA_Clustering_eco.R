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
library(ggridges)
library(forcats)

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

################################ Per cluster ###################################
# calculate mean values per Cluster
d_eco_mean = d_all %>%
  group_by(Cluster, Year) %>%
  summarize(mean_nuptake = mean(Nuptake, na.rm = T),
            mean_nuptake_total = mean(Nuptake_total, na.rm = T)) %>%
  ungroup %>%
  mutate(Cluster = as.factor(Cluster))

# Nuptake
ggplot(d_eco_mean) + 
  geom_line(data = d_eco_mean, linewidth = 1, 
            aes(x = Year, y = mean_nuptake, color = Cluster)) +
  scale_x_continuous(name = "Year after disturbance", expand = c(0,0), 
                     breaks = seq(2020, 2120, by = 20)) +
  scale_y_continuous(name = "Nitrogen uptake on the grid cell level", 
                     expand = c(0,0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Nitrogen uptake per cluster averaged over disturbed grid cells") +
  theme_bw() + 
  theme(
    text = element_text(size = 15), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/nuptake.pdf", width = 10, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/nuptake.png", width = 10, height = 5)

# Nuptake total
ggplot(d_eco_mean) + 
  geom_line(data = d_eco_mean, linewidth = 1, 
            aes(x = Year, y = mean_nuptake_total, color = Cluster)) +
  scale_x_continuous(name = "Year after disturbance", expand = c(0,0), 
                     breaks = seq(2020, 2120, by = 20)) +
  scale_y_continuous(name = "Nitrogen uptake", 
                     expand = c(0,0)) +
  scale_color_manual(name = "Cluster", values = c("1" = "#F8766D", "2" = "#7CAE00" , "3" = "#00BFC4", "4" = "#C77CFF", "5" = "darkgrey"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")) +
  ggtitle("Total nitrogen uptake") +
  theme_bw() + 
  theme(
    text = element_text(size = 15), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/nuptake_total.pdf", width = 6, height = 3.5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/nuptake_total.png", width = 6, height = 3.5)

###

d_recruit = d_all %>%
  distinct(Lon,Lat,PFT,initial_recruitment, recruitment_ten_years, previous_state, time_since_dist, Scenario, Cluster) %>%
  mutate(PFT = factor(PFT, levels = c("Needleleaf evergreen", "Pioneering broadleaf", "Conifers (other)", "Temperate broadleaf", "Tundra")),
         Cluster = factor(Cluster, levels = c("1", "2", "3", "4"), labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")))

############################# Initial recruitment ##############################

d_recruit <- d_recruit %>%
  mutate(Cluster = fct_rev(Cluster))

# For all scenarios and PFTs together
ggplot(d_recruit, aes(x = initial_recruitment, y = Cluster, fill = Cluster)) + 
  geom_density_ridges(aes(height = after_stat(density)), stat = "density", scale = 0.75) + 
  theme_bw() + 
  ylab("Cluster") + xlab("Number of new seedlings") + 
  #xlim(0, 100) +
  scale_fill_manual(
    name = "Cluster", 
    values = c("Cluster 1" = "#F8766D", "Cluster 2" = "#7CAE00", "Cluster 3" = "#00BFC4", "Cluster 4" = "#C77CFF"),
    guide = "none"  # Remove legend
  ) +
  theme(
    text = element_text(size = 15), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  ggtitle("Number of new seedlings per cluster right after disturbance ")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/initial_recruitment.pdf", width = 10, height = 8)
ggsave("Scripts/Plots/MFPCA/Clusters/png/initial_recruitment.png", width = 10, height = 8)

# Cut x-axis
ggplot(d_recruit, aes(x = initial_recruitment, y = Cluster, fill = Cluster)) + 
  geom_density_ridges(aes(height = after_stat(density)), stat = "density", scale = 0.75) + 
  theme_bw() + 
  ylab("") + xlab("Number of new seedlings") + 
  xlim(0, 100) +
  scale_fill_manual(
    name = "Cluster", 
    values = c("Cluster 1" = "#F8766D", "Cluster 2" = "#7CAE00", "Cluster 3" = "#00BFC4", "Cluster 4" = "#C77CFF"),
    guide = "none"  # Remove legend
  ) +
  theme(
    text = element_text(size = 15), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  ggtitle("Number of new seedlings right after disturbance ")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/initial_recruitment_cut.pdf", width = 6, height = 3.5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/initial_recruitment_cut.png", width = 6, height = 3.5)


############################# Recruitment 10 years #############################

# For all scenarios and PFTs together

ggplot(d_recruit, aes(x = recruitment_ten_years, y = Cluster, fill = Cluster)) + 
  geom_density_ridges(aes(height = after_stat(density)), stat = "density", scale = 0.75) + 
  theme_bw() + 
  ylab("Cluster") + xlab("Number of new seedlings") + 
  #xlim(0, 800) +
  scale_fill_manual(
    name = "Cluster", 
    values = c("Cluster 1" = "#F8766D", "Cluster 2" = "#7CAE00", "Cluster 3" = "#00BFC4", "Cluster 4" = "#C77CFF"),
    guide = "none"  # Remove legend
  ) +
  theme(
    text = element_text(size = 15), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  ggtitle("Number of new seedling per cluster in the ten years after disturbance (summed up)")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/recruitment_10.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/recruitment_10.png", width = 8, height = 5)

# Cut x-axis
ggplot(d_recruit, aes(x = recruitment_ten_years, y = Cluster, fill = Cluster)) + 
  geom_density_ridges(aes(height = after_stat(density)), stat = "density", scale = 0.75) + 
  theme_bw() + 
  ylab("") + xlab("Number of new seedlings") + 
  xlim(0, 800) +
  scale_fill_manual(
    name = "Cluster", 
    values = c("Cluster 1" = "#F8766D", "Cluster 2" = "#7CAE00", "Cluster 3" = "#00BFC4", "Cluster 4" = "#C77CFF"),
    guide = "none"  # Remove legend
  ) +
  theme(
    text = element_text(size = 15), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  ggtitle("Number of new seedlings after ten years")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/recruitment_10_cut.pdf", width = 6, height = 3.5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/recruitment_10_cut.png", width = 6, height = 3.5)


############################## previous state ##################################

# For all scenarios and PFTs together

ggplot(d_recruit, aes(x = previous_state, y = Cluster, fill = Cluster)) + 
  geom_density_ridges(aes(height = after_stat(density)), stat = "density", scale = 0.75) + 
  theme_bw() + 
  ylab("") + xlab(bquote("Aboveground carbon in kg/m"^2)) + 
  #xlim(0, 2.5) +
  scale_fill_manual(
    name = "Cluster", 
    values = c("Cluster 1" = "#F8766D", "Cluster 2" = "#7CAE00", "Cluster 3" = "#00BFC4", "Cluster 4" = "#C77CFF"),
    guide = "none"  # Remove legend
  ) +
  theme(
    text = element_text(size = 15), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  ggtitle("Vegetation composition per cluster before the disturbance")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/previous_state.pdf", width = 10, height = 8)
ggsave("Scripts/Plots/MFPCA/Clusters/png/previous_state.png", width = 10, height = 8)

# Cut x-axis
ggplot(d_recruit, aes(x = previous_state, y = Cluster, fill = Cluster)) + 
  geom_density_ridges(aes(height = after_stat(density)), stat = "density", scale = 0.75) + 
  theme_bw() + 
  ylab("") + xlab(bquote("Aboveground carbon in kg/m"^2)) + 
  xlim(0, 2.5) +
  scale_fill_manual(
    name = "Cluster", 
    values = c("Cluster 1" = "#F8766D", "Cluster 2" = "#7CAE00", "Cluster 3" = "#00BFC4", "Cluster 4" = "#C77CFF"),
    guide = "none"  # Remove legend
  ) +
  theme(
    text = element_text(size = 15), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  ggtitle("Vegetation composition before the disturbance")

ggsave("Scripts/Plots/MFPCA/Clusters/pdf/previous_state_cut.pdf", width = 6, height = 3.5)
ggsave("Scripts/Plots/MFPCA/Clusters/png/previous_state_cut.png", width = 6, height = 3.5)
