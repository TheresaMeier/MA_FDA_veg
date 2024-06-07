################################################################################
############################ Master's Thesis ###################################
################################################################################

############################## Ecological data #################################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/01_Description/utils.R")
source("Scripts/MA_FDA_veg/02_FPCA/functions.R")

## Load libraries
library(dplyr)
library(data.table)
library(stringr)
library(cowplot)
library(ggplot2)
library(ggmap)

scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")

## Load data
pid = 1

for (scen in scenarios){
  d_scen_2015_2040 = fread(paste0("Data/data_", scen, "_2015_2040.csv"))
  
  d_eco_scen = d_scen_2015_2040 %>%
    filter(PID == pid) %>%
    distinct(Lon,Lat,Year,age,PFT,relative,initial_recruitment,
             recruitment_ten_years, previous_state, time_since_dist,
             Nuptake_total,Nuptake) %>%
    mutate(Scenario = long_names_scenarios(scen),
           PFT = long_names_pfts(tolower(PFT)))
  
  assign(paste0("d_eco_", scen), d_eco_scen)
}

d_eco = rbind(d_eco_picontrol, d_eco_ssp126, d_eco_ssp370, d_eco_ssp585)

# calculate mean temperatures
d_eco_mean = d_eco %>%
  group_by(Scenario, PFT, Year) %>%
  summarize(mean_nuptake = mean(Nuptake, na.rm = T),
            mean_nuptake_total = mean(Nuptake_total, na.rm = T)) 
# Nuptake
ggplot(d_eco) + 
  geom_line(data = d_eco, linewidth = .05, alpha = .25,
            aes(x = Year, y = Nuptake, color = PFT, group = interaction(Lon, Lat, PFT))) +
  # geom_line(data = d_eco_mean, aes(x = Year, y = mean_nuptake, color = PFT, group = PFT), linewidth = 2) +
  geom_smooth(data = d_eco_mean, aes(x = Year, y = mean_nuptake, color = PFT, group = PFT),
              method = "loess", se = FALSE, linewidth = 2.5, color = "black", show.legend = FALSE) +
  geom_smooth(data = d_eco_mean, aes(x = Year, y = mean_nuptake, color = PFT, group = PFT), 
              method = "loess", se = FALSE, linewidth = 1.5) +  
  facet_grid(rows = vars(Scenario)) +
  scale_x_continuous(name = "Year after disturbance", expand = c(0,0), 
                     breaks = seq(2020, 2140, by = 10)) +
  scale_y_continuous(name = "Nitrogen uptake per PFT on the gridcell level", 
                     expand = c(0,0)) +
  scale_color_manual(name = "Dominant vegetation", drop = TRUE,
                     values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                                "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  ggtitle("Nitrogen uptake per PFT for each disturbed grid cell and scenario") +
  theme_bw() + 
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )
ggsave("Scripts/Plots/Descriptive/Eco/pdf/nuptake.pdf", width = 10, height = 8)
ggsave("Scripts/Plots/Descriptive/Eco/png/nuptake.png", width = 10, height = 8)

# cut plot
ggplot(d_eco) + 
  geom_line(data = d_eco, linewidth = .05, alpha = .25,
            aes(x = Year, y = Nuptake, color = PFT, group = interaction(Lon, Lat, PFT))) +
  # geom_line(data = d_eco_mean, aes(x = Year, y = mean_nuptake, color = PFT, group = PFT), linewidth = 2) +
  geom_smooth(data = d_eco_mean, aes(x = Year, y = mean_nuptake, color = PFT, group = PFT),
              method = "loess", se = FALSE, linewidth = 2.5, color = "black", show.legend = FALSE) +
  geom_smooth(data = d_eco_mean, aes(x = Year, y = mean_nuptake, color = PFT, group = PFT), 
              method = "loess", se = FALSE, linewidth = 1.5) +  
  facet_grid(rows = vars(Scenario)) +
  scale_x_continuous(name = "Year after disturbance", expand = c(0,0), 
                     breaks = seq(2020, 2140, by = 10)) +
  scale_y_continuous(name = "Nitrogen uptake per PFT on the gridcell level", 
                     expand = c(0,0), limits = c(0,25)) +
  scale_color_manual(name = "Dominant vegetation", drop = TRUE,
                     values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                                "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  ggtitle("Nitrogen uptake per PFT for each disturbed grid cell and scenario") +
  theme_bw() + 
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )
ggsave("Scripts/Plots/Descriptive/Eco/pdf/nuptake_cut.pdf", width = 10, height = 8)
ggsave("Scripts/Plots/Descriptive/Eco/png/nuptake_cut.png", width = 10, height = 8)

# Nuptake total
ggplot(d_eco) + 
  geom_line(data = d_eco, linewidth = 0.05, alpha = .25,
            aes(x = Year, y = Nuptake_total, col = Scenario, group = interaction(Lon, Lat))) +
  geom_smooth(data = d_eco_mean, aes(x = Year, y = mean_nuptake_total, group = Scenario), 
              method = "loess", se = FALSE, linewidth = 2.5, color = "black", show.legend = FALSE) +
  geom_smooth(data = d_eco_mean, aes(x = Year, y = mean_nuptake_total, color = Scenario), 
              method = "loess", se = FALSE, linewidth = 1.5) +  
  # geom_hline(yintercept = 20, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
  # geom_hline(yintercept = 30, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
  scale_x_continuous(name = "Year", expand = c(0, 0), breaks = seq(2020, 2140, by = 10)) +
  scale_y_continuous(name = "Total nitrogen uptake of the gridcell", expand = c(0, 0)) +
  scale_color_manual(name = "Scenario", drop = TRUE,
                     values = c("Control" = "darkorange", "SSP1-RCP2.6" = "green", "SSP3-RCP7.0" = "turquoise", "SSP5-RCP8.5" = "magenta3")) +
  ggtitle("Total nitrogen uptake of the gridcell") +
  theme_bw() + 
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )
ggsave("Scripts/Plots/Descriptive/Eco/pdf/nuptake_total.pdf", width = 10, height = 8)
ggsave("Scripts/Plots/Descriptive/Eco/png/nuptake_total.png", width = 10, height = 8)

