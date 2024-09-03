################################################################################
############################ Master's Thesis ###################################
################################################################################

######################### Description: Climate data ############################

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
  
  d_climate_scen = d_scen_2015_2040 %>%
    filter(PID == pid) %>%
    distinct(Lon,Lat,Year,tas_yearlymin,tas_yearlymeam, tas_yearlymax,pr_yearlysum) %>%
    mutate(tas_yearlymin = tas_yearlymin - 273.15,
           tas_yearlymeam = tas_yearlymeam - 273.15,
           tas_yearlymax = tas_yearlymax - 273.15,
           Scenario = long_names_scenarios(scen),
           region = classify_region(Lat,Lon))
  
  assign(paste0("d_climate_", scen), d_climate_scen)
}

d_climate = rbind(d_climate_picontrol, d_climate_ssp126, d_climate_ssp370, d_climate_ssp585)

# calculate mean values
d_climate_mean = d_climate %>%
  group_by(Scenario, Year) %>%
  summarize(mean_tas_min = mean(tas_yearlymin, na.rm = T),
            mean_tas_mean = mean(tas_yearlymeam, na.rm = T),
            mean_tas_max = mean(tas_yearlymax, na.rm = T),
            mean_precip = mean(pr_yearlysum, na.rm = T)) 

# Yearly mean temperature
ggplot(d_climate_mean) + 
  geom_line(data = d_climate_mean, linewidth = 1,
            aes(x = Year, y = mean_tas_mean, col = Scenario)) +
  geom_hline(yintercept = 0, color = "black", size = 0.25, alpha = .75, lty = "dashed") +
  scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by = 10)) +
  scale_y_continuous(name = "Yearly mean temperature in °C", expand = c(0,0), limits = c(-3,6)) +
  scale_color_manual(name = "Scenario", drop = TRUE,
                     values = c("Control" = "darkorange", "SSP1-RCP2.6" = "green", "SSP3-RCP7.0" = "turquoise", "SSP5-RCP8.5" = "magenta3")) +
  ggtitle("Annual mean temperature") +
  theme_bw() + theme(
    text = element_text(size = 20), 
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust=1)
  )

ggsave("Scripts/Plots/Descriptive/Climate/pdf/temp_mean.pdf", width = 8, height = 6)
ggsave("Scripts/Plots/Descriptive/Climate/png/temp_mean.png", width = 8, height = 6)

# Daily minimum temperature per year
ggplot(d_climate_mean) + 
  geom_line(data = d_climate_mean, linewidth = 1,
            aes(x = Year, y = mean_tas_min, col = Scenario)) +
  scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by=10)) +
  scale_y_continuous(name = "Yearly minimum temperature in °C", expand = c(0,0), limits = c(-45,-20)) +
  scale_color_manual(name = "Scenario", drop = TRUE,
                     values = c("Control" = "darkorange", "SSP1-RCP2.6" = "green", "SSP3-RCP7.0" = "turquoise", "SSP5-RCP8.5" = "magenta3")) +
  ggtitle("Annual minimum temperature") +
  theme_bw() + theme(
    text = element_text(size = 20), 
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust=1)
  )

ggsave("Scripts/Plots/Descriptive/Climate/pdf/temp_min.pdf", width = 8, height = 6)
ggsave("Scripts/Plots/Descriptive/Climate/png/temp_min.png", width = 8, height = 6)

# Daily maximum temperature per year
ggplot(d_climate_mean) + 
  geom_line(data = d_climate_mean, linewidth = 1,
            aes(x = Year, y = mean_tas_max, col = Scenario)) +
  scale_x_continuous(name = "Year", expand = c(0,0), breaks = seq(2020,2140, by=10)) +
  scale_y_continuous(name = "Yearly maximum temperature in °C", expand = c(0,0), limits = c(20,28)) +
  scale_color_manual(name = "Scenario", drop = TRUE,
                     values = c("Control" = "darkorange", "SSP1-RCP2.6" = "green", "SSP3-RCP7.0" = "turquoise", "SSP5-RCP8.5" = "magenta3")) +
  ggtitle("Annual maximum temperature") +
  theme_bw() + theme(
    text = element_text(size = 20), 
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust=1)
  )

ggsave("Scripts/Plots/Descriptive/Climate/pdf/temp_max.pdf", width = 8, height = 6)
ggsave("Scripts/Plots/Descriptive/Climate/png/temp_max.png", width = 8, height = 6)

# Precipitation
ggplot(d_climate_mean) + 
  geom_line(data = d_climate_mean, linewidth = 1,
            aes(x = Year, y = mean_precip, col = Scenario)) +
  scale_x_continuous(name = "Year", expand = c(0, 0), breaks = seq(2020, 2140, by = 10)) +
  scale_y_continuous(name = bquote("Yearly precipitation in kg/m"^2), expand = c(0, 0), limits = c(550,910)) +
  scale_color_manual(name = "Scenario", drop = TRUE,
                     values = c("Control" = "darkorange", "SSP1-RCP2.6" = "green", "SSP3-RCP7.0" = "turquoise", "SSP5-RCP8.5" = "magenta3")) +
  ggtitle("Annual precipitation") +
  theme_bw() + theme(
    text = element_text(size = 20), 
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30,hjust=1)
  )
ggsave("Scripts/Plots/Descriptive/Climate/pdf/precip.pdf", width = 8, height = 6)
ggsave("Scripts/Plots/Descriptive/Climate/png/precip.png", width = 8, height = 6)


