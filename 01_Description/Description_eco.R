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

# calculate mean values
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
  scale_y_continuous(name = "Nitrogen uptake per PFT on the grid cell level", 
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
  scale_y_continuous(name = "Nitrogen uptake per PFT on the grid cell level", 
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
  scale_y_continuous(name = "Total nitrogen uptake of the grid cell", expand = c(0, 0)) +
  scale_color_manual(name = "Scenario", drop = TRUE,
                     values = c("Control" = "darkorange", "SSP1-RCP2.6" = "green", "SSP3-RCP7.0" = "turquoise", "SSP5-RCP8.5" = "magenta3")) +
  ggtitle("Total nitrogen uptake of the grid cell") +
  theme_bw() + 
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )
ggsave("Scripts/Plots/Descriptive/Eco/pdf/nuptake_total.pdf", width = 10, height = 8)
ggsave("Scripts/Plots/Descriptive/Eco/png/nuptake_total.png", width = 10, height = 8)

d_recruit = d_eco %>%
  distinct(Lon,Lat,PFT,initial_recruitment, recruitment_ten_years, previous_state, time_since_dist, Scenario) %>%
  mutate(PFT = factor(PFT, levels = c("Needleleaf evergreen", "Pioneering broadleaf", "Conifers (other)", "Temperate broadleaf", "Tundra")))

## Initial recruitment
ggplot(d_recruit, aes(x=Scenario, y=initial_recruitment, fill = PFT)) + 
  geom_boxplot() + theme_bw() + ylab("Number of new seedlings") + #ylim(0,100) +
  scale_fill_manual(name = "Dominant vegetation", drop = TRUE,
                    values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                               "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Number of new seedlings per PFT right after disturbance ")

ggsave("Scripts/Plots/Descriptive/Eco/pdf/initial_recruitment.pdf", width = 10, height = 6)
ggsave("Scripts/Plots/Descriptive/Eco/png/initial_recruitment.png", width = 10, height = 6)

# cut plot
ggplot(d_recruit, aes(x=Scenario, y=initial_recruitment, fill = PFT)) + 
  geom_boxplot() + theme_bw() + ylab("Number of new seedlings") + ylim(0,100) +
  scale_fill_manual(name = "Dominant vegetation", drop = TRUE,
                    values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                               "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Number of new seedlings per PFT right after disturbance ")

ggsave("Scripts/Plots/Descriptive/Eco/pdf/initial_recruitment_cut.pdf", width = 10, height = 6)
ggsave("Scripts/Plots/Descriptive/Eco/png/initial_recruitment_cut.png", width = 10, height = 6)


## Recruitment 10 years
ggplot(d_recruit, aes(x=Scenario, y=recruitment_ten_years, fill = PFT)) + 
  geom_boxplot() + theme_bw() + ylab("Number of new seedlings") + #ylim(0,100) +
  scale_fill_manual(name = "Dominant vegetation", drop = TRUE,
                    values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                               "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Number of new seedling per PFT in the ten years after disturbance (summed up)")

ggsave("Scripts/Plots/Descriptive/Eco/pdf/recruitment_ten_years.pdf", width = 10, height = 6)
ggsave("Scripts/Plots/Descriptive/Eco/png/recruitment_ten_years.png", width = 10, height = 6)

# cut plot
ggplot(d_recruit, aes(x=Scenario, y=recruitment_ten_years, fill = PFT)) + 
  geom_boxplot() + theme_bw() + ylab("Number of new seedlings") + ylim(0,1000) +
  scale_fill_manual(name = "Dominant vegetation", drop = TRUE,
                    values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                               "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Number of new seedling per PFT in the ten years after disturbance (summed up)")

ggsave("Scripts/Plots/Descriptive/Eco/pdf/recruitment_ten_years_cut.pdf", width = 10, height = 6)
ggsave("Scripts/Plots/Descriptive/Eco/png/recruitment_ten_years_cut.png", width = 10, height = 6)


## previous state
ggplot(d_recruit, aes(x=Scenario, y=previous_state, fill = PFT)) + 
  geom_boxplot() + theme_bw() + ylab(bquote("Aboveground carbon in kg/m"^2)) + #ylim(0,100) +
  scale_fill_manual(name = "Dominant vegetation", drop = TRUE,
                    values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                               "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Vegetation composition before the disturbance")

ggsave("Scripts/Plots/Descriptive/Eco/pdf/previous_state.pdf", width = 10, height = 6)
ggsave("Scripts/Plots/Descriptive/Eco/png/previous_state.png", width = 10, height = 6)

# cut plot
ggplot(d_recruit, aes(x=Scenario, y=previous_state, fill = PFT)) + 
  geom_boxplot() + theme_bw() + ylab(bquote("Aboveground carbon in kg/m"^2)) + ylim(0,10) +
  scale_fill_manual(name = "Dominant vegetation", drop = TRUE,
                    values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                               "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Vegetation composition before the disturbance")

ggsave("Scripts/Plots/Descriptive/Eco/pdf/previous_state_cut.pdf", width = 10, height = 6)
ggsave("Scripts/Plots/Descriptive/Eco/png/previous_state_cut.png", width = 10, height = 6)

############################## Spatial Distribution ############################

register_google(key = "AIzaSyATIWAZ4gtgJhdH6GP_E8iSubdFh6XQ32Y")

# Get map
worldmap <- get_map(location = c(lon = -4.068561, lat = 58.87355), zoom = 1)

ggmap(worldmap) +
  geom_point(data = d_recruit, 
             aes(x = Lon, y = Lat, color = initial_recruitment), 
             size = 1) +
  facet_grid(Scenario ~ PFT) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("Number of new seedlings per PFT right after disturbance") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        strip.text.x = element_text(size = 15, face = "bold", angle = 0, hjust = 0.5),
        strip.text.y = element_text(size = 15, face = "bold", angle = 0, hjust = 0.5)) +
  labs(color = "Number of new seedlings") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 50, limits = c(0, 100))

ggsave("Scripts/Plots/Descriptive/Eco/pdf/initial_recruitment_map.pdf", width = 25, height = 10)
ggsave("Scripts/Plots/Descriptive/Eco/png/initial_recruitment_map.png", width = 25, height = 10)

ggmap(worldmap) +
  geom_point(data = d_recruit, 
             aes(x = Lon, y = Lat, color = recruitment_ten_years), 
             size = 1) +
  facet_grid(Scenario ~ PFT) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("Number of new seedling per PFT in the ten years after disturbance (summed up)") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        strip.text.x = element_text(size = 15, face = "bold", angle = 0, hjust = 0.5),
        strip.text.y = element_text(size = 15, face = "bold", angle = 0, hjust = 0.5)) +
  labs(color = "Number of new seedlings") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 500, limits = c(0, 1000))

ggsave("Scripts/Plots/Descriptive/Eco/pdf/recruitment_ten_years_map.pdf", width = 25, height = 10)
ggsave("Scripts/Plots/Descriptive/Eco/png/recruitment_ten_years_map.png", width = 25, height = 10)

ggmap(worldmap) +
  geom_point(data = d_recruit, 
             aes(x = Lon, y = Lat, color = previous_state), 
             size = 1) +
  facet_grid(Scenario ~ PFT) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("Vegetation composition before the disturbance") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        strip.text.x = element_text(size = 15, face = "bold", angle = 0, hjust = 0.5),
        strip.text.y = element_text(size = 15, face = "bold", angle = 0, hjust = 0.5)) +
  labs(color = bquote("Aboveground carbon in kg/m"^2)) +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 5, limits = c(0, 10))

ggsave("Scripts/Plots/Descriptive/Eco/pdf/previous_state_map.pdf", width = 25, height = 10)
ggsave("Scripts/Plots/Descriptive/Eco/png/previous_state_map.png", width = 25, height = 10)

