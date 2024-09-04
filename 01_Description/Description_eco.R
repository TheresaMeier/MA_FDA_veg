################################################################################
############################ Master's Thesis ###################################
################################################################################

######################## Description: Ecological data ##########################

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
library(ggridges)
library(forcats)

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
  
  d_eco_scen = d_eco_scen %>%
    group_by(age, Lon, Lat) %>%
    mutate(relative_previous = previous_state/sum(previous_state))  %>% # we calculate relative composition
    ungroup()
  
  assign(paste0("d_eco_", scen), d_eco_scen)
}

d_eco = rbind(d_eco_picontrol, d_eco_ssp126, d_eco_ssp370, d_eco_ssp585) %>%
  mutate(PFT = factor(PFT, levels = c("Needleleaf evergreen", "Pioneering broadleaf", "Conifers (other)", "Temperate broadleaf", "Tundra")))

# calculate mean values
d_eco_mean = d_eco %>%
  group_by(Scenario, PFT, Year) %>%
  summarize(mean_nuptake = mean(Nuptake, na.rm = T),
            mean_nuptake_total = mean(Nuptake_total, na.rm = T)) 
# Nuptake
ggplot(d_eco_mean) + 
  geom_line(data = d_eco_mean, linewidth = 1, aes(x = Year, y = mean_nuptake, color = PFT)) +
  facet_grid(rows = vars(Scenario)) +
  scale_x_continuous(name = "Year", expand = c(0,0), 
                     breaks = seq(2020, 2140, by = 10)) +
  scale_y_continuous(name = "Nitrogen uptake per PFT on the grid cell level", 
                     expand = c(0,0)) +
  scale_color_manual(name = "Dominant vegetation", drop = TRUE,
                     values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                                "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  ggtitle("Nitrogen uptake per PFT") +
  theme_bw() + theme(
    text = element_text(size = 11), 
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust=1)
  )
ggsave("Scripts/Plots/Descriptive/Eco/pdf/nuptake.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/Descriptive/Eco/png/nuptake.png", width = 8, height = 5)

# Nuptake total
ggplot(d_eco_mean) + 
  geom_line(data = d_eco_mean, linewidth = 1,
            aes(x = Year, y = mean_nuptake_total, col = Scenario)) +
  scale_x_continuous(name = "Year", expand = c(0, 0), breaks = seq(2020, 2140, by = 10)) +
  scale_y_continuous(name = "Total nitrogen uptake of the grid cell", expand = c(0, 0), limits = c(27,52)) +
  scale_color_manual(name = "Scenario", drop = TRUE,
                     values = c("Control" = "darkorange", "SSP1-RCP2.6" = "green", "SSP3-RCP7.0" = "turquoise", "SSP5-RCP8.5" = "magenta3")) +
  ggtitle("Total nitrogen uptake of the grid cell") +
  theme_bw() + theme(
    text = element_text(size = 11), 
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust=1)
  )
ggsave("Scripts/Plots/Descriptive/Eco/pdf/nuptake_total.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/Descriptive/Eco/png/nuptake_total.png", width = 8, height = 5)

d_recruit = d_eco %>%
  distinct(Lon,Lat,PFT,initial_recruitment, recruitment_ten_years, previous_state, relative_previous, time_since_dist, Scenario) %>%
  mutate(PFT = factor(PFT, levels = c("Needleleaf evergreen", "Pioneering broadleaf", "Conifers (other)", "Temperate broadleaf", "Tundra")))

## Initial recruitment
d_recruit <- d_recruit %>%
  mutate(Scenario = fct_rev(Scenario))

# cut plot
ggplot(d_recruit, aes(x = initial_recruitment, y = Scenario, fill = PFT)) + 
  geom_density_ridges(aes(height = after_stat(density)), stat = "density", scale = 0.75) + 
  theme_bw() + 
  ylab("Scenario") + xlab("Number of new seedlings") + 
  xlim(0, 120) +
  facet_grid(~PFT, scales = "free_x") + 
  scale_fill_manual(name = "Dominant vegetation", drop = TRUE,
                    values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                               "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(size = 8),
    legend.position = "none"
  ) +
  ggtitle("Number of new seedlings per PFT right after disturbance")

ggsave("Scripts/Plots/Descriptive/Eco/pdf/initial_recruitment_cut.pdf", width = 7, height = 6)
ggsave("Scripts/Plots/Descriptive/Eco/png/initial_recruitment_cut.png", width = 7, height = 6)


## Recruitment 10 years
# cut plot
ggplot(d_recruit, aes(x = recruitment_ten_years, y = Scenario, fill = PFT)) + 
  geom_density_ridges(aes(height = after_stat(density)), stat = "density", scale = 0.75) + 
  theme_bw() + 
  ylab("Scenario") + xlab("Number of new seedlings (summed up)") + 
  xlim(0, 800) +
  facet_grid(~PFT, scales = "free_x") + 
  scale_fill_manual(name = "Dominant vegetation", drop = TRUE,
                    values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                               "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(size = 8),
    legend.position = "none"
  ) +
  ggtitle("Number of new seedling per PFT ten years after disturbance")

ggsave("Scripts/Plots/Descriptive/Eco/pdf/recruitment_ten_years_cut.pdf", width = 7, height = 6)
ggsave("Scripts/Plots/Descriptive/Eco/png/recruitment_ten_years_cut.png", width = 7, height = 6)


## previous state

ggplot(d_recruit, aes(x = relative_previous, y = Scenario, fill = PFT)) + 
  geom_density_ridges(aes(height = after_stat(density)), stat = "density", scale = 0.75) + 
  theme_bw() + 
  ylab("Scenario") + xlab(bquote("Share of aboveground carbon in kg/m"^2)) +
  #xlim(0, 800) +
  facet_grid(~PFT, scales = "free_x") + 
  scale_fill_manual(name = "Dominant vegetation", drop = TRUE,
                    values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                               "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(size = 8),
    legend.position = "none"
  ) +
  ggtitle("Vegetation composition before the disturbance")

ggsave("Scripts/Plots/Descriptive/Eco/pdf/previous_state.pdf", width = 7, height = 6)
ggsave("Scripts/Plots/Descriptive/Eco/png/previous_state.png", width = 7, height = 6)

############################## Spatial Distribution ############################

register_google(key = "AIzaSyA_eUpOhj7hoPyzyynWvyMqcGEA1Z_SZVY")

# Get map
worldmap <- get_map(location = c(lon = -4.068561, lat = 58.87355), zoom = 1)

# Initial recruitment
d_initial_recruitment <- d_recruit %>%
  group_by(Lon, Lat, Scenario) %>%
  filter(initial_recruitment == max(initial_recruitment)) %>%
  ungroup()

d_initial_recruitment$Scenario <- factor(d_initial_recruitment$Scenario, levels = c("Control", "SSP1-RCP2.6", "SSP3-RCP7.0", "SSP5-RCP8.5"))

ggmap(worldmap) +
  geom_point(data = d_initial_recruitment, 
             aes(x = Lon, y = Lat, color = PFT), 
             size = 2) +
  facet_grid(rows = vars(Scenario)) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") +
  scale_color_manual(name = "Dominant vegetation", drop = TRUE,
                     values = c("Needleleaf evergreen" = "#0072B2", "Pioneering broadleaf" = "#E69F00", "Conifers (other)" = "#56B4E9",
                                "Temperate broadleaf" = "#D55E00", "Tundra" = "#009E73"),
                     labels = c("Needleleaf evergreen","Pioneering broadleaf", "Conifers (other)",
                                "Temperate broadleaf","Tundra")) +
  theme_bw() + ggtitle("New seedlings right after disturbance") +
  theme(text = element_text(size = 12),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        strip.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        strip.text.y = element_text(size = 12, angle = 0, hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  labs(color = "Dominant vegetation")

ggsave("Scripts/Plots/Descriptive/Eco/pdf/initial_recruitment_map.pdf", width = 10, height = 7.5)
ggsave("Scripts/Plots/Descriptive/Eco/png/initial_recruitment_map.png", width = 10, height = 7.5)


# Recruitment ten years
d_recruitment_ten_years <- d_recruit %>%
  group_by(Lon, Lat, Scenario) %>%
  filter(recruitment_ten_years == max(recruitment_ten_years)) %>%
  ungroup()

d_recruitment_ten_years$Scenario <- factor(d_recruitment_ten_years$Scenario, levels = c("Control", "SSP1-RCP2.6", "SSP3-RCP7.0", "SSP5-RCP8.5"))

ggmap(worldmap) +
  geom_point(data = d_recruitment_ten_years, 
             aes(x = Lon, y = Lat, color = PFT), 
             size = 2) +
  facet_grid(rows = vars(Scenario)) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") +
  scale_color_manual(name = "Dominant vegetation", drop = TRUE,
                     values = c("Needleleaf evergreen" = "#0072B2", "Pioneering broadleaf" = "#E69F00", "Conifers (other)" = "#56B4E9",
                                "Temperate broadleaf" = "#D55E00", "Tundra" = "#009E73"),
                     labels = c("Needleleaf evergreen","Pioneering broadleaf", "Conifers (other)",
                                "Temperate broadleaf","Tundra")) +
  theme_bw() + ggtitle("New seedlings right after disturbance") +
  theme(text = element_text(size = 12),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        strip.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        strip.text.y = element_text(size = 12, angle = 0, hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  labs(color = "Dominant vegetation")

ggsave("Scripts/Plots/Descriptive/Eco/pdf/recruitment_ten_years_map.pdf", width = 10, height = 7.5)
ggsave("Scripts/Plots/Descriptive/Eco/png/recruitment_ten_years_map.png", width = 10, height = 7.5)


# Previous State
d_previous_state <- d_recruit %>%
  group_by(Lon, Lat, Scenario) %>%
  filter(previous_state  == max(previous_state )) %>%
  ungroup()

d_previous_state $Scenario <- factor(d_previous_state$Scenario, levels = c("Control", "SSP1-RCP2.6", "SSP3-RCP7.0", "SSP5-RCP8.5"))

ggmap(worldmap) +
  geom_point(data = d_previous_state, 
             aes(x = Lon, y = Lat, color = PFT), 
             size = 2) +
  facet_grid(rows = vars(Scenario)) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") +
  scale_color_manual(name = "Dominant vegetation", drop = TRUE,
                     values = c("Needleleaf evergreen" = "#0072B2", "Pioneering broadleaf" = "#E69F00", "Conifers (other)" = "#56B4E9",
                                "Temperate broadleaf" = "#D55E00", "Tundra" = "#009E73"),
                     labels = c("Needleleaf evergreen","Pioneering broadleaf", "Conifers (other)",
                                "Temperate broadleaf","Tundra")) +
  theme_bw() + ggtitle("Dominant vegetation prio to the disturbance") +
  theme(text = element_text(size = 12),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        strip.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        strip.text.y = element_text(size = 12, angle = 0, hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  labs(color = "Dominant vegetation")

ggsave("Scripts/Plots/Descriptive/Eco/pdf/previous_state_map.pdf", width = 10, height = 7.5)
ggsave("Scripts/Plots/Descriptive/Eco/png/previous_state_map.png", width = 10, height = 7.5)

