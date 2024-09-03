################################################################################
############################ Master's Thesis ###################################
################################################################################

########################### Description: Soil data #############################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/01_Description/utils.R")

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
  
  d_soil_scen = d_scen_2015_2040 %>%
    distinct(Lon,Lat,sand_fraction,silt_fraction, clay_fraction, bulkdensity_soil, ph_soil, soilcarbon) %>%
    mutate(Scenario = long_names_scenarios(scen))
  
  assign(paste0("d_soil_", scen), d_soil_scen)
}

d_soil = rbind(d_soil_picontrol, d_soil_ssp126, d_soil_ssp370, d_soil_ssp585) %>%
  distinct(Lon,Lat,sand_fraction, silt_fraction,clay_fraction,bulkdensity_soil,ph_soil,soilcarbon)

d_soil_long = melt(setDT(d_soil[,c(1:5)]), id.vars = c("Lon", "Lat"), variable.name = "property")
d_soil_long$property = gsub("_", " ", d_soil_long$property)
d_soil_long$property = gsub("sand", "Sand ", d_soil_long$property)
d_soil_long$property = gsub("clay", "Clay ", d_soil_long$property)
d_soil_long$property = gsub("silt", "Silt ", d_soil_long$property)

register_google(key = "AIzaSyA_eUpOhj7hoPyzyynWvyMqcGEA1Z_SZVY")

# Get map
worldmap <- get_map(location = c(lon =-4.068561, lat = 58.87355), zoom = 1)

########################### Soil composition ####################################
ggmap(worldmap) +
  geom_point(data = d_soil_long, 
             aes(x = Lon, y = Lat, color = value), 
             size = 1.5) +
  facet_wrap(~ property, ncol = 1) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("Soil composition for all disturbed grid cells") +
  theme(text = element_text(size = 15),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      strip.text.x = element_text(size = 15, face = "bold", angle = 0, hjust = 0.5),
      strip.text.y = element_text(size = 15, face = "bold", angle = 0, hjust = 0.5)) +
  labs(color = "Fraction") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 0.5, limits = c(0, 1))

ggsave("Scripts/Plots/Descriptive/Soil/pdf/soil_map.pdf", width = 9, height = 6)
ggsave("Scripts/Plots/Descriptive/Soil/png/soil_map.png", width = 9, height = 6)

# As ridge plot
ggplot(d_soil_long, aes(x = value)) + 
  geom_density_ridges(aes(height = after_stat(density)), stat = "density", scale = 0.75) + 
  theme_bw() +  facet_wrap(~property) +
  ylab("Density") + xlab("Fraction") + 
   theme(
    text = element_text(size = 20), 
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave("Scripts/Plots/Descriptive/Soil/pdf/soil_ridgeplot.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/Descriptive/Soil/png/soil_ridgeplot.png", width = 8, height = 5)


################################# Bulkdensity ##################################

ggmap(worldmap) +
  geom_point(data = d_soil, 
             aes(x = Lon, y = Lat, color = bulkdensity_soil), 
             size = 1.5) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("Bulk density for all disturbed grid cells") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(color = "Value") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 0.85, limits = c(0.1, 1.7))

ggsave("Scripts/Plots/Descriptive/Soil/pdf/bulk_map.pdf", width = 9, height = 6)
ggsave("Scripts/Plots/Descriptive/Soil/png/bulk_map.png", width = 9, height = 6)

################################## PH value ####################################

ggmap(worldmap) +
  geom_point(data = d_soil, 
             aes(x = Lon, y = Lat, color = ph_soil), 
             size = 1.5) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("pH value in water for all disturbed grid cells") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(color = "Value") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 6.4, limits = c(4.5, 8.2))

ggsave("Scripts/Plots/Descriptive/Soil/pdf/ph_map.pdf", width = 9, height = 6)
ggsave("Scripts/Plots/Descriptive/Soil/png/ph_map.png", width = 9, height = 6)

################################# Soil carbon ##################################

ggmap(worldmap) +
  geom_point(data = d_soil, 
             aes(x = Lon, y = Lat, color = soilcarbon), 
             size = 1.5) +
  xlab("Longitude") + ylim(45,70) + xlim(-200,200) +
  ylab("Latitude") + theme_bw() + ggtitle("Organic carbon content for all disturbed grid cells") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  labs(color = "Value") +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12)) +
  scale_color_gradient2(low = "slateblue4", mid = "floralwhite", high = "firebrick1", midpoint = 9.4, limits = c(0.5, 17.3))

ggsave("Scripts/Plots/Descriptive/Soil/pdf/soilcarbon_map.pdf", width = 9, height = 6)
ggsave("Scripts/Plots/Descriptive/Soil/png/soilcarbon_map.png", width = 9, height = 6)

# As ridge plots

d_soil_long_2 = melt(setDT(d_soil[,c(1,2,6:8)]), id.vars = c("Lon", "Lat"), variable.name = "property")
d_soil_long_2$property = ifelse(d_soil_long_2$property == "bulkdensity_soil", "Bulk density", 
                 ifelse(d_soil_long_2$property == "ph_soil", "pH in water",
                         ifelse(d_soil_long_2$property == "soilcarbon", "Organic carbon content", NA)))

ggplot(d_soil_long_2, aes(x = value)) + 
  geom_density_ridges(aes(height = after_stat(density)), stat = "density", scale = 0.75) + 
  theme_bw() +  facet_wrap(~property, scales = "free") +
  ylab("Density") + xlab("Value") + 
  theme(
    text = element_text(size = 18), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave("Scripts/Plots/Descriptive/Soil/pdf/soil_attributes_ridgeplot.pdf", width = 8, height = 5)
ggsave("Scripts/Plots/Descriptive/Soil/png/soil_attributes_ridgeplot.png", width = 8, height = 5)
