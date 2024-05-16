# Install and load necessary libraries
library(leaflet)

setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/Description/utils.R")

# Load libraries
library(duckdb)
library(tidyverse)
library(ggplot2)

# Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)

data_loc <- dbGetQuery(con, paste("SELECT Lon,Lat FROM 'picontrol_d150_cmass' WHERE PID == 0")) %>% unique()

# Create a data frame
data <- data.frame(lon = data_loc$Lon, lat = data_loc$Lat)

# Create a leaflet map
leaflet(data = data) %>%
  addTiles() %>%
  addMarkers(~lon, ~lat)
