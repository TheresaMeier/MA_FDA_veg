setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit") # specific working dir
source("Scripts/utils.R")

library(duckdb)
library(tidyverse)
library(ggplot2)

con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) #create the database
dbListTables(con)

table_name <- "ssp585_d150_cmass"

# Fetch the column names and data types for a specific table
overview <- dbGetQuery(con2, paste("SELECT column_name, data_type FROM INFORMATION_SCHEMA.COLUMNS WHERE table_name = '", table_name, "'", sep = ""))

# Print the result
print(overview)

# Fetch the number of rows in the table
row_count <- dbGetQuery(con2, paste("SELECT COUNT(*) FROM ", table_name))
print(row_count)


# Tests
# 1) Load all data for year 2099
data_2099 <- dbGetQuery(con, paste("SELECT * FROM ", table_name, " WHERE Year == 2099"))

# 2) Load one patch time series for the location (-164.25, 54.75)
data_patch <- dbGetQuery(con, paste("SELECT * FROM ", table_name, " WHERE PID == 1 AND Lon == -164.25 AND Lat == 54.75"))

# 3) Load the mean cmass per gridcell over the period 2075 - 2100
data_meanCmass <- dbGetQuery(con, paste("SELECT Lon, Lat, AVG(cmass) AS cmass_avg FROM ", table_name, " WHERE year BETWEEN 2075 AND 2100 GROUP BY Lon, Lat"))

# 4) Get the total number of disturbances per year
data_ndist <- dbGetQuery(con, paste("SELECT Year, SUM(CASE WHEN dhist=1 THEN 1 ELSE 0 END) AS ndist FROM ", table_name, " GROUP BY Year")) %>%
  arrange(Year)


# Assuming your dataset is named "data_ndist"
ggplot(data_ndist, aes(x = Year, y = ndist)) +
  geom_line() +  # Use geom_line() for a line plot
  geom_point() + # Add points for each data point
  labs(x = "Year", y = "ndist", title = "Distribution over Years") +
  theme_minimal()


ggplot(data_ndist, aes(x = Year, y = ndist)) +
  geom_line(color = "blue") +  # Use geom_line() for a line plot
  geom_point(color = "blue") + # Add points for each data point
  geom_smooth(method = "loess", color = "red", se = FALSE) + # Add a smoother
  labs(x = "Year", y = "ndist", title = "Distribution over Years") +
  theme_minimal() 



# Descriptive Analysis

get_mean_cmass = function(scenario){
  data_cmass <- dbGetQuery(con, paste("SELECT PID, Year, PFT, AVG(cmass) AS cmass_avg FROM ", scenario, " GROUP BY PID, Year, PFT")) %>%
    arrange(PID,Year)%>%
    mutate(PFT = long_names_pfts(tolower(PFT)), # make pft names pretty
            name = paste0(long_names_scenarios(scenario)))
  
  return(data_cmass)
}


data_control = get_mean_cmass("picontrol_d150_cmass")
data_scen_ssp126 = get_mean_cmass("ssp126_d150_cmass")

data_mean <- purrr::reduce(list(data_control, data_scen_ssp126), bind_rows)

data_cmass_avg <- data_mean %>%
  group_by(Year, PFT, name) %>%
  summarize(cmass_avg = mean(cmass_avg, na.rm = T)) 


ggplot(data = data_mean) + 
  geom_rect(xmin = 2015, xmax = 2100, ymin = -Inf, ymax = Inf, fill = "lightgrey", alpha = 0.9) + # Add grey box
  geom_line(aes(x = Year, y = cmass_avg, color = PFT, group = interaction(PID, PFT)), 
            linewidth = 0.1, alpha = 0.25) +
  geom_line(data = data_cmass_avg, aes(x = Year, y = cmass_avg, color = PFT, group = PFT), linewidth = 0.5) +
  facet_grid(rows = vars(name)) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Average aboveground biomass") +
  scale_color_manual(name = "Dominant vegetation", drop = TRUE,
                     values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  
                                "Needleleaf evergreen" = "#0072B2", "Conifers (other)" = "#56B4E9", 
                                "Tundra" = "#009E73")) +
  add_common_layout() +
  theme(text = element_text(size = 25))
