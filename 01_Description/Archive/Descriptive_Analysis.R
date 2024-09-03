################################################################################
############################ Master's Thesis ###################################
################################################################################

######################## Descriptive Analysis ##################################

# Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/01_Description/utils.R")

# Load libraries
library(duckdb)
library(tidyverse)
library(ggplot2)

# Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)

############################# Annual Disturbances ###############################

# Idea: get a table of all four scenarios with summary statistics (number of total disturbances, mean/total cmass in intervals, number of gridcells etc.)

data_ndist_picontrol <- dbGetQuery(con, paste("SELECT Year, SUM(CASE WHEN dhist=1 THEN 1 ELSE 0 END) AS ndist FROM 'picontrol_d150_cmass' GROUP BY Year")) %>%
  arrange(Year) %>%
  mutate(ndist = ndist/(5*25)) %>% 
  mutate(name = 'picontrol')

data_ndist_ssp126 <- dbGetQuery(con, paste("SELECT Year, SUM(CASE WHEN dhist=1 THEN 1 ELSE 0 END) AS ndist FROM 'ssp126_d150_cmass' GROUP BY Year")) %>%
  arrange(Year) %>%
  mutate(ndist = ndist/(5*25))%>% 
  mutate(name = 'ssp126')

data_ndist_ssp370 <- dbGetQuery(con, paste("SELECT Year, SUM(CASE WHEN dhist=1 THEN 1 ELSE 0 END) AS ndist FROM 'ssp370_d150_cmass' GROUP BY Year")) %>%
  arrange(Year) %>%
  mutate(ndist = ndist/(5*25))%>% 
  mutate(name = 'ssp370')

data_ndist_ssp585 <- dbGetQuery(con, paste("SELECT Year, SUM(CASE WHEN dhist=1 THEN 1 ELSE 0 END) AS ndist FROM 'ssp585_d150_cmass' GROUP BY Year")) %>%
  arrange(Year) %>%
  mutate(ndist = ndist/(5*25))%>% 
  mutate(name = 'ssp585')

df_ndist = purrr::reduce(list(data_ndist_picontrol, data_ndist_ssp126, data_ndist_ssp370, data_ndist_ssp585), bind_rows) %>%
  mutate(name = long_names_scenarios(name))

# calculate means
df_mean = data.frame(name = c('Control','SSP1-RCP2.6', 'SSP3-RCP7.0', 'SSP5-RCP8.5'), mean_ndist = c(mean(data_ndist_picontrol$ndist), mean(data_ndist_ssp126$ndist), mean(data_ndist_ssp370$ndist), mean(data_ndist_ssp585$ndist)))

df_ndist = merge(df_ndist, df_mean, by = 'name', all.x = TRUE)

ggplot(df_ndist, aes(x = ndist)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  geom_vline(aes(xintercept = mean_ndist), data = df_mean, color = "red", linetype = "dashed", size = 1) +
  geom_text(aes(x = 36, y = 150, label = round(mean_ndist, 2)), data = df_mean, vjust = -1, color = "red") +
  facet_wrap(~ name, scales = "free") +
  labs(title = "Histogram of Annual Number of Disturbances (1800-2299)", x = "Number of disturbances", y = "Frequency") +
  theme_minimal() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10))

ggsave("Scripts/Plots/Descriptive/ndist_dist.png")
ggsave("Scripts/Plots/Descriptive/ndist_dist.pdf")


############################# Maximum Age ######################################

data_age_picontrol <- dbGetQuery(con, paste("SELECT Lon, Lat, PID, MAX(age) AS age FROM 'picontrol_d150_cmass' WHERE Year BETWEEN 2015 AND 2300 AND PFT == 'BNE' AND PID == 1 GROUP BY Lon, Lat, PID")) %>% mutate(name = 'picontrol')
data_age_ssp126 <- dbGetQuery(con, paste("SELECT Lon, Lat, PID, MAX(age) AS age FROM 'ssp126_d150_cmass' WHERE Year BETWEEN 2015 AND 2300 AND PFT == 'BNE' AND PID == 1 GROUP BY Lon, Lat, PID"))%>% mutate(name = 'ssp126')
data_age_ssp370 <- dbGetQuery(con, paste("SELECT Lon, Lat, PID, MAX(age) AS age FROM 'ssp370_d150_cmass' WHERE Year BETWEEN 2015 AND 2300 AND PFT == 'BNE' AND PID == 1 GROUP BY Lon, Lat, PID"))%>% mutate(name = 'ssp370')
data_age_ssp585 <- dbGetQuery(con, paste("SELECT Lon, Lat, PID, MAX(age) AS age FROM 'ssp585_d150_cmass' WHERE Year BETWEEN 2015 AND 2300 AND PFT == 'BNE' AND PID == 1 GROUP BY Lon, Lat, PID"))%>% mutate(name = 'ssp585')

df_age = purrr::reduce(list(data_age_picontrol, data_age_ssp126, data_age_ssp370, data_age_ssp585), bind_rows) %>%
  mutate(name = long_names_scenarios(name))

# calculate mean ages
df_mean = data.frame(name = c('Control','SSP1-RCP2.6', 'SSP3-RCP7.0', 'SSP5-RCP8.5'), mean_age = c(mean(data_age_picontrol$age), mean(data_age_ssp126$age), mean(data_age_ssp370$age), mean(data_age_ssp585$age)))

df_age = merge(df_age, df_mean, by = 'name', all.x = TRUE)

ggplot(df_age, aes(x = age)) +
      geom_histogram(binwidth = 50, fill = "skyblue", color = "black") +
      geom_vline(aes(xintercept = mean_age), data = df_mean, color = "red", linetype = "dashed", size = 1) +
      geom_text(aes(x = 480, y = 810, label = round(mean_age, 2)), data = df_mean, vjust = -1, color = "red") +
      facet_wrap(~ name, scales = "free") +
      labs(title = "Histogram of Maximum Age (2015-2300) for one Patch", x = "Age", y = "Frequency") +
      theme_minimal() +
      theme(text = element_text(size = 10),plot.title = element_text(size = 10))

ggsave("Scripts/Plots/Descriptive/age_dist.png")
ggsave("Scripts/Plots/Descriptive/age_dist.pdf")

######################## Mean above ground biomass ##############################

get_mean_var = function(scenario, variable){
  data_var <- dbGetQuery(con, paste0("SELECT PID, Year, PFT, AVG(cmass) AS cmass_avg FROM '", scenario, "_d150_", variable, 
                                       "' _d150_cmass GROUP BY PID, Year, PFT")) %>%
    arrange(PID,Year) %>%
    mutate(PFT = long_names_pfts(tolower(PFT)), # make pft names pretty
           name = paste0(long_names_scenarios(scenario)))
  
  return(data_var)
}

data_control = get_mean_var("picontrol", "cmass")
data_scen_ssp126 = get_mean_var("ssp126", "cmass")
data_scen_ssp370 = get_mean_var("ssp370", "cmass")
data_scen_ssp585 = get_mean_var("ssp585", "cmass")

data_cmass <- purrr::reduce(list(data_control, data_scen_ssp126, data_scen_ssp370, data_scen_ssp585), bind_rows)

data_cmass_avg <- data_cmass %>%
  group_by(Year, PFT, name) %>%
  summarize(cmass_avg = mean(cmass_avg, na.rm = T)) 


ggplot(data = data_cmass) + 
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
  add_common_layout() + ggtitle("Average annual aboveground biomass for all patches over the whole area")+
  theme(text = element_text(size = 10),plot.title = element_text(size = 10))

ggsave("Scripts/Plots/Descriptive/mean_cmass_dist.pdf")
ggsave("Scripts/Plots/Descriptive/mean_cmass_dist.png")

######################## Regeneration after Disturbances #######################
# Idea: Plot regeneration of plants after a pre-defined number of disturbances

# start year: 2015
get_data_regeneration = function(scenario, variable, start_year, end_year, num_dist) {

  locations_disturbed = dbGetQuery(con, paste0("SELECT PID, Lon, Lat, Year, ndist FROM '", scenario, "_d150_", variable, "' WHERE Year BETWEEN ", start_year, " AND ", end_year, "AND dhist = 1 AND ndist = ", num_dist)) %>% unique()
    
  dbWriteTable(con, "locations_disturbed", locations_disturbed, overwrite = T)
  
  locations_disturbed_num = dbGetQuery(con, paste0("SELECT d.PID, d.Lon, d.Lat, d.ndist FROM '", scenario, "_d150_", variable, 
                                                   "' AS d INNER JOIN locations_disturbed AS l ON d.PID = l.PID AND d.Lon = l.Lon AND d.Lat = l.Lat  WHERE d.Year BETWEEN ", 
                                                   start_year + 100, " AND ", end_year + 100, " AND age = 100 AND PFT = 'BNE'")) 
  
  # now we want to retrieve the whole time series, but only for these patches (again an inner join)
  # we write `locations_disturbed_once` to the database:
  dbWriteTable(con, "locations_disturbed_num", locations_disturbed_num, overwrite = T)
  
  df = dbGetQuery(con, paste0("SELECT d.Year, d.PFT, d.PID, d.Lon, d.Lat, d.", variable, ", d.age, d.ndist FROM '", scenario, "_d150_", variable, 
                              "' AS d INNER JOIN locations_disturbed_num AS l ON d.PID = l.PID AND d.Lon = l.Lon AND d.Lat = l.Lat AND d.ndist = l.ndist WHERE d.Year BETWEEN ", 
                              start_year, " AND ", end_year + 100, " AND d.PID BETWEEN 0 AND 10 AND d.ndist >= ", num_dist)) %>%
    group_by(age, Lon, Lat, PID) %>%
    mutate(relative = !!rlang::sym(variable)/sum(!!rlang::sym(variable)))  %>% # we calculate relative composition
    mutate(across(everything(), ~ifelse(is.na(.), 0, .)), #if sum(variable) = 0, this will be NA (can happen in the first years after a disturbance)
           scenario = scenario,
           name = paste0(long_names_scenarios(scenario), " (", start_year, " - ", end_year, ", ", num_dist, ")")) 
  
  # we remove `locations_disturbed` again from database
  dbExecute(con, "DROP TABLE locations_disturbed")
  
  # we remove `locations_disturbed_num` again from database
  dbExecute(con, "DROP TABLE locations_disturbed_num")
  
  return(df)
}

data_reg_ssp585_1 = get_data_regeneration('ssp585', 'cmass',2015, 2040, 1)
data_reg_ssp585_2 = get_data_regeneration('ssp585', 'cmass',2015, 2040, 2)
data_reg_ssp585_3 = get_data_regeneration('ssp585', 'cmass',2015, 2040, 3)
data_reg_ssp585_4 = get_data_regeneration('ssp585', 'cmass',2015, 2040, 4)

data_reg <- purrr::reduce(list(data_reg_ssp585_1, data_reg_ssp585_2,data_reg_ssp585_3,data_reg_ssp585_4), bind_rows)%>%
  mutate(PFT = long_names_pfts(tolower(PFT))) # make pft names pretty

data_reg_avg <- data_reg %>%
  group_by(age, PFT, name) %>%
  summarize(relative = mean(relative, na.rm = T)) 


ggplot(data = data_reg) + 
  geom_line(aes(x = age, y = relative, color = PFT, group = interaction(PID,Lon,Lat,PFT)), 
            linewidth = 0.05, alpha = 0.25) +
  geom_line(data = data_reg_avg, aes(x = age, y = relative, color = PFT, group = PFT), linewidth = 0.5) +
  facet_grid(rows = vars(name)) +
  scale_x_continuous(name = "Year after disturbance", expand = c(0,0), limits = c(0, 100)) +
  scale_y_continuous(name = "Share of aboveground biomass", expand = c(0,0),
                     breaks = c(0.25, 0.50, 0.75, 1.00)) +
  scale_color_manual(name = "Dominant vegetation", drop = TRUE,
                     values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  
                                "Needleleaf evergreen" = "#0072B2", "Conifers (other)" = "#56B4E9", 
                                "Tundra" = "#009E73")) +
  add_common_layout() +
  theme(text = element_text(size = 15)) #+ 
  # ggtitle(paste("Vegetation after ", num_disturbances, " disturbances"))

