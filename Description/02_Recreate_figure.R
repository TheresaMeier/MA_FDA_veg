setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit") # specific working dir
source("Scripts/Description/utils.R")

library(duckdb)
library(tidyverse)

con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) #create the database
dbListTables(con)

scenario = "picontrol"
variable = "cmass"
end_year = 2040
start_year = 2015

# this function will get data for one scenario. I wrapped it in a function to
# make it tidy, you can also set the variables manually (as above) and go through
# it step by step if this helps understanding.

# PID = patch id
get_data_scenario = function(scenario, start_year, end_year, variable) {
  # ok, so we want all recovery trajectories where a disturbance happened between
  # `start_year` and `and_year` and that were recovering for at last 100 years.
  
  # first we need an unique identifier of all patches disturbed between `start_year` and `end_year`
  # we filter for dhist = 1 to only get disturbed patches and for PFT = BNE, as it will make the table smaller and for now we only want the  (Lon, Lat, PID) identifiers
  locations_disturbed = dbGetQuery(con, paste0("SELECT PID, Lon, Lat, Year, ndist FROM '", scenario, "_d150_", variable, "' WHERE Year BETWEEN ", start_year, " AND ", 
                                               end_year, " AND dhist = 1 AND PFT = 'BNE';")) %>% unique()
  
  # to select only patches who where able to recover for at least 100 years, we 
  # want to inner join this table with patches in the final recovery period `start_year` + 100 - `end_year`  + 100
  # for this we need to write `locations disturbed` to the database
  dbWriteTable(con, "locations_disturbed", locations_disturbed, overwrite = T)
  
  # and perform an inner join (remember that each patch is uniquely defined by Lon, Lat and PID):
  # the resulting table `locations_disturbed_once` should have less rows than `locations_disturbed`
  locations_disturbed_once = dbGetQuery(con, paste0("SELECT d.PID, d.Lon, d.Lat, d.ndist FROM '", scenario, "_d150_", variable, 
                              "' AS d INNER JOIN locations_disturbed AS l ON d.PID = l.PID AND d.Lon = l.Lon AND d.Lat = l.Lat  WHERE d.Year BETWEEN ", 
                              start_year + 100, " AND ", end_year + 100, " AND age = 100 AND PFT = 'BNE'")) 
  
  # now we want to retrieve the whole time series, but only for these patches (again an inner join)
  # we write `locations_disturbed_once` to the database:
  dbWriteTable(con, "locations_disturbed_once", locations_disturbed_once, overwrite = T)
  
  # and join. We additionally join by `ndist`to make sure we only get that recovery trajectory.
  # for example, if a patch is disturbed in year 2030, we otherwise would additionally get the years 2015 - 2029 that we are not interested in
  df = dbGetQuery(con, paste0("SELECT d.Year, d.PFT, d.PID, d.Lon, d.Lat, d.", variable, ", d.age FROM '", scenario, "_d150_", variable, 
                              "' AS d INNER JOIN locations_disturbed_once AS l ON d.PID = l.PID AND d.Lon = l.Lon AND d.Lat = l.Lat AND d.ndist = l.ndist WHERE d.Year BETWEEN ", 
                              start_year, " AND ", end_year + 100, " AND d.PID = 1")) %>%
    group_by(age, Lon, Lat, PID) %>%
    mutate(relative = !!rlang::sym(variable)/sum(!!rlang::sym(variable)))  %>% # we calculate relative composition
    mutate(across(everything(), ~ifelse(is.na(.), 0, .)), #if sum(variable) = 0, this will be NA (can happen in the first years after a disturbance)
           scenario = scenario,
           name = paste0(long_names_scenarios(scenario))) 
 # name = paste0(long_names_scenarios(scenario), " (", start_year, " - ", end_year, ")")) 

  # we remove `locations_disturbed` again from database
  dbExecute(con, "DROP TABLE locations_disturbed")
  # we remove `locations_disturbed_once` again from database
  dbExecute(con, "DROP TABLE locations_disturbed_once")
  
  return(df)
}

trajectories_picontrol = get_data_scenario("picontrol", 2100, 2150, "cmass")
trajectories_ssp585 = get_data_scenario("ssp585", 2100, 2150, "cmass")
trajectories_ssp126 = get_data_scenario("ssp126", 2100, 2150, "cmass")
trajectories_ssp370 = get_data_scenario("ssp370", 2100, 2150, "cmass")

# stick scenarios together and make names pretty
df_trajectories = purrr::reduce(list(trajectories_picontrol, trajectories_ssp126, trajectories_ssp370, trajectories_ssp585), bind_rows) %>%
  mutate(PFT = long_names_pfts(tolower(PFT))) # make pft names pretty

# calculate mean trajectories
df_mean = df_trajectories %>%
  group_by(age, PFT, name) %>%
  summarize(relative = mean(relative, na.rm = T)) 

# plot
ggplot() + 
  geom_line(data = df_trajectories, linewidth = .05, alpha = .25,
            aes(x = age, y = relative, color = PFT, group = interaction(Lon, Lat, PID, PFT))) +
  geom_line(data = df_mean, aes(x = age, y = relative, color = PFT, group = PFT), linewidth = 2) +
  facet_grid(rows = vars(name)) +
  scale_x_continuous(name = "Year after disturbance", expand = c(0,0), limits = c(0, 100)) +
  scale_y_continuous(name = "Share of aboveground biomass", expand = c(0,0),
                     breaks = c(0.25, 0.50, 0.75, 1.00)) +
  scale_color_manual(name = "Dominant vegetation", drop = TRUE,
                    values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                               "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  ggtitle("Recovery trajectories with disturbances between 2100 and 2150") +
  add_common_layout() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10))
  

ggsave("Scripts/Plots/Descriptive/g_regeneration_2100_2150.pdf")
ggsave("Scripts/Plots/Descriptive/g_regeneration_2100_2150.png")

