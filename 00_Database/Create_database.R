################################################################################
############################ Master's Thesis ###################################
################################################################################

######################### Create database (by Lucia) ###########################

setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit") #specify working directory. data should lie in this folder too, in a folder called data

library(duckdb)
library(tidyverse)
library(purrr)

con = dbConnect(duckdb(), dbdir = "patches_all.duckdb", read_only = FALSE) #create the database
dbListTables(con) #check, should be empty

#function to load one run folder of one scenario of one variable
read_write_dataset_var = function(configuration, run, variable) {
  print(paste(configuration, "run", run, " for variable ", variable))
  
  configuration = as.character(configuration) # purrr::pmap turns these into factors, but we want actual strings
  variable = as.character(variable) # purrr::pmap turns these into factors, but we want actual strings
  
  df = readr::read_table(paste0("Data/", configuration, "/run", run, "/vegstruct_patch.out"), show_col_types = F) %>%
    filter(!is.na(Year)) %>%
    select(Year, age, Lon, Lat, PID, PFT, !!rlang::sym(variable)) %>%
    pivot_wider(names_from = PFT, values_from = !!rlang::sym(variable)) %>%
    mutate(across(everything(), ~replace_na(.x, 0))) %>% # add in 0 where PFT is NA because not present
    mutate(Tundra = C3G +  HSE + HSS + LSS + LSE + GRT + EPDS + SPDS + CLM, #aggregate Tundra PFTs
           otherC = BNS + BINE + TeNE) %>%  #aggregate other conifers
    select(-c("C3G",  "HSE", "HSS", "LSS", "LSE", "GRT", "EPDS", "SPDS", "CLM", "BNS", "BINE", "TeNE")) %>%
    mutate(dhist = !sign(age)) %>% #calculate dhist
    group_by(Lon, Lat, PID) %>%
    mutate(ndist = cumsum(dhist)) %>% #calculate ndist
    pivot_longer(cols = -c(Year, Lon, Lat, PID, dhist, ndist, age), names_to = "PFT", values_to = variable) #bring back in long format
  
  dbWriteTable(con, paste0(configuration, "_", variable), df, append = T) 
  
  rm(df)
  
  gc()
  
  return("sucessful")
  
} 

s = c("ssp370", "ssp126", "ssp585", "picontrol") #specify scenarios
d = c("150") #specify disturbance of scenario
v = c("cmass", "lai", "dens", "exp_est", "fpc") #specify variables 

configurations = paste0(s, "_d", d)  #create scenario tags

input = expand.grid(scenario = configurations, run = seq(1, 160), vars = v) #create a table with all unique combinations of scenario, run folder and variable

purrr::pmap(input, ~ read_write_dataset_var(..1, ..2, ..3)) #give table and function to purrr to map over all combinations
