################################################################################
############################ Master's Thesis ###################################
################################################################################

########################## Useful functions ####################################

get_data_scenario = function(scenario, start_year, end_year, variable,pid) {
  locations_disturbed = dbGetQuery(con, paste0("SELECT PID, Lon, Lat, Year, ndist FROM '", scenario, "_d150_", variable, "' WHERE Year BETWEEN ", start_year, " AND ", 
                                               end_year, " AND dhist = 1 AND PFT = 'BNE';")) %>% unique()
  
  dbWriteTable(con, "locations_disturbed", locations_disturbed, overwrite = T)
  
  locations_disturbed_once = dbGetQuery(con, paste0("SELECT d.PID, d.Lon, d.Lat, d.ndist FROM '", scenario, "_d150_", variable, 
                                                    "' AS d INNER JOIN locations_disturbed AS l ON d.PID = l.PID AND d.Lon = l.Lon AND d.Lat = l.Lat  WHERE d.Year BETWEEN ", 
                                                    start_year + 100, " AND ", end_year + 100, " AND age = 100 AND PFT = 'BNE'")) 
  
  dbWriteTable(con, "locations_disturbed_once", locations_disturbed_once, overwrite = T)
  
  df = dbGetQuery(con, paste0("SELECT d.Year, d.PFT, d.PID, d.Lon, d.Lat, d.", variable, ", d.age FROM '", scenario, "_d150_", variable, 
                              "' AS d INNER JOIN locations_disturbed_once AS l ON d.PID = l.PID AND d.Lon = l.Lon AND d.Lat = l.Lat AND d.ndist = l.ndist WHERE d.Year BETWEEN ", 
                              start_year, " AND ", end_year + 100, ifelse(pid=="all","",paste0("AND d.PID = ", pid)))) %>%
    group_by(age, Lon, Lat, PID) %>%
    mutate(relative = !!rlang::sym(variable)/sum(!!rlang::sym(variable)))  %>% # we calculate relative composition
    mutate(across(everything(), ~ifelse(is.na(.), 0, .)), #if sum(variable) = 0, this will be NA (can happen in the first years after a disturbance)
           scenario = scenario,
           name = paste0(long_names_scenarios(scenario))) %>%
    arrange(Year,PID,Lon,Lat)
  
  dbExecute(con, "DROP TABLE locations_disturbed")
  dbExecute(con, "DROP TABLE locations_disturbed_once")
  
  return(df)
}


classify_region <- function(lat, lon) {
  
  # Set bounds
  europe_bounds <- list(lat_min = 36, lat_max = 70, lon_min = -11, lon_max = 40)
  asia_bounds <- list(lat_min = -11, lat_max = 70, lon_min = 40, lon_max = 180)
  america_bounds <- list(lat_min = -56, lat_max = 70, lon_min = -180, lon_max = -11)
  
  region <- character(length(lat))
  for (i in seq_along(lat)) {
    if (lat[i] >= europe_bounds$lat_min && lat[i] <= europe_bounds$lat_max &&
        lon[i] >= europe_bounds$lon_min && lon[i] <= europe_bounds$lon_max) {
      region[i] <- "Europe"
    } else if (lat[i] >= asia_bounds$lat_min && lat[i] <= asia_bounds$lat_max &&
               lon[i] >= asia_bounds$lon_min && lon[i] <= asia_bounds$lon_max) {
      region[i] <- "Asia"
    } else if (lat[i] >= america_bounds$lat_min && lat[i] <= america_bounds$lat_max &&
               lon[i] >= america_bounds$lon_min && lon[i] <= america_bounds$lon_max) {
      region[i] <- "America"
    } else {
      region[i] <- "Other"
    }
  }
  return(region)
}


# Get basis representation
get_basis_rep = function(start_year, end_year, data, lambda = 1, norder = 6, nderiv = 3) {
  yearrange = seq(1,100,by=1)
  yearbasis = create.bspline.basis(c(1,100), ncol(data), norder)
  
  WfdPar = fdPar(yearbasis,nderiv, lambda) # 2 means it penalizes the second derivative
  
  yearrange.5 = c(1,yearrange[-1]-0.5)
  
  # fit = smooth.basis(yearrange.5, data, WfdPar)
  fit = smooth.pos(yearrange.5, data, WfdPar, dbglev = 0)
  fit$Wfdobj$fdnames = list('Year after Disturbance' = yearrange, 'Location/PID' = colnames(data), 'Share of aboveground carbon')
  
  #fit$fd$fdnames = list('Year after Disturbance' = yearrange, 'Location/PID' = colnames(data), 'Share of aboveground carbon')
  return(fit)
}

get_data_fpca = function(scenario = "picontrol", start_year, end_year, pid = 1, pft = "Tundra") {
  
  data <- get_data_scenario(scenario, start_year, end_year,'cmass',pid)
  data = data[!duplicated(data),]
  
  data_pft <- data[data$PFT == pft, c("Lon", "Lat", "age", "relative")] %>%
    filter(age <= 100 & age >0) # Filter for 100 years of recovery
  
  data_loc = data_pft[,c(1,2)] %>%
    distinct(Lon,Lat)
  
  pivot_data_pft <- data_pft %>%
    pivot_wider(names_from = c(Lon,Lat), values_from = relative) %>%
    arrange(age) %>%
    rename_with(~ gsub("_", "/", .), everything()) %>%
    ungroup()
  
  # pivot_data_pft <- pivot_data_pft %>%
  #   unnest_wider(everything(), names_sep = "/")

  return(list(as.data.frame(pivot_data_pft), data_loc))
}

get_data_fpca_cmass = function(scenario = "picontrol", start_year, end_year, pid = 1, pft = "Tundra") {
  
  data <- get_data_scenario(scenario, start_year, end_year,'cmass',pid)
  data = data[!duplicated(data),]
  
  data_pft <- data[data$PFT == pft, c("Lon", "Lat", "age", "cmass")] %>%
    filter(age <= 100 & age >0) # Filter for 100 years of recovery
  
  data_loc = data_pft[,c(1,2)] %>%
    distinct(Lon,Lat)
  
  pivot_data_pft <- data_pft %>%
    pivot_wider(names_from = c(Lon,Lat), values_from = cmass) %>%
    arrange(age) %>%
    rename_with(~ gsub("_", "/", .), everything()) %>%
    ungroup()
  
  # pivot_data_pft <- pivot_data_pft %>%
  #   unnest_wider(everything(), names_sep = "/")
  
  return(list(as.data.frame(pivot_data_pft), data_loc))
}
