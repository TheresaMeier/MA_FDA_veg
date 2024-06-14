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
  fit = smooth.pos(yearrange.5, data, WfdPar, dbglev = 1)
  fit$Wfdobj$fdnames = list('Year after Disturbance' = yearrange, 'Location/PID' = colnames(data), 'Share of aboveground carbon')
  
  #fit$fd$fdnames = list('Year after Disturbance' = yearrange, 'Location/PID' = colnames(data), 'Share of aboveground carbon')
  return(fit)
}

get_data_fpca = function(scenario = "picontrol", start_year, end_year, pid = 1, pft = "Tundra") {
  
  data <- get_data_scenario(scenario, start_year, end_year,'cmass',pid)
  
  data_pft <- data[data$PFT == pft, c("Year", "Lon", "Lat", "PID", "relative")]
  
  data_loc = data_pft %>%
    distinct(Lon,Lat,PID)
  
  pivot_data_pft <- data_pft %>%
    pivot_wider(names_from = c(Lon,Lat,PID), values_from = relative, values_fn = list) %>%
    arrange(Year) %>%
    rename_with(~ gsub("_", "/", .), everything()) 

  max_length <- max(sapply(pivot_data_pft[,-1], function(x) length(unlist(x))))
  
  # Pad shorter lists with NA values
  pivot_data_pft <- as.data.frame(lapply(pivot_data_pft, function(x) {
    if(length(unlist(x)) < max_length) {
      c(unlist(x),rep(NA, max_length - length(unlist(x))))
      #c(rep(NA, max_length - length(unlist(x))),unlist(x))
    } else {
      unlist(x)
    }
  }))
  
  pivot_data_pft = pivot_data_pft %>%
    filter(!is.na(Year)) %>%
    mutate_all(~ifelse(is.na(.), 0, .))
  
  return(list(pivot_data_pft, data_loc))
}


# Plotting 2D MFPCA results
plot.MFPCAfit2D <- function(x, plotPCs = seq_len(nObs(x$functions)), stretchFactor = NULL, combined = FALSE, ylim, ylab, cex = cex, ...)
{
  if(!inherits(x, "MFPCAfit"))
    stop("Argument is not of class 'MFPCAfit'.")
  
  if(!(is.numeric(plotPCs) & 0 < length(plotPCs) & length(plotPCs) <= length(x$values) & all(0 < plotPCs, plotPCs <= length(x$values))))
    stop("Parameter 'plotPCs' must be a vector with values between 1 and ", length(x$values), ".")
  
  if(!(is.null(stretchFactor) | all(is.numeric(stretchFactor), length(stretchFactor) == 1, stretchFactor > 0)))
    stop("Parameter 'stretchFactor' must be either NULL or a positive number.")
  
  if(!is.logical(combined))
    stop("Parameter 'combined' must be passed as a logical.")
  
  # check dimensions
  dims <- funData::dimSupp(x$functions)
  
  if(any(dims > 2))
    stop("Cannot plot principal components having a 3- or higher dimensional domain.")
  
  # Wait for user input for each new eigenfunction
  oldPar <- graphics::par(no.readonly = TRUE)
  graphics::par(ask = TRUE)
  
  # set number of rows:
  # 1: all dimensions from left to right, "+" and "-" in one plot
  # 2: all dimensions from left to right, upper row for "+", lower for "-"
  nRows <- if(combined == TRUE)
  {
    if(all(dims == 1)) {1} else {
      warning("Cannot combine plots for two-dimensional elements. Will use separate plots (combined = FALSE)")
      2
    } 
  } else{2}
  
  graphics::par(mfrow = c(nRows, length(x$functions)))
  
  
  for(ord in plotPCs) # for each order
  {
    # calculate stretch factor if not given
    if(is.null(stretchFactor))
      stretchFactor <- stats::median(abs(x$scores[,ord]))
    
    PCplus <-  stretchFactor * x$functions #x$meanFunction + stretchFactor * x$functions
    PCminus <-  -stretchFactor * x$functions #x$meanFunction - stretchFactor * x$functions
    
    min_values_X <- numeric()
    max_values_X <- numeric()    
    
    min_values_Y <- numeric()
    max_values_Y <- numeric()
    
    # Iterate over each list element in PCplus@.Data
    for (i in 1:length(PCplus@.Data)) {
      # Extract the @X component from the ith list element
      X_i <- PCplus@.Data[[i]]@X[ord,,]
      Y_i <- PCminus@.Data[[i]]@X[ord,,]
      
      # Append the minimum and maximum values of X_i to the respective vectors
      min_values_X <- c(min_values_X, min(X_i))
      max_values_X <- c(max_values_X, max(X_i))      
      
      min_values_Y <- c(min_values_Y, min(Y_i))
      max_values_Y <- c(max_values_Y, max(Y_i))
    }
    
    for(rows in seq_len(nRows))
    {
      for(i in seq_len(length(x$functions))) # for each element
      {
        #yRange <- range(PCplus[[i]]@X, PCminus[[i]]@X)
        main <- paste("PC", ord, "(explains", round(x$values[ord]/sum(x$values)*100, 2), "% of total variability)")
        
        if(dims[i] == 1)
        {
          funData::plot(x$meanFunction[[i]], lwd = 2, col = "black", 
                        main = main, ylim = ylim, ...)
          if(rows == 1)
            funData::plot(PCplus[[i]], obs = ord, 
                          add = TRUE, type = "p", pch = "+", col = "grey50", ...)
          if(rows == 2 | combined == TRUE)
            funData::plot( PCminus[[i]], obs = ord,
                           add = TRUE, type = "p", pch = "-", col = "grey50", ...)
        }  
        else # dims[i] == 2 (higher dimensional domains are not supported)
        {
          if(rows == 1)
            funData::plot(PCplus[[i]], obs = ord, main = main, ylim = ylim, zlim = c(min(min_values_X),max(max_values_X)), ylab = ylab, cex = cex, col = c(rev(colorRampPalette(brewer.pal(9, "Blues"))(ceiling(abs(min(min_values_X))))), colorRampPalette(brewer.pal(9, "Reds"))(ceiling(abs(max(max_values_X))))), ...)

          else
            funData::plot(PCminus[[i]], obs = ord, main = main, ylim = ylim, zlim = c(min(min_values_Y),max(max_values_Y)), ylab = ylab, col = c(rev(colorRampPalette(brewer.pal(9, "Blues"))(ceiling(abs(min(min_values_Y))))), colorRampPalette(brewer.pal(9, "Reds"))(ceiling(abs(max(max_values_Y))))), ...)

          
        }
      }
    }
  }
  
  graphics::par(oldPar)
  
  invisible(NULL)
}
