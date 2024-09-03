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


plot.pca.fd_2 <- function(x, nx = 128, pointplot = TRUE, harm = 0,
                          expand = 0, cycle = FALSE, main.user = "", ...)
{
  #
  # Credits go to J.Ramsay who developed the package 'fda' 
  # (http://www.functionaldata.org)
  #
  # Differences to the original: main and axis labels, coloring of curves
  #
  #	Plots the harmonics produced by PCA.FD.
  #
  #  If pointplot=TRUE, then the harmonics are plotted as + and -
  #  otherwise lines are used.	Another thing that needs doing is an
  #		 arrowplot option.
  #
  # If harm = 0 (the default) then all the computed harmonics are plotted.
  #	 Otherwise those in jharm are plotted.
  # If expand =0 then effect of +/- 2 standard deviations of each pc are given
  #	 otherwise the factor expand is used.
  # If cycle=TRUE and there are 2 variables then a cycle plot will be drawn
  #	If the number of variables is anything else, cycle will be ignored.
  #
  
  # last modified 2007 May 3 by Spencer Graves	
  #	previously modified 20 March 2006
  
  pcafd <- x
  if (!(inherits(pcafd, "pca.fd"))) 
    stop("Argument 'x' is not a pca.fd object.")
  
  harmfd	<- pcafd[[1]]
  basisfd <- harmfd$basis
  rangex	<- basisfd$rangeval
  {
    if(length(nx)>1){
      x <- nx
      nx <- length(x)
    }
    else    
      x <- seq(rangex[1], rangex[2], length = nx)
  }
  fdmat	 <- eval.fd(x, harmfd)
  meanmat <- eval.fd(x, pcafd$meanfd)
  dimfd	 <- dim(fdmat)
  nharm	 <- dimfd[2]
  #
  # check number of panels
  plotsPerPg <- sum(par("mfrow"))
  #   
  harm	<- as.vector(harm)
  if(harm[1] == 0) harm <- (1:nharm)
  #  
  if(length(dimfd) == 2) {
    for(jharm in 1:length(harm)) {
      #    for(iharm in harm) {
      if(jharm==2){
        op <- par(ask=TRUE)
        on.exit(par(op)) 
      }
      #        
      iharm <- harm[jharm] 
      if(expand == 0) {
        fac <- sqrt(pcafd$values[iharm])
      } else {
        fac <- expand
      }
      #      
      vecharm <- fdmat[, iharm]
      pcmat <- cbind(meanmat + fac * vecharm, meanmat - fac * vecharm)
      if (pointplot) plottype <- "p" else plottype <- "l"
      percentvar <- round(100 * pcafd$varprop[iharm], 1)
      plot(x, meanmat, type = "l", ylim=c(min(pcmat),max(pcmat)),
           ylab=paste("Share of aboveground carbon"),               
           xlab = "Year after Disturbance",
           main=paste("PC", iharm, "-", main.user),
           ...)
      if (pointplot) {
        points(x, pcmat[,1], pch="+", col = "darkred")
        points(x, pcmat[,2], pch="-", col = "cornflowerblue")
      } else {
        lines(x, pcmat[,1], lty=2)
        lines(x, pcmat[,2], lty=3)
      }
    }
  } else {
    if(cycle && dimfd[3] == 2) {
      meanmat <- drop(meanmat)
      for(jharm in 1:length(harm)) {
        #      for(iharm in harm) {
        if(jharm==2){
          op <- par(ask=TRUE)
          on.exit(par(op)) 
        }
        iharm <- harm[jharm]
        #
        {
          if(expand == 0) fac <- 2 * sqrt(pcafd$values[iharm])
          else fac <- expand
        }
        matharm <- fdmat[, iharm,	]
        mat1 <- meanmat + fac * matharm
        mat2 <- meanmat - fac * matharm
        if (pointplot) plottype <- "p" else plottype <- "l"
        percentvar <- round(100 * pcafd$varprop[iharm],1)
        plot(meanmat[,1], meanmat[,2], type=plottype,
             xlim=c(min(c(mat1[,1],mat2[,1])),max(c(mat1[,1],mat2[,1]))),
             ylim=c(min(c(mat1[,2],mat2[,2])),max(c(mat1[,2],mat2[,2]))),
             main=paste("PC", iharm, "-", main.user),
             ...)
        if (pointplot) {
          points(mat1[, 1], mat1[, 2], pch="+", col = "darkred")
          points(mat2[, 1], mat2[, 2], pch="-", col = "cornflowerblue")
        }
        else {
          lines (mat1[, 1], mat1[, 2], lty=2)
          lines (mat2[, 1], mat2[, 2], lty=3)
        }
      }
    }
    else {
      for(jharm in 1:length(harm)) {
        #      for(iharm in harm) {
        if(jharm==2){
          op <- par(ask=TRUE)
          on.exit(par(op)) 
        }
        iharm <- harm[jharm]
        #        
        fac <- {
          if (expand == 0) sqrt(pcafd$values[iharm]) 
          else expand
        }
        
        meanmat <- drop(meanmat)
        matharm <- fdmat[, iharm, ]
        nvar <- dim(matharm)[2]
        for (jvar in 1:nvar) {
          pcmat <- cbind(meanmat[, jvar] + fac * matharm[, jvar],
                         meanmat[, jvar] - fac * matharm[, jvar])
          if (pointplot) plottype <- "p" else plottype <- "l"
          percentvar <- round(100 * pcafd$varprop[iharm], 1)
          plot(x, meanmat[,jvar], type=plottype,
               ylab=paste("Share of aboveground carbon"),
               xlab = "Year after Disturbance",
               sub = paste("Control"),
               main = dimnames(fdmat)[[3]][jvar],
               ...)
          if (pointplot) {
            points(x, pcmat[,1], pch="+", col = "darkred")
            points(x, pcmat[,2], pch="-", col = "cornflowerblue")
          }
          else {
            lines (x, pcmat[,1], lty=2)
            lines (x, pcmat[,2], lty=3)
          }
        }
      }
    }
  }
  invisible(NULL)
}

