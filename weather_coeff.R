## Read in libraries
library(rgdal)
library(raster)
library(ncdf4)
library(sp)
#library(googledrive)

weather_coeff <- function(directory, output_directory, start, end, time_step, states_of_interest= c('California'), pest, 
                          prcp_index = 'NO', prcp_method = "threshold",  prcp_a0 = 0, prcp_a1 = 0, prcp_a2 = 0, prcp_a3 = 0, 
                          prcp_thresh = 0, prcp_x1mod = 0, prcp_x2mod = 0, prcp_x3mod = 0,
                          temp_index = 'YES', temp_method = "polynomial", temp_a0 = 0, temp_a1 = 0, temp_a2 = 0, temp_a3 = 0, 
                          temp_thresh = 0, temp_x1mod = 0, temp_x2mod = 0, temp_x3mod = 0){
  
  ## create time range
  time_range <- seq(start, end, 1)
  
  ## read in list of daymet files to choose from later 
  if(prcp_index == 'YES'){
    precip_files <- list.files(directory,pattern='prcp', full.names = TRUE)
    prcp <- stack() # Create raster stack for the area of interest and years of interest from Daymet data
    dates <- substr(precip_files,28,31) # Assumes daymet data is saved in the exact naming format that it is downloaded as
    precip_files <- precip_files[dates %in% time_range]
  } 
  
  if(temp_index == 'YES'){
    tmax_files <- list.files(directory,pattern='tmax', full.names = TRUE)
    tmin_files <- list.files(directory,pattern='tmin', full.names = TRUE)
    if(prcp_index == 'NO'){
      dates <- substr(tmax_files,28,31) # Assumes daymet data is saved in the exact naming format that it is downloaded as
    }
    tmin_files <- tmin_files[dates %in% time_range]
    tmax_files <- tmax_files[dates %in% time_range]
    ## Create raster stacks for the area of interest and years of interest from Daymet data
    tmin_s <- stack()
    tmax_s <- stack()
    tavg_s <- stack()
  }
  
  ## reference shapefile used to clip, project, and resample 
  states <- readOGR("C:/Users/Chris/Desktop/California/us_states_lccproj.shp") # link to your local copy
  reference_area <- states[states@data$STATE_NAME %in% states_of_interest,]
  rm(states)
  
  for (i in 1:length(precip_files)) {
    ## Precipitation 
    if(prcp_index == 'YES'){
      precip <- stack(precip_files[[i]], varname = "prcp")
      precip <- crop(precip, reference_area)
      precip <- mask(precip, reference_area)
      if (i>1 && compareCRS(precip,prcp) == FALSE) { precip@crs <- crs(prcp) }
      prcp <- stack(prcp, precip)
      rm(precip)
    }

    ## Temperature
    if(temp_index == 'YES'){
      tmin <- stack(tmin_files[[i]], varname = "tmin")
      tmin <- crop(tmin, reference_area)
      tmin <- mask(tmin, reference_area)
      if (i>1 && compareCRS(tmin,tmin_s) == FALSE) { tmin@crs <- crs(tmin_s) }
      tmax <- stack(tmax_files[[i]], varname = "tmax")
      tmax <- crop(tmax, reference_area)
      tmax <- mask(tmax, reference_area)
      if (i>1 && compareCRS(tmax,tmax_s) == FALSE) { tmax@crs <- crs(tmax_s) }
      tavg <- tmax
      for (j in 1:nlayers(tmax)){
        tavg[[j]] <- overlay(tmax[[j]], tmin[[j]], fun = function(r1, r2){return((r1+r2)/2)})
        print(j)
      }
      tmin_s <- stack(tmin, tmin_s)
      tmax_s <- stack(tmax, tmax_s)
      tavg_s <- stack(tavg, tavg_s)
    }
    print(i)
  }
  
  ## create indices based on timestep
  if(prcp_index =='YES'){
    if(time_step == "daily"){
      indices <- format(as.Date(names(prcp), format = "X%Y.%m.%d"), format = "%d")
      indices <- as.numeric(indices)
    } else if(time_step == "weekly"){
      indices <- rep(seq(1,((nlayers(prcp)/7)+1),1),7)
      indices <- indices[1:nlayers(prcp)]
      indices <- indices[order(indices)]
    } else if(time_step == "monthly"){
      indices <- format(as.Date(names(prcp), format = "X%Y.%m.%d"), format = "%Y%m")
      indices <- as.numeric(as.factor(indices))
    }
  }
  
  ## Create temperature and/or precipitation indices if method is threshold
  if (prcp_index == 'YES' && prcp_method == "threshold"){
    pm <- c(0, prcp_thresh, 0,  prcp_thresh, Inf, 1)
    prcp_reclass <- matrix(pm, ncol=3, byrow=TRUE)
    prcp_coeff <- reclassify(prcp,prcp_reclass)
    prcp_coeff <- stackApply(prcp_coeff, indices, fun=mean)
  }
  
  if (temp_index == 'YES' && temp_method == "threshold"){
    tm <- c(0, temp_thresh, 0,  temp_thresh, Inf, 1)
    temp_reclass <- matrix(tm, ncol=3, byrow=TRUE)
    temp_coeff <- reclassify(temp,temp_reclass)
    temp_coeff <- stackApply(temp_coeff, indices, fun=mean)
  }
  
  ## create temperature and/or precipitation indices from daymet data based on time-step and variables of interest
  if (prcp_index == 'YES' && prcp_method == "polynomial"){
    prcp_coeff <- stackApply(prcp, indices, fun=mean)
    prcp_coeff <- prcp_a0 + (prcp_a1 * (prcp_coef  + prcp_x1mod)) + (prcp_a2 * (prcp_coef + prcp_x2mod)**2) + (prcp_a3 * (prcp_coef + prcp_x3mod)**3)
    prcp_coeff[prcp_coeff < 0] <- 0 # restrain lower limit to 0
    prcp_coeff[prcp_coeff > 1] <- 1 # restrain upper limit to 1
  }
  
  if (temp_index == 'YES' && temp_method == "polynomial"){
    temp_coeff <- stackApply(tavg_s, indices, fun=mean)
    temp_coeff <- temp_a0 + (temp_a1 * (temp_coeff + temp_x1mod)) + (temp_a2 * (temp_coeff + temp_x2mod)**2) + (temp_a3 * (temp_coeff + temp_x3mod)**3)
    temp_coeff[temp_coeff < 0] <- 0 # restrain lower limit to 0
    temp_coeff[temp_coeff > 1] <- 1 # restrain upper limit to 1
  }
  
  ## create directory for writing files
  dir.create(output_directory)
  
  ## Write outputs as raster format to output directory
  writeRaster(x=prcp_coeff, filename = paste(output_directory, "/prcp_coeff_", start, "_", end, "_", pest, ".tif", sep = ""), overwrite=TRUE, format = 'GTiff')
  writeRaster(x=temp_coeff, filename = paste(output_directory, "/temp_coeff_", start, "_", end, "_", pest, ".tif", sep = ""), overwrite=TRUE, format = 'GTiff')
  
}
