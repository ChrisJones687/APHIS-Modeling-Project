## install packages (only needed on first run of the model)
#install.packages(c("rgdal","raster","lubridate","CircStats","Rcpp", "rgrass7", "optparse", "plotrix", "ncdf4", "dismo", "sp"))

## load packages:
suppressPackageStartupMessages(library(raster))    # Raster operation and I/O. Depends R (≥ 2.15.0)
suppressPackageStartupMessages(library(rgdal))     # Geospatial data abstraction library. Depends R (≥ 2.14.0)
suppressPackageStartupMessages(library(lubridate)) # Make dealing with dates a little easier. Depends R (≥ 3.0.0)
suppressPackageStartupMessages(library(CircStats)) # Circular Statistics - Von Mises distribution
suppressPackageStartupMessages(library(Rcpp))      # Seamless R and C++ Integration. Depends R (≥ 3.0.0)
suppressPackageStartupMessages(library(ncdf4))     # work with NetCDF datasets
suppressPackageStartupMessages(library(dismo))     # Regression for ecological datasets
suppressPackageStartupMessages(library(sp))        # Classes and methods for spatial data

pest <- function(host1_rast, host1_score = NULL, host2_rast=NULL, host2_score=NULL, host3_rast=NULL, host3_score=NULL, host4_rast=NULL, host4_score=NULL,
                 host5_rast=NULL, host5_score=NULL, host6_rast=NULL, host6_score=NULL, host7_rast=NULL, host7_score=NULL, host8_rast=NULL, host8_score=NULL,
                 host9_rast=NULL, host9_score=NULL, host10_rast=NULL, host10_score=NULL, allTrees, initialPopulation, start, end, seasonality = 'NO',
                 s1 = 1 , s2 = 12, sporeRate, windQ, windDir, tempQ, tempData, precipQ, precipData, kernelType ='Cauchy', kappa = 2, number_of_hosts = 1, 
                 scale1 = 20.57, scale2 = NULL, gamma = 1, seed_n = 42, time_step = "weeks"){
  
## Define the main working directory based on the current script path (un commment next line if used outside of shiny framework)
## setwd("C:\\Users\\chris\\Dropbox\\Projects\\Code\\APHIS-Modeling-Project2")

## Use an external source file w/ all modules (functions) used within this script. 
## Use FULL PATH if source file is not in the same folder w/ this script
# source('scripts/myfunctions_SOD.r') # loads custom functions for dispersal using R
sourceCpp("scripts/myCppFunctions2.cpp") # load custom functions dispersal that use C++ (Faster)
source("scripts/myfunctions_SOD.r")

host_score <- c(host1_score, host2_score, host3_score, host4_score, host5_score, host6_score, host7_score, host8_score, host9_score, host10_score)
host_score[(number_of_hosts+1):10] <-0
host_score <- host_score/10

## All live trees (for calculating the proportion of infected) (tree density per hectare)
all_trees_rast <- allTrees
all_trees_rast[is.na(all_trees_rast)]<- 0

## raster resolution
res_win <- res(host1_rast)[1]
res_win2 <- res(host1_rast)[2]
res_area <- res_win*res_win2
n_cols <- as.numeric(ncol(host1_rast))
n_rows <- as.numeric(nrow(host1_rast))

### INFECTED AND SUSCEPTIBLES ####
## Initial infection:
initialPopulation[is.na(initialPopulation)]<- 0
initial_infection <- as.matrix(initialPopulation)
number_of_hosts = number_of_hosts

## define matrices for infected species of interest
## This is a set of nested if Statements (not properly tabbed)
if (number_of_hosts>0){
  I_host1 <- matrix(0, nrow=n_rows, ncol=n_cols)
  host1_rast[is.na(host1_rast)]<- 0
  S_host1 <- as.matrix(host1_rast)
  if(any(S_host1[initial_infection > 0] > 0)) I_host1[initial_infection > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y*2)), x), S_host1[initial_infection > 0], initial_infection[initial_infection > 0])
  S_host1 <- S_host1 - I_host1 
  I_host1_rast <- initialPopulation
  I_host1_rast[] <- I_host1
  I_host1_stack <- stack(I_host1_rast)
  stack_list <- list(I_host1_stack)
  S_matrix_list <- list(S_host1)
  I_matrix_list <- list(I_host1)

if (number_of_hosts>1) {
  I_host2 <- matrix(0, nrow=n_rows, ncol=n_cols)
  host2_rast[is.na(host2_rast)]<- 0
  S_host2 <- as.matrix(host2_rast)
  if(any(S_host2[initial_infection > 0] > 0)) I_host2[initial_infection > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y*2)), x), S_host2[initial_infection > 0], initial_infection[initial_infection > 0]) 
  S_host2 <- S_host2 - I_host2 
  I_host2_rast <- initialPopulation
  I_host2_rast[] <- I_host2
  I_host2_stack <- stack(I_host2_rast)
  stack_list <- c(stack_list,I_host2_stack)
  S_matrix_list[[2]] <- S_host2
  I_matrix_list[[2]] <- I_host2

if (number_of_hosts>2) {
  I_host3 <- matrix(0, nrow=n_rows, ncol=n_cols)
  host3_rast[is.na(host3_rast)]<- 0
  S_host3 <- as.matrix(host3_rast)
  if(any(S_host3[initial_infection > 0] > 0)) I_host3[initial_infection > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y*2)), x), S_host3[initial_infection > 0], initial_infection[initial_infection > 0])
  S_host3 <- S_host3 - I_host3 
  I_host3_rast <- initialPopulation
  I_host3_rast[] <- I_host3
  I_host3_stack <- stack(I_host3_rast)
  stack_list <- c(stack_list,I_host3_stack)
  S_matrix_list[[3]] <- S_host3
  I_matrix_list[[3]] <- I_host3
  
if (number_of_hosts>3) {
  I_host4 <- matrix(0, nrow=n_rows, ncol=n_cols)
  host4_rast[is.na(host4_rast)]<- 0
  S_host4 <- as.matrix(host4_rast)
  if(any(S_host4[initial_infection > 0] > 0)) I_host4[initial_infection > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y*2)), x), S_host4[initial_infection > 0], initial_infection[initial_infection > 0])
  S_host4 <- S_host4 - I_host4 
  I_host4_rast <- initialPopulation
  I_host4_rast[] <- I_host4
  I_host4_stack <- stack(I_host4_rast)
  stack_list <- c(stack_list,I_host4_stack)
  S_matrix_list[[4]] <- S_host4
  I_matrix_list[[4]] <- I_host4

if (number_of_hosts>4) {
  I_host5 <- matrix(0, nrow=n_rows, ncol=n_cols)
  host5_rast[is.na(host5_rast)]<- 0
  S_host5 <- as.matrix(host5_rast)
  if(any(S_host5[initial_infection > 0] > 0)) I_host5[initial_infection > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y*2)), x), S_host5[initial_infection > 0], initial_infection[initial_infection > 0])
  S_host5 <- S_host5 - I_host5 
  I_host5_rast <- initialPopulation
  I_host5_rast[] <- I_host5
  I_host5_stack <- stack(I_host5_rast)
  stack_list <- c(stack_list,I_host5_stack)
  S_matrix_list[[5]] <- S_host5
  I_matrix_list[[5]] <- I_host5

if (number_of_hosts>5) {
  I_host6 <- matrix(0, nrow=n_rows, ncol=n_cols)
  host6_rast[is.na(host6_rast)]<- 0
  S_host6 <- as.matrix(host6_rast)
  if(any(S_host6[initial_infection > 0] > 0)) I_host6[initial_infection > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y*2)), x), S_host6[initial_infection > 0], initial_infection[initial_infection > 0])
  S_host6 <- S_host6 - I_host6 
  I_host6_rast <- initialPopulation
  I_host6_rast[] <- I_host6
  I_host6_stack <- stack(I_host6_rast)
  stack_list <- c(stack_list,I_host6_stack)
  S_matrix_list[[6]] <- S_host6
  I_matrix_list[[6]] <- I_host6
  
if (number_of_hosts>6) {
  I_host7 <- matrix(0, nrow=n_rows, ncol=n_cols)
  host7_rast[is.na(host7_rast)]<- 0
  S_host7 <- as.matrix(host7_rast)
  if(any(S_host7[initial_infection > 0] > 0)) I_host7[initial_infection > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y*2)), x), S_host7[initial_infection > 0], initial_infection[initial_infection > 0])
  S_host7 <- S_host7 - I_host7 
  I_host7_rast <- initialPopulation
  I_host7_rast[] <- I_host7
  I_host7_stack <- stack(I_host7_rast)
  stack_list <- c(stack_list,I_host7_stack)
  S_matrix_list[[7]] <- S_host7
  I_matrix_list[[7]] <- I_host7
  
if (number_of_hosts>7) {
  I_host8 <- matrix(0, nrow=n_rows, ncol=n_cols)
  host8_rast[is.na(host8_rast)]<- 0
  S_host8 <- as.matrix(host8_rast)
  if(any(S_host8[initial_infection > 0] > 0)) I_host8[initial_infection > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y*2)), x), S_host8[initial_infection > 0], initial_infection[initial_infection > 0])
  S_host8 <- S_host8 - I_host8 
  I_host8_rast <- initialPopulation
  I_host8_rast[] <- I_host8
  I_host8_stack <- stack(I_host8_rast)
  stack_list <- c(stack_list,I_host8_stack)
  S_matrix_list[[8]] <- S_host8
  I_matrix_list[[8]] <- I_host8
  
if (number_of_hosts>8) {
  I_host9 <- matrix(0, nrow=n_rows, ncol=n_cols)
  host9_rast[is.na(host9_rast)]<- 0
  S_host9 <- as.matrix(host9_rast)
  if(any(S_host9[initial_infection > 0] > 0)) I_host9[initial_infection > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y*2)), x), S_host9[initial_infection > 0], initial_infection[initial_infection > 0])
  S_host9 <- S_host9 - I_host9 
  I_host9_rast <- initialPopulation
  I_host9_rast[] <- I_host9
  I_host9_stack <- stack(I_host9_rast)
  stack_list <- c(stack_list,I_host9_stack)
  S_matrix_list[[9]] <- S_host9
  I_matrix_list[[9]] <- I_host9
  
if (number_of_hosts>9) {
  I_host10 <- matrix(0, nrow=n_rows, ncol=n_cols)
  host10_rast[is.na(host10_rast)]<- 0
  S_host10 <- as.matrix(host10_rast)
  if(any(S_host10[initial_infection > 0] > 0)) I_host10[initial_infection > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y*2)), x), S_host10[initial_infection > 0], initial_infection[initial_infection > 0])
  S_host10 <- S_host10 - I_host10
  I_host10_rast <- initialPopulation
  I_host10_rast[] <- I_host10
  I_host10_stack <- stack(I_host10_rast)
  stack_list <- c(stack_list,I_host10_stack)
  S_matrix_list[[10]] <- S_host10
  I_matrix_list[[10]] <- I_host10
}}}}}}}}}} 

for (i in (number_of_hosts+1):10){
  S_matrix_list[[i]] <- matrix(0, nrow=n_rows, ncol=n_cols)
  I_matrix_list[[i]] <- matrix(0, nrow=n_rows, ncol=n_cols)
}


## define matrix for all live trees (for calculating the percentage of infected)
all_trees <- as.matrix(all_trees_rast)

## Start-End date:
start = start
end = end
if (start > end) stop('start date must precede end date!!')

## build time series for simulation steps:
dd_start <- as.POSIXlt(as.Date(paste(start,'-01-01',sep='')))
dd_end <- as.POSIXlt(as.Date(paste(end,'-12-31',sep='')))
tstep <- as.character(seq(dd_start, dd_end, time_step))

## create list for yearly output
if (time_step == "weeks") {
  split_date2 = unlist(strsplit(tstep, '-'))
  split_date2 = as.data.frame(as.numeric(split_date2[seq(2,length(split_date2),3)]))
  listvar = 1
  yearlyoutputlist = 0
  for (i in 2:nrow(split_date2)) {
    if (split_date2[i,1] > split_date2[i-1,1] && split_date2[i,1] == 10) {
      yearlyoutputlist[listvar] <- i-1
      listvar = listvar +1
    } 
  }
} else if (time_step == "months") {
  n_years <- end-start+1
  yearlyoutputlist <- seq(s2, s2+(n_years-1)*12,12)
}

## Create data frame for infected host data 
years = seq(start, end, 1)
dataForOutput <- data.frame(years = years, infectedHost1Individuals = 0, infectedHost1Area = 0, infectedHost2Individuals = 0, infectedHost2Area = 0) # replace infected host with actual host names
yearTracker = 0

### WEATHER SUITABILITY: read and stack weather suitability raster BEFORE running the simulation ### 
## weather coefficients
if (tempQ == "YES" && precipQ == "YES") {
  if (extension(precipData)==".nc"){
    mcf.array <- ncvar_get(nc_open(precipData),  varid = "Mcoef") #M = moisture;
    ccf.array <- ncvar_get(nc_open(tempData),  varid = "Ccoef") #C = temperature;
  } else {
    temp_data <- stack(tempData)
    temp_data[is.na(temp_data)] <- 0
    precip_data <- stack(precipData)
    precip_data[is.na(precip_data)] <- 0
    mcf.array <- as.array(stack(precip_data))
    ccf.array <- as.array(stack(temp_data))
  }
} else if (tempQ == "YES" && precipQ == "NO") {
  if (extension(tempData)==".nc"){
    ccf.array <- ncvar_get(nc_open(tempData),  varid = "Ccoef") #C = temperature;
  } else {
    temp_data <- stack(tempData)
    temp_data[is.na(temp_data)] <- 0
    ccf.array <- as.array(temp_data)
  }
} else if (tempQ == "NO" && precipQ == "YES") {
  if (extension(precipData)==".nc"){
    mcf.array <- ncvar_get(nc_open(precipData),  varid = "Mcoef") #M = moisture;
  } else {
    precip_data <- stack(precipData)
    precip_data[is.na(precip_data)] <- 0
    mcf.array <- as.array(stack(precipData))
  }}

## Seasonality: Do you want the spread to be limited to certain months?
seasonality <- seasonality   #'YES' or 'NO'
if (seasonality == 'YES') months_msk <- formatC(s1:s2, width = 2, format = "d", flag = "0") # 1=January 12=December(Default to 1-12)

## Wind: Do you want the spread to be affected by wind?
wind <- windQ #'YES' or 'NO'
if (wind == "YES"){
  pwdir <- windDir
}

spore_rate <- sporeRate

#time counter to access pos index in weather raster stacks
cnt <- 0 
crit_cnt <- 0

## ----> MAIN SIMULATION LOOP (weekly time steps) <------
for (tt in tstep){
  
    ## check if there are any susceptible host left on the landscape (IF NOT continue LOOP till the end)
    if(!any(S_host1 > 0)) break
    
    ## update counter
    cnt <- cnt + 1
    
    ## add in removal of infected hosts when some critical limit is met (note: to generalize will need to replace of -12.87 and "-01-01". Maybe other steps as well)
    if (substr(tt,5,10) == "-01-01"){
      crit_cnt <- crit_cnt + 1
      if(any(I_matrix_list[[1]][crit_temp[,,crit_cnt] < -12.87] > 0)) {
        S_matrix_list[[1]][crit_temp[,,crit_cnt] < -12.87] <- S_matrix_list[[1]][crit_temp[,,crit_cnt] < -12.87] + I_matrix_list[[1]][crit_temp[,,crit_cnt] < -12.87]
        I_matrix_list[[1]][crit_temp[,,crit_cnt] < -12.87] <- 0
      } 
    }
    
    ## is current time step within a spread month (as defined by input parameters)?
    if (seasonality == 'YES' & !any(substr(tt,6,7) %in% months_msk)) next
    
    ## Total weather suitability:
    if (tempQ == "YES" && precipQ == "YES") {
      weather_suitability <- mcf.array[,,cnt] * ccf.array[,,cnt]
    } else if (tempQ == "YES" && precipQ == "NO") {
      weather_suitability <- ccf.array[,,cnt]
    } else if (tempQ == "NO" && precipQ == "YES") {
      weather_suitability <- mcf.array[,,cnt]
    } else if (tempQ =="NO" && precipQ=="NO"){
      weather_suitability <-  matrix(1, nrow=n_rows, ncol=n_cols)
    }
    
    ## GENERATE SPORES:  
    set.seed(seed_n)
    infected_matrix <- matrix(0, nrow=n_rows, ncol=n_cols)
    for (i in 1:number_of_hosts){
      infected_matrix <- infected_matrix + (I_matrix_list[[i]]*(host_score[i]))
    }
    spores_mat <- SporeGenCpp(infected_matrix, weather_suitability, rate = spore_rate) # rate spores/week
    
    ##SPORE DISPERSAL:  
    #'List'
    if (wind == 'YES') {
      
      #Check if predominant wind direction has been specified correctly:
      if (!(pwdir %in% c('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'))) stop('A predominant wind direction must be specified: N, NE, E, SE, S, SW, W, NW')
      if (kernelType == "Cauchy Mixture") {
        out <- SporeDispCppWind_mh(spores_mat, 
                                   S_host1_mat=S_matrix_list[[1]],S_host2_mat=S_matrix_list[[2]],S_host3_mat=S_matrix_list[[3]],S_host4_mat=S_matrix_list[[4]],S_host5_mat=S_matrix_list[[5]],
                                   S_host6_mat=S_matrix_list[[6]],S_host7_mat=S_matrix_list[[7]],S_host8_mat=S_matrix_list[[8]],S_host9_mat=S_matrix_list[[9]],S_host10_mat=S_matrix_list[[10]],
                                   I_host1_mat=I_matrix_list[[1]],I_host2_mat=I_matrix_list[[2]],I_host3_mat=I_matrix_list[[3]],I_host4_mat=I_matrix_list[[4]],I_host5_mat=I_matrix_list[[5]],
                                   I_host6_mat=I_matrix_list[[6]],I_host7_mat=I_matrix_list[[7]],I_host8_mat=I_matrix_list[[8]],I_host9_mat=I_matrix_list[[9]],I_host10_mat=I_matrix_list[[10]],
                                   N_LVE=all_trees, weather_suitability, rs=res_win, rtype=kernelType, scale1=20.57, wdir=pwdir, kappa=kappa, host_score = host_score, scale2 = scale2, gamma = gamma)
      }else{
        out <- SporeDispCppWind_mh(spores_mat, 
                                   S_host1_mat=S_matrix_list[[1]],S_host2_mat=S_matrix_list[[2]],S_host3_mat=S_matrix_list[[3]],S_host4_mat=S_matrix_list[[4]],S_host5_mat=S_matrix_list[[5]],
                                   S_host6_mat=S_matrix_list[[6]],S_host7_mat=S_matrix_list[[7]],S_host8_mat=S_matrix_list[[8]],S_host9_mat=S_matrix_list[[9]],S_host10_mat=S_matrix_list[[10]],
                                   I_host1_mat=I_matrix_list[[1]],I_host2_mat=I_matrix_list[[2]],I_host3_mat=I_matrix_list[[3]],I_host4_mat=I_matrix_list[[4]],I_host5_mat=I_matrix_list[[5]],
                                   I_host6_mat=I_matrix_list[[6]],I_host7_mat=I_matrix_list[[7]],I_host8_mat=I_matrix_list[[8]],I_host9_mat=I_matrix_list[[9]],I_host10_mat=I_matrix_list[[10]],
                                   N_LVE=all_trees, weather_suitability, rs=res_win, rtype=kernelType, scale1=20.57, wdir=pwdir, kappa=kappa, host_score = host_score)
      }
      
    }else{
      if (kernelType == "Cauchy Mixture") {
        out <- SporeDispCpp_mh(spores_mat, 
                               S_host1_mat=S_matrix_list[[1]],S_host2_mat=S_matrix_list[[2]],S_host3_mat=S_matrix_list[[3]],S_host4_mat=S_matrix_list[[4]],S_host5_mat=S_matrix_list[[5]],
                               S_host6_mat=S_matrix_list[[6]],S_host7_mat=S_matrix_list[[7]],S_host8_mat=S_matrix_list[[8]],S_host9_mat=S_matrix_list[[9]],S_host10_mat=S_matrix_list[[10]],
                               I_host1_mat=I_matrix_list[[1]],I_host2_mat=I_matrix_list[[2]],I_host3_mat=I_matrix_list[[3]],I_host4_mat=I_matrix_list[[4]],I_host5_mat=I_matrix_list[[5]],
                               I_host6_mat=I_matrix_list[[6]],I_host7_mat=I_matrix_list[[7]],I_host8_mat=I_matrix_list[[8]],I_host9_mat=I_matrix_list[[9]],I_host10_mat=I_matrix_list[[10]],
                               N_LVE=all_trees, weather_suitability, rs=res_win, rtype=kernelType, scale1=20.57, host_score = host_score, gamma = gamma, scale2 = scale2) ##TO DO
      }else{
        out <- SporeDispCpp_mh(spores_mat, 
                               S_host1_mat=S_matrix_list[[1]],S_host2_mat=S_matrix_list[[2]],S_host3_mat=S_matrix_list[[3]],S_host4_mat=S_matrix_list[[4]],S_host5_mat=S_matrix_list[[5]],
                               S_host6_mat=S_matrix_list[[6]],S_host7_mat=S_matrix_list[[7]],S_host8_mat=S_matrix_list[[8]],S_host9_mat=S_matrix_list[[9]],S_host10_mat=S_matrix_list[[10]],
                               I_host1_mat=I_matrix_list[[1]],I_host2_mat=I_matrix_list[[2]],I_host3_mat=I_matrix_list[[3]],I_host4_mat=I_matrix_list[[4]],I_host5_mat=I_matrix_list[[5]],
                               I_host6_mat=I_matrix_list[[6]],I_host7_mat=I_matrix_list[[7]],I_host8_mat=I_matrix_list[[8]],I_host9_mat=I_matrix_list[[9]],I_host10_mat=I_matrix_list[[10]],
                               N_LVE=all_trees, weather_suitability, rs=res_win, rtype=kernelType, scale1=20.57, host_score = host_score) ##TO DO
      }
    }  
    
    ## update R matrices: ## Note this is a set of nested if statements
    if (number_of_hosts>0){
    S_matrix_list[[1]] <- out$S_host1_mat
    I_matrix_list[[1]] <- out$I_host1_mat
    if (number_of_hosts>1){
    S_matrix_list[[2]] <- out$S_host2_mat
    I_matrix_list[[2]] <- out$I_host2_mat
    if (number_of_hosts>2){
    S_matrix_list[[3]] <- out$S_host3_mat
    I_matrix_list[[3]] <- out$I_host3_mat
    if (number_of_hosts>3){
    S_matrix_list[[4]] <- out$S_host4_mat
    I_matrix_list[[4]] <- out$I_host4_mat
    if (number_of_hosts>4){
    S_matrix_list[[5]] <- out$S_host5_mat
    I_matrix_list[[5]] <- out$I_host5_mat
    if (number_of_hosts>5){
    S_matrix_list[[6]] <- out$S_host6_mat
    I_matrix_list[[6]] <- out$I_host6_mat
    if (number_of_hosts>6){
    S_matrix_list[[7]] <- out$S_host7_mat
    I_matrix_list[[7]] <- out$I_host7_mat
    if (number_of_hosts>7){
    S_matrix_list[[8]] <- out$S_host8_mat
    I_matrix_list[[8]] <- out$I_host8_mat
    if (number_of_hosts>8){
    S_matrix_list[[9]] <- out$S_host9_mat
    I_matrix_list[[9]] <- out$I_host9_mat
    if (number_of_hosts>9){
    S_matrix_list[[10]] <- out$S_host10_mat
    I_matrix_list[[10]] <- out$I_host10_mat
    }}}}}}}}}}
    
    ## CALCULATE OUTPUT TO PLOT:
    #I_host1_rast[] <- I_matrix_list[[1]]
    #I_host2_rast[] <- I_matrix_list[[2]]
    
    # 1) values as % infected
    #I_host2_rast[] <- ifelse(I_host2_rast[] == 0, NA, I_host2_rast[]/host2_rast[])
    
    # 2) values as number of infected per cell
    #I_host2_rast[] <- ifelse(I_host2_rast[] == 0, NA, I_host2_rast[])
    #I_host1_rast[] <- ifelse(I_host1_rast[] == 0, NA, I_host1_rast[])
    
    # 3) values as 0 (non infected) and 1 (infected) cell
    #I_host2_rast[] <- ifelse(I_host2_rast[] > 0, 1, 0) 
    #I_host2_rast[] <- ifelse(I_host2_rast[] > 0, 1, NA)
    
      
    if (cnt %in% yearlyoutputlist){
      yearTracker = yearTracker+1
      ## This is a set of nested if Statements
      if (number_of_hosts>0){
        I_host1_rast[] <- I_matrix_list[[1]]
        I_host1_rast[] <- ifelse(I_host1_rast[] == 0, NA, I_host1_rast[])
        I_host1_stack <- stack(I_host1_rast, I_host1_stack)
        dataForOutput$infectedHost1Individuals[yearTracker] <- sum(na.omit(I_host1_rast@data@values))/1000
        dataForOutput$infectedHost1Area[yearTracker] <- ncell(na.omit(I_host1_rast@data@values))*res_area
      if (number_of_hosts>1){
        I_host2_rast[] <- I_matrix_list[[2]]
        I_host2_rast[] <- ifelse(I_host2_rast[] == 0, NA, I_host2_rast[])
        I_host2_stack <- stack(I_host2_rast, I_host2_stack)
        dataForOutput$infectedHost2Individuals[yearTracker] <- sum(na.omit(I_host2_rast@data@values))/1000
        dataForOutput$infectedHost2Area[yearTracker] <- ncell(na.omit(I_host2_rast@data@values))*res_area
      if (number_of_hosts>2){
        I_host3_rast[] <- I_matrix_list[[3]]
        I_host3_rast[] <- ifelse(I_host3_rast[] == 0, NA, I_host3_rast[])
        I_host3_stack <- stack(I_host3_rast, I_host3_stack)
        dataForOutput$infectedHost3Individuals[yearTracker] <- sum(na.omit(I_host3_rast@data@values))/1000
        dataForOutput$infectedHost3Area[yearTracker] <- ncell(na.omit(I_host3_rast@data@values))*res_area
      if (number_of_hosts>3){
        I_host4_rast[] <- I_matrix_list[[4]]
        I_host4_rast[] <- ifelse(I_host4_rast[] == 0, NA, I_host4_rast[])
        I_host4_stack <- stack(I_host4_rast, I_host4_stack)
        dataForOutput$infectedHost4Individuals[yearTracker] <- sum(na.omit(I_host4_rast@data@values))/1000
        dataForOutput$infectedHost4Area[yearTracker] <- ncell(na.omit(I_host4_rast@data@values))*res_area
      if (number_of_hosts>4){
        I_host5_rast[] <- I_matrix_list[[5]]
        I_host5_rast[] <- ifelse(I_host5_rast[] == 0, NA, I_host5_rast[])
        I_host5_stack <- stack(I_host5_rast, I_host5_stack)
        dataForOutput$infectedHost5Individuals[yearTracker] <- sum(na.omit(I_host5_rast@data@values))/1000
        dataForOutput$infectedHost5Area[yearTracker] <- ncell(na.omit(I_host5_rast@data@values))*res_area
      if (number_of_hosts>5){
        I_host6_rast[] <- I_matrix_list[[6]]
        I_host6_rast[] <- ifelse(I_host6_rast[] == 0, NA, I_host6_rast[])
        I_host6_stack <- stack(I_host6_rast, I_host6_stack)
        dataForOutput$infectedHost6Individuals[yearTracker] <- sum(na.omit(I_host6_rast@data@values))/1000
        dataForOutput$infectedHost6Area[yearTracker] <- ncell(na.omit(I_host6_rast@data@values))*res_area
      if (number_of_hosts>6){
        I_host7_rast[] <- I_matrix_list[[7]]
        I_host7_rast[] <- ifelse(I_host7_rast[] == 0, NA, I_host7_rast[])
        I_host7_stack <- stack(I_host7_rast, I_host7_stack)
        dataForOutput$infectedHost7Individuals[yearTracker] <- sum(na.omit(I_host7_rast@data@values))/1000
        dataForOutput$infectedHost7Area[yearTracker] <- ncell(na.omit(I_host7_rast@data@values))*res_area 
      if (number_of_hosts>7){
        I_host8_rast[] <- I_matrix_list[[8]]
        I_host8_rast[] <- ifelse(I_host8_rast[] == 0, NA, I_host8_rast[])
        I_host8_stack <- stack(I_host8_rast, I_host8_stack)
        dataForOutput$infectedHost8Individuals[yearTracker] <- sum(na.omit(I_host8_rast@data@values))/1000
        dataForOutput$infectedHost8Area[yearTracker] <- ncell(na.omit(I_host8_rast@data@values))*res_area 
      if (number_of_hosts>8){
        I_host9_rast[] <- I_matrix_list[[9]]
        I_host9_rast[] <- ifelse(I_host9_rast[] == 0, NA, I_host9_rast[])
        I_host9_stack <- stack(I_host9_rast, I_host9_stack)
        dataForOutput$infectedHost9Individuals[yearTracker] <- sum(na.omit(I_host9_rast@data@values))/1000
        dataForOutput$infectedHost9Area[yearTracker] <- ncell(na.omit(I_host9_rast@data@values))*res_area 
      if (number_of_hosts>9){
        I_host10_rast[] <- I_matrix_list[[10]]
        I_host10_rast[] <- ifelse(I_host10_rast[] == 0, NA, I_host10_rast[])
        I_host10_stack <- stack(I_host10_rast, I_host10_stack)
        dataForOutput$infectedHost10Individuals[yearTracker] <- sum(na.omit(I_host10_rast@data@values))/1000
        dataForOutput$infectedHost10Area[yearTracker] <- ncell(na.omit(I_host10_rast@data@values))*res_area 
                        }}}}}}}}}}

    }
  
}

## switch the order so that it is from start of simulation to end and label the bands and add to output list
if (number_of_hosts>0){
  I_host1_stack <- subset(I_host1_stack, order(seq(nlayers(I_host1_stack)-1, 1, -1)))
  names(I_host1_stack) <- years
  I_totalhost_stack <- I_host1_stack
  data <- list(dataForOutput, I_totalhost_stack, I_host1_stack)
if (number_of_hosts>1){
  I_host2_stack <- subset(I_host2_stack, order(seq(nlayers(I_host2_stack)-1, 1, -1)))
  names(I_host2_stack) <- years
  data[[2]] <- data[[2]]+I_host2_stack
  data[[4]] <- I_host2_stack
if (number_of_hosts>2){
  I_host3_stack <- subset(I_host3_stack, order(seq(nlayers(I_host3_stack)-1, 1, -1)))
  names(I_host3_stack) <- years
  data[[2]] <- data[[2]]+I_host3_stack
  data[[5]] <- I_host3_stack
if (number_of_hosts>3){
  I_host4_stack <- subset(I_host4_stack, order(seq(nlayers(I_host4_stack)-1, 1, -1)))
  names(I_host4_stack) <- years
  data[[2]] <- data[[2]]+I_host4_stack
  data[[6]] <- I_host4_stack
if (number_of_hosts>4){
  I_host5_stack <- subset(I_host5_stack, order(seq(nlayers(I_host5_stack)-1, 1, -1)))
  names(I_host5_stack) <- years
  data[[2]] <- data[[2]]+I_host5_stack
  data[[7]] <- I_host5_stack
if (number_of_hosts>5){
  I_host6_stack <- subset(I_host6_stack, order(seq(nlayers(I_host6_stack)-1, 1, -1)))
  names(I_host6_stack) <- years
  data[[2]] <- data[[2]]+I_host6_stack
  data[[8]] <- I_host6_stack
if (number_of_hosts>6){
  I_host7_stack <- subset(I_host7_stack, order(seq(nlayers(I_host7_stack)-1, 1, -1)))
  names(I_host7_stack) <- years
  data[[2]] <- data[[2]]+I_host7_stack
  data[[9]] <- I_host7_stack
if (number_of_hosts>7){
  I_host8_stack <- subset(I_host8_stack, order(seq(nlayers(I_host8_stack)-1, 1, -1)))
  names(I_host8_stack) <- years
  data[[2]] <- data[[2]]+I_host8_stack
  data[[10]] <- I_host8_stack
if (number_of_hosts>8){
  I_host9_stack <- subset(I_host9_stack, order(seq(nlayers(I_host9_stack)-1, 1, -1)))
  names(I_host9_stack) <- years
  data[[2]] <- data[[2]]+I_host9_stack
  data[[11]] <- I_host9_stack
if (number_of_hosts>9){
  I_host10_stack <- subset(I_host10_stack, order(seq(nlayers(I_host10_stack)-1, 1, -1)))
  names(I_host10_stack) <- years
  data[[2]] <- data[[2]]+I_host10_stack
  data[[12]] <- I_host10_stack
}}}}}}}}}}

return(data)
}



