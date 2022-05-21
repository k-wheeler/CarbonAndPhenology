##General Functions
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

thinMCMCOutput <- function(out.burn){
  if(typeof(out.burn)!=typeof(FALSE)){
    out.mat <- as.matrix(out.burn$params)
    thinAmount <- round(nrow(out.mat)/5000,digits=0)
    out.burn2 <- list()
    out.burn2$params <- window(out.burn$params,thin=thinAmount)
    out.burn2$predict <- window(out.burn$predict,thin=thinAmount)
    out.burn <- out.burn2
  }
  return(out.burn)
}

calculateDayLengths <- function(dates,lat=42.5351, long = -72.1744){
  dayLengths <- numeric()
  sunRises <- rep(Sys.time(),length(dates))
  sunSets <- rep(Sys.time(),length(dates))
  for(d in 1:length(dates)){
    suntimes <- getSunlightTimes(date=dates[d],
                                 lat=lat,lon=long,keep=c("nauticalDawn","nauticalDusk"),
                                 tz = "GMT") #GMT because I only care about difference
    dayLengths <- c(dayLengths,as.numeric(suntimes$nauticalDusk-suntimes$nauticalDawn))
    sunRises[d] <- (suntimes$nauticalDawn)
    sunSets[d] <- (suntimes$nauticalDusk)
  }
  return(cbind(dayLengths,sunRises,sunSets))
}


Licor_QAQC <- function(files,dat){
  if(files %in% c("Data/Licor_Measurements/HarvardForest_Licor_Wheeler-06SEP21-B1A",
                  "Data/Licor_Measurements/HarvardForest_Licor_Wheeler-20AUG21-B2B")){
    dat <- dat[dat$PARi<1490,]
  }else if(files=="Data/Licor_Measurements/HarvardForest_Licor_Wheeler-10AUG21-B1A"){
    dat <- dat[-c(16,17),]
  }else if(files=="Data/Licor_Measurements/HarvardForest_Licor_Wheeler-17SEP21-B3B"){
    dat <- dat[-14,]
  }else if(files=="Data/Licor_Measurements/HarvardForest_Licor_Wheeler-20AUG21-B2E"){
    dat <- dat[-5,]
  }else if(files%in% c("Data/Licor_Measurements/HarvardForest_Licor_Wheeler-06SEP21-O4B",
                       "Data/Licor_Measurements/HarvardForest_Licor_Wheeler-12AUG21-O2B",
                       "Data/Licor_Measurements/HarvardForest_Licor_Wheeler-18SEP21-O4B")){
    dat <- dat[-(1:nrow(dat)),]
  }else if(files=="Data/Licor_Measurements/HarvardForest_Licor_Wheeler-07SEP21-O4B"){
    dat <- dat[-13,]
  }else if(files=="Data/Licor_Measurements/HarvardForest_Licor_Wheeler-08SEP21-O2A"){
    dat <- dat[-13,]
  }else if(files=="Data/Licor_Measurements/HarvardForest_Licor_Wheeler-19SEP21-O1A"){ #Very high photosynthesis values...
    dat <- dat[-15,]
  }
  return(dat)
}
