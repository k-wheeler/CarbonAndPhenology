calculateDailyMet <- function(dates,isEMS=FALSE,lat=42.5351,long= -72.1744){
  dayLengthMat <- calculateDayLengths(dates,lat=lat,long = long)
  dayLengths <- dayLengthMat[,1]
  sunRises <- dayLengthMat[,2]
  sunSets <- dayLengthMat[,3]
  
  if(isEMS){
    metDat <- read.csv('Data/hf004-01-final.csv',stringsAsFactors = FALSE)
    
    subMet <- metDat[lubridate::date(metDat$datetime)%in%dates,
                     c('datetime','ta.27.9m',"par.29",
                       "rh.27.9m")]
    subMet$datetime <- stringr::str_replace(subMet$datetime,"T"," ")
    subMet$datetime <- as.POSIXct(subMet$datetime,format='%Y-%m-%d %H:%M')
    subMet$date <- lubridate::date(subMet$datetime)
    subMet$ta.27.9m <- as.numeric(subMet$ta.27.9m)
    subMet$par.29 <- as.numeric(subMet$par.29)
    subMet$rh.27.9m<- as.numeric(subMet$rh.27.9m)
   
    dailyMet <- as.data.frame(matrix(nrow=length(dates),ncol=6))
    colnames(dailyMet) <- c('date','ta.27.9m',"par.29",
                              "rh.27.9m","co2","co2_sd")
    dailyMet$date <- dates
    #subMet$date <- as.Date(subMet$date)
    splitMet <- split(subMet,subMet$date)
    splitDates <- as.Date(names(splitMet))
    for(d in 1:length(dates)){
      ind <- which(splitDates==dates[d])
      if(length(ind)>0){
        dateMet <- splitMet[ind][[1]]
        dayMet <- dateMet[dateMet$datetime>sunRises[d]&dateMet$datetime<sunSets[d],]
        dailyMet$ta.27.9m[ind] <- mean(dayMet$ta.27.9m,na.rm=TRUE)
        dailyMet$par.29[ind] <- mean(dayMet$par.29,na.rm=TRUE)
        dailyMet$rh.27.9m[ind] <- mean(dayMet$rh.27.9m,na.rm=TRUE)
      }
    }
  }else{
  metHeights <- c(28,18,9,1)*100 #Heights in mm of met sensor locations
  metHeightNames <- c('airt.ac','airt.mid','airt.low','airt.us1')
  metDat <- read.csv('Data/hf282-01-hdwd-tower.csv',stringsAsFactors = FALSE)
  
  subMet <- metDat[lubridate::date(metDat$datetime)%in%dates,
                   c('datetime','airt.ac','airt.mid','airt.low',"parac.up","parmid","parlow",
                     "rh.ac","rh.mid","rh.low")]
  subMet$datetime <- stringr::str_replace(subMet$datetime,"T"," ")
  subMet$datetime <- as.POSIXct(subMet$datetime,format='%Y-%m-%d %H:%M')
  subMet$date <- lubridate::date(subMet$datetime)
  subMet$airt.ac <- as.numeric(subMet$airt.ac)
  subMet$airt.mid <- as.numeric(subMet$airt.mid)
  subMet$airt.low <- as.numeric(subMet$airt.low)
  subMet$parac.up <- as.numeric(subMet$parac.up)
  subMet$parmid <- as.numeric(subMet$parmid)
  subMet$parlow <- as.numeric(subMet$parlow)
  subMet$rh.ac <- as.numeric(subMet$rh.ac)
  subMet$rh.mid <- as.numeric(subMet$rh.mid)
  subMet$rh.low <- as.numeric(subMet$rh.low)
  dailyMet <- as.data.frame(matrix(nrow=length(dates),ncol=14))
  colnames(dailyMet) <- c('date','airt.ac','airt.mid','airt.low','airt.us1',
                          "parac","parmid","parlow",'parus1',"rh.ac","rh.mid","rh.low","co2","co2_sd")
  dailyMet$date <- dates
  #subMet$date <- as.Date(subMet$date)
  splitMet <- split(subMet,subMet$date)
  splitDates <- as.Date(names(splitMet))
  for(d in 1:length(dates)){
    ind <- which(splitDates==dates[d])
    if(length(ind)>0){
      dateMet <- splitMet[ind][[1]]
      dayMet <- dateMet[dateMet$datetime>sunRises[d]&dateMet$datetime<sunSets[d],]
      dailyMet$airt.ac[ind] <- mean(dayMet$airt.ac,na.rm=TRUE)
      dailyMet$airt.mid[ind] <- mean(dayMet$airt.mid,na.rm=TRUE)
      dailyMet$airt.low[ind] <- mean(dayMet$airt.low,na.rm=TRUE)
      dailyMet$parac[ind] <- mean(dayMet$parac.up,na.rm = TRUE)
      dailyMet$parmid[ind] <- mean(dayMet$parmid,na.rm = TRUE)
      dailyMet$parlow[ind] <- mean(dayMet$parlow,na.rm = TRUE)
      dailyMet$rh.ac[ind] <- mean(dayMet$rh.ac,na.rm=TRUE)
      dailyMet$rh.mid[ind] <- mean(dayMet$rh.mid,na.rm=TRUE)
      dailyMet$rh.low[ind] <- mean(dayMet$rh.low,na.rm=TRUE)
    }
  }
  
  metDat2 <- read.csv('Data/hf282-02-hdwd-tower-understory.csv',stringsAsFactors = FALSE)
  subMet2 <- metDat2[lubridate::date(metDat2$datetime)%in%dates,
                     c('datetime','airt.us1','parus1')]
  subMet2$date <- lubridate::date(subMet2$datetime)
  subMet2$datetime <- stringr::str_replace(subMet2$datetime,"T"," ")
  subMet2$datetime <- as.POSIXct(subMet2$datetime,format='%Y-%m-%d %H:%M')
  splitMet <- split(subMet2,subMet2$date)
  splitDates <- as.Date(names(splitMet))
  for(d in 1:length(dates)){
    ind <- which(splitDates==dates[d])
    if(length(ind)>0){
      dateMet <- splitMet[ind][[1]]
      dayMet <- dateMet[dateMet$datetime>sunRises[d]&dateMet$datetime<sunSets[d],]
      dailyMet$airt.us1[ind] <- mean(dayMet$airt.us1,na.rm=TRUE)
      dailyMet$parus1[ind] <- mean(dayMet$parus1,na.rm=TRUE)
    } 
  }
  }
  
  co2Met <- read.table("Data/co2_mlo_surface-insitu_1_ccgg_DailyData.txt",header=TRUE,skip = 150,sep=" ")
  co2Select <- numeric()
  co2SelectSD <- numeric()
  for(d in 1:length(dates)){
    subDat <- co2Met[co2Met$month==lubridate::month(dates[d]) & co2Met$day==lubridate::day(dates[d]) &co2Met$year==lubridate::year(dates[d])
                     ,c('value','value_std_dev',"qcflag")]
    if(substr(subDat$qcflag,1,1)!="."){
      co2Select <- c(co2Select,NA)
      co2SelectSD <- c(co2SelectSD,NA)
    }else{
      co2Select <- c(co2Select,subDat$value)
      co2SelectSD <- c(co2SelectSD,subDat$value_std_dev)
    }
  }
  while(sum(is.na(co2Select))>0){
    co2Select[is.na(co2Select)] <- co2Select[which(is.na(co2Select))-1]
  }
  while(sum(is.na(co2SelectSD))>0){
    co2SelectSD[is.na(co2SelectSD)] <- co2SelectSD[which(is.na(co2SelectSD))-1]
  }
  dailyMet$co2 <- co2Select
  dailyMet$co2_sd <- co2SelectSD
  return(dailyMet)
}



#From: https://github.com/PecanProject/pecan/blob/ee38544e1c995f47073019bd31b3179e7e64991a/modules/data.atmosphere/R/metutils.R

##' Calculate saturation vapor pressure
##'
##' @title get es
##' @param temp temperature in degrees C
##' @return saturation vapor pressure in mb
##' @export
##' @author David LeBauer
##' @examples
##' temp <- -30:30
##' plot(temp, get.es(temp))
get.es <- function(temp) {
  return(6.11 * exp((2500000/461) * (1/273 - 1/(273 + temp))))
} # get.es

##' Calculate VPD
##'
##' Calculate vapor pressure deficit from relative humidity and temperature.
##' @title VPD
##' @param rh relative humidity, in percent
##' @param temp temperature, degrees celsius
##' @return vpd: vapor pressure deficit, in mb
##' @export
##' @author David LeBauer
##' @examples
##' temp <- -30:30
##' plot(temp, get.vpd(0, temp))
get.vpd <- function(rh, temp) {
  ## calculate saturation vapor pressure
  es <- get.es(temp)
  ## calculate vapor pressure deficit
  return(((100 - rh)/100) * es)
} # get.vpd
