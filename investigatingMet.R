library(suncalc)
lat <- 42.5351
long <- -72.1744
#Read in met data
dates <- seq(as.Date('20-01'),as.Date('2021-12-31'),'day')
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

metHeights <- c(28,18,9,1)*100 #Heights in mm of met sensor locations
metHeightNames <- c('airt.ac','airt.mid','airt.low','airt.us1')
metDat <- read.csv('Data/hf282-01-hdwd-tower.csv',stringsAsFactors = FALSE)

subMet <- metDat[lubridate::date(metDat$datetime)%in%dates,
                 c('datetime','airt.ac','airt.mid','airt.low')]
subMet$datetime <- stringr::str_replace(subMet$datetime,"T"," ")
subMet$datetime <- as.POSIXct(subMet$datetime,format='%Y-%m-%d %H:%M')
subMet$date <- lubridate::date(subMet$datetime)
subMet$airt.ac <- as.numeric(subMet$airt.ac)
subMet$airt.mid <- as.numeric(subMet$airt.mid)
subMet$airt.low <- as.numeric(subMet$airt.low)
dailyMet <- as.data.frame(matrix(nrow=length(dates),ncol=5))
colnames(dailyMet) <- c('date','airt.ac','airt.mid','airt.low','airt.us1')
dailyMet$date <- dates
#subMet$date <- as.Date(subMet$date)
splitMet <- split(subMet,subMet$date)
splitDates <- as.Date(names(splitMet))
for(d in 1:length(dates)){
  ind <- which(splitDates==dates[d])
  if(length(ind)>0){
    dateMet <- splitMet[ind][[1]]
    dayMet <- dateMet[dateMet$datetime>sunRises[d]&dateMet$datetime<sunSets[d],]
    dailyMet$airt.ac[d] <- mean(dayMet$airt.ac,na.rm=TRUE)
    dailyMet$airt.mid[d] <- mean(dayMet$airt.mid,na.rm=TRUE)
    dailyMet$airt.low[d] <- mean(dayMet$airt.low,na.rm=TRUE)
  }
}

metDat2 <- read.csv('Data/hf282-02-hdwd-tower-understory.csv',stringsAsFactors = FALSE)
subMet2 <- metDat2[lubridate::date(metDat2$datetime)%in%dates,
                   c('datetime','airt.us1')]
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
    dailyMet$airt.us1[d] <- mean(dayMet$airt.us1)
  } 
}

load("Data/bbc1_dataFinal_includeJuly.RData")
#dataFinal$TairMuDay[,dataFinal$N]

plot(dates,dailyMet$airt.ac,pch=20,ylim=c(-10,30),col="#0868ac",cex=2,xlab="Date",ylab="Tair (C)",main="Average Day Time Temperature")
points(dates,dailyMet$airt.mid,pch=20,col=scales::alpha("#43a2ca"),cex=2)
points(dates,dailyMet$airt.low,pch=20,col=scales::alpha("#7bccc4"),cex=2)
points(dates,dailyMet$airt.us1,pch=20,col=scales::alpha("#bae4bc"),cex=2)
#points(dates,dataFinal$TairMuDay[,dataFinal$N],pch=20,cex=2,col="red")
legend('bottomleft',legend=metHeights,col=c("#0868ac","#43a2ca","#7bccc4","#bae4bc"),pch=rep(20,4),cex=1.5)

plot(dates,dailyMet$airt.ac-dailyMet$airt.mid,pch=20)
abline(h=0,col="gray",lty=2)

plot(dates,dailyMet$airt.mid-dailyMet$airt.low,pch=20)
abline(h=0,col="gray",lty=2)

plot(dates,dailyMet$airt.low-dailyMet$airt.us1,pch=20)
abline(h=0,col="gray",lty=2)

plot(dates,dailyMet$airt.ac-dailyMet$airt.us1,pch=20)
abline(h=0,col="gray",lty=2)
abline(v=as.Date("2019-01-01"),col="red")

allDat <- merge(subMet,subMet2,'datetime')
plot(allDat$datetime,allDat$airt.ac-allDat$airt.us1,pch=20,cex=0.25,xlab="Time",ylab="Difference",
     main="Difference in Temperature Between Levels")
#abline(h=0,col="gray",lty=2,lwd=2)

points(allDat$datetime,allDat$airt.ac-allDat$airt.low,pch=20,cex=0.25,col=scales::alpha("cyan",0.5))
abline(h=0,col="red",lty=2,lwd=3)
legend('bottomleft',legend=c("28m - 1m","28m - 9m"),col=c("black","cyan"),pch=rep(20,2))
