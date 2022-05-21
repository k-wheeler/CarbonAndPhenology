library(suncalc)
source('metFunctions.R')
source('sharedVariables.R')

##Read in and calculate daily met for full bbc1 (HWD tower) PhenoCam time-series (fall only) ----
dates <- seq(as.Date('2016-01-01'),as.Date('2020-12-31'),'day')
dates <- dates[format(dates,"%j")%in%seq(213,335)]
dailyMet <- calculateDailyMet(dates)
write.csv(dailyMet,file="dailyMetPhenoCam.csv",quote=FALSE,row.names = FALSE)

##Read in and calculate daily met for full harvard (EMS tower) PhenoCam time-series (fall only) ----
dates <- seq(as.Date('2010-01-01'),as.Date('2020-12-31'),'day')
dates <- dates[format(dates,"%j")%in%seq(213,335)]
latEMS=42.5378
longEMS=-72.1715
dailyMet <- calculateDailyMet(dates,isEMS=TRUE,lat = latEMS,long = longEMS)
write.csv(dailyMet,file="dailyMetPhenoCam_harvardEMS.csv",quote=FALSE,row.names = FALSE)

#Read in and calculate met data for fall field season only ----
dates <- seq(as.Date('2021-08-01'),as.Date('2021-12-01'),'day')
dailyMet <- calculateDailyMet(dates)
write.csv(dailyMet,file="dailyMetFieldSeasonOnly.csv",quote=FALSE,row.names = FALSE)

#Calculate Vertical Profile of Temperature ----
allMat <- cbind(dailyMet$date,dailyMet$airt.ac,dailyMet$airt.mid,dailyMet$airt.low)
metHeightsM <- metHeights[1:3]

pdfName <- "dailyMetMeasurements.pdf"
pdf(file=pdfName,height=8,width=10)
par(mfcol=c(2,2))
dtes <- seq(as.Date("2021-08-01"),as.Date("2021-12-01"),"day")
fittedMet <- as.data.frame(matrix(ncol=4,nrow=length(dtes)))
colnames(fittedMet) <- c("date","intercept","slope","R2")
fittedMet$date <- rep(Sys.Date(),nrow(fittedMet))
for(d in 1:length(dtes)){
  dte <- dtes[d]
  #dte <- as.Date("2021-10-01")
  dayMat <- allMat[allMat[,1]==dte,2:4]
  plot(metHeightsM[1:3],dayMat,pch=20,main=(dte),xlim=c(0,3000),ylab="Temperature",xlab="Height (m)")
  
  mdl <- lm(dayMat~metHeightsM)
  sm <- summary(mdl)
  abline(mdl,col="red")
  text(500,dayMat[3]+0.2,paste("R2:",round(sm$r.squared,digits=3)))
  fittedMet$date[d] <- dte
  fittedMet$R2[d] <- sm$r.squared
  fittedMet$slope[d] <- mdl$coefficients[2]
  fittedMet$intercept[d] <- mdl$coefficients[1]
}
plot(density(fittedMet$R2),main="Density of R^2s")
plot(fittedMet$date,fittedMet$R2,pch=20,xlab="Date",ylab="R^2",main="R^2 Over Time")
plot(fittedMet$date,fittedMet$slope,pch=20,xlab="Date",ylab="Slope",main="Slope Over Time")
plot(fittedMet$date,fittedMet$intercept,pch=20,xlab="Date",ylab="Intercept",main="Intercept Over Time")
dev.off()
write.csv(fittedMet,file="DailyDaytimeVerticalTemperatureProfiles.csv",quote=FALSE,row.names = FALSE)

#Calculate Vertical Profile of PAR ----
allMat <- cbind(dailyMet$date,dailyMet$parac,dailyMet$parmid,dailyMet$parlow,dailyMet$parus1)
metHeightsM <- metHeights[1:4]

pdfName <- "dailyMetMeasurements_par.pdf"
pdf(file=pdfName,height=8,width=10)
par(mfcol=c(2,2))
dtes <- seq(as.Date("2021-08-01"),as.Date("2021-12-01"),"day")
fittedMet <- as.data.frame(matrix(ncol=4,nrow=length(dtes)))
colnames(fittedMet) <- c("date","intercept","slope","R2")
fittedMet$date <- rep(Sys.Date(),nrow(fittedMet))
for(d in 1:length(dtes)){
  dte <- dtes[d]
  #dte <- as.Date("2021-10-01")
  dayMat <- allMat[allMat[,1]==dte,2:5]
  plot(metHeightsM[1:4],log(dayMat),pch=20,main=(dte),xlim=c(0,3000),ylab="Log(PAR)",xlab="Height (m)")
  
  mdl <- lm(log(dayMat)~metHeightsM)
  sm <- summary(mdl)
  abline(mdl,col="red")
  text(500,log(dayMat[3])+0.2,paste("R2:",round(sm$r.squared,digits=3)))
  fittedMet$date[d] <- dte
  fittedMet$R2[d] <- sm$r.squared
  fittedMet$slope[d] <- mdl$coefficients[2]
  fittedMet$intercept[d] <- mdl$coefficients[1]
  plot(metHeightsM[1:4],(dayMat),pch=20,main=(dte),xlim=c(0,3000),ylab="PAR",xlab="Height (m)")
  lines(seq(0,3000),exp(fittedMet$slope[d] *seq(0,3000)+fittedMet$intercept[d]))
  
}
plot(density(fittedMet$R2),main="Density of R^2s")
plot(fittedMet$date,fittedMet$R2,pch=20,xlab="Date",ylab="R^2",main="R^2 Over Time")
plot(fittedMet$date,fittedMet$slope,pch=20,xlab="Date",ylab="Slope",main="Slope Over Time")
plot(fittedMet$date,fittedMet$intercept,pch=20,xlab="Date",ylab="Intercept",main="Intercept Over Time")
dev.off()
write.csv(fittedMet,file="DailyDaytimeVerticalParProfiles.csv",quote=FALSE,row.names = FALSE)

#Calculate Vertical Profile of Relative Humidity ----
allMat <- cbind(dailyMet$date,dailyMet$rh.ac,dailyMet$rh.mid,dailyMet$rh.low)
metHeightsM <- metHeights[1:3]

pdfName <- "dailyMetMeasurements_rh.pdf"
pdf(file=pdfName,height=8,width=10)
par(mfcol=c(2,2))
dtes <- seq(as.Date("2021-08-01"),as.Date("2021-12-01"),"day")
fittedMet <- as.data.frame(matrix(ncol=6,nrow=length(dtes)))
colnames(fittedMet) <- c("date","intercept","slope","R2","average","var")
fittedMet$date <- rep(Sys.Date(),nrow(fittedMet))
for(d in 1:length(dtes)){
  dte <- dtes[d]
  #dte <- as.Date("2021-10-01")
  dayMat <- allMat[allMat[,1]==dte,2:4]
  plot(metHeightsM[1:3],dayMat,pch=20,main=(dte),xlim=c(0,3000),ylab="Daily Average RH (%)",xlab="Height (cm)")
  
  mdl <- lm(dayMat~metHeightsM)
  sm <- summary(mdl)
  abline(mdl,col="red")
  text(500,dayMat[1]+0.2,paste("R2:",round(sm$r.squared,digits=3)))
  fittedMet$date[d] <- dte
  fittedMet$R2[d] <- sm$r.squared
  fittedMet$slope[d] <- mdl$coefficients[2]
  fittedMet$intercept[d] <- mdl$coefficients[1]
  fittedMet$average[d] <- mean(dayMat,na.rm=TRUE)
  fittedMet$var[d] <- var(dayMat,na.rm=TRUE)
  
}
plot(density(fittedMet$R2),main="Density of R^2s")
plot(fittedMet$date,fittedMet$R2,pch=20,xlab="Date",ylab="R^2",main="R^2 Over Time")
plot(fittedMet$date,fittedMet$slope,pch=20,xlab="Date",ylab="Slope",main="Slope Over Time")
plot(fittedMet$date,fittedMet$intercept,pch=20,xlab="Date",ylab="Intercept",main="Intercept Over Time")
plot(fittedMet$var,fittedMet$R2,pch=20,xlab="Variance in RH",ylab="R^2",main="R^2 vs. Variance in RH")
abline(h=0.8,col="red",lwd=2)
plot(fittedMet$average,fittedMet$R2,pch=20,xlab="Average RH",ylab="R^2",main="R^2 vs. Average RH")
dev.off()
write.csv(fittedMet,file="DailyDaytimeVerticalRHProfiles.csv",quote=FALSE,row.names = FALSE)

