#Organize Data
library(suncalc)
library('dplyr')
source('metFunctions.R')

##Organize Leaf-Level Data ----
heightData <- read.csv('Data/leafHeights.csv',stringsAsFactors = FALSE)
heightData$Height <- heightData$Height*2.54 #Convert from inches to cm
leafNames <- heightData$LeafName[order(heightData$Height,decreasing = FALSE)]
output <- matrix(nrow=0,ncol=4)
for(t in 1:length(trees)){
  tree <- trees[t]
  files <- dir(path=CCIdataDirectory,pattern=tree)
  for(f in 1:length(files)){
    dat <- read.csv(paste0(CCIdataDirectory,files[f]),stringsAsFactors = FALSE,header=TRUE)
    dte <- strsplit(as.character(dat$Time.Date[2])," ")[[1]][2]
    sds <- numeric()
    means <- numeric()
    vls <- numeric()
    lfNums <- character()
    for(i in 1:nrow(dat)){
      #for(i in 1:5){
      if(as.character(dat$Units[i])==" CCI"){
        sds <- c(sds,sd(vls))
        means <- c(means,mean(vls,na.rm=TRUE))
        vls <- numeric()
      }else if((is.na(dat$Sample[i]) | dat$Sample[i]=="")){
        vls <- c(vls,as.numeric(dat$Reading[i]))
      }else if(substr(dat$Sample[i],1,1)%in%c("A","B","C","D","E","F","M")){
        lfNums <- c(lfNums,substr(dat$Sample[i],1,1))
      }
    }
      output <- rbind(output,
                      cbind(paste0(tree,lfNums),rep(dte,length(means)),means,sds))
  }
}
output <- as.data.frame(output)
colnames(output) <- c("Leaf","Date","CCI_mean","CCI_sd")
output$Date <- as.Date(output$Date,format="%m/%d/%Y")

dates <- seq(as.Date('2021-08-01'),as.Date('2021-11-30'),'day')

dayLengthMat <- calculateDayLengths(dates)
dayLengths <- dayLengthMat[,1]
sunRises <- dayLengthMat[,2]
sunSets <- dayLengthMat[,3]
fittedMet <- read.csv("DailyDaytimeVerticalTemperatureProfiles.csv")
fittedParMet <- read.csv("DailyDaytimeVerticalParProfiles.csv")
fittedRHMet <- read.csv('DailyDaytimeVerticalRHProfiles.csv')
dailyMet <- read.csv("dailyMetFallOnly.csv",header=TRUE)

for(l in 1:length(leafNames)){
  lfName <- leafNames[l]
  print(lfName)
  finalMat <- as.data.frame(matrix(ncol=1,nrow=length(dates)))
  colnames(finalMat) <- c("Date")
  finalMat$Date <- dates
  subDat <- subset(output,Leaf==lfName)
  subDat <- subDat[!duplicated(subDat[,2]),]
  finalMat <- merge(finalMat,subDat,"Date",all=TRUE)
  finalMat[,3:4] <- scales::rescale(c(as.numeric(as.character(finalMat$CCI_mean)),
                                      as.numeric(as.character(finalMat$CCI_sd))),to=c(0.0001,0.9999))
  finalMat$CCI_sd[is.na(finalMat$CCI_sd)] <- mean(finalMat$CCI_sd,na.rm=TRUE)
  finalData <- list(leaf=lfName,CCI_means=as.numeric(as.character(finalMat$CCI_mean)),
                    CCI_precs=1/(as.numeric(as.character(finalMat$CCI_sd)))**2)
  finalData$CCI_precs[is.infinite(finalData$CCI_precs)] <- getmode(finalData$CCI_precs)
  finalData$CCI_precs[finalData$CCI_precs>50000] <- getmode(finalData$CCI_precs)
  finalData$n <- length(dates)
  finalData$dates <- dates
  finalData$D <- dayLengths
  finalData$height <- heightData[heightData$LeafName==lfName,2]
  finalData$Tair <- fittedMet$intercept + finalData$height * fittedMet$slope
  finalData$par <- exp(fittedParMet$intercept + finalData$height * fittedParMet$slope)
  finalData$RH <- fittedRHMet$intercept + finalData$height * fittedRHMet$slope
  finalData$RH[fittedRHMet$R2<0.9] <- fittedRHMet$average[fittedRHMet$R2<0.9]
  finalData$vpd <- get.vpd(finalData$RH,finalData$Tair)*0.1 #Convert to kPa
  for(d in 1:finalData$n){
    if(finalData$par[d]<=0){
      finalData$par[d] <- 0.001
    }
    if(finalData$vpd[d]<=0){
      finalData$vpd[d] <- 0.001
    }
  }
  finalData$CCI_means[seq(which.min(finalData$CCI_means),finalData$n)] <- NA
  finalData$co2 <- dailyMet$co2
  finalData$co2_sd <- dailyMet$co2_sd
  
  save(finalData,file=paste0('Data/finalData/',lfName,"_finalData.RData"))
  #plot(dates,finalData$CCI_means,pch=20)
}

#Organize PhenoCam Data ----
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(suncalc)
library(rnoaa)
library(doParallel)
library(ncdf4)
source('/projectnb/dietzelab/kiwheel/chlorophyllCycling/load_ERA5.R')

dataDirectory <- "/projectnb/dietzelab/kiwheel/chlorophyllCycling/data/"
siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)
s=6 #For bbc1 site
s=17 #For harvard EMS site
siteName <- as.character(siteData$siteName[s])
print(siteName)
lat <- as.numeric(siteData[s,2])
long <- as.numeric(siteData[s,3])
startDate <- (as.Date(siteData[s,7]))
library(lubridate)
if(s==17){
  startDate <- startDate %m+% lubridate::years(1) #Not enough met data for 2009
}

endDate <- as.character(siteData$endDate[s])
URL <- as.character(siteData$URL[s])
URL2 <- as.character(siteData$URL2[s])
URL3 <- as.character(siteData$URL3[s])
if(!is.na(URL2)){
  URL <- c(URL,URL2)
  if(!is.na(URL3)){
    URL <- c(URL,URL3)
  }
}
TZ <- as.numeric(siteData[s,6])
URLs <- URL
load(file=paste(dataDirectory,siteName,"_phenopixOutputs.RData",sep=""))
fittedDat=allDat
ERA5dataFolder <- paste("/projectnb/dietzelab/kiwheel/ERA5/Data/",siteName,"/",sep="")

phenoData <- matrix(nrow=0,ncol=32)
print(URLs[1])
for(u in 1:length(URLs)){
  phenoDataSub <- download.phenocam(URLs[u])
  phenoData <- rbind(phenoData,phenoDataSub)
}

##Order and remove duplicate PC data
phenoData2 <- phenoData[order(phenoData$date),]
phenoData3 <- phenoData2[!duplicated(phenoData2$date),]
phenoData <- phenoData3

phenoData <- phenoData[phenoData$date<endDate,]
p.old <- phenoData$gcc_90
time.old <-  as.Date(phenoData$date)
days <- seq(as.Date(startDate),(as.Date(endDate)),"day")
p <- rep(NA,length(days))

for(i in 1:length(p.old)){
  p[which(days==time.old[i])] <- p.old[i]
}

months <- lubridate::month(days)
years <- lubridate::year(days)

dat2 <- data.frame(dates=days,years=years,months=months,p=p)
# calFileName <- paste0(siteName,"_",startDate,"_",endDate,"_era5TemperatureMembers.nc")
# datTairEns <- load_ERA5(ERA5dataFolder=ERA5dataFolder,calFileName=calFileName,TZ_offset=TZ,variable="Tair")
# datTairEnsDay <- load_ERA5_daytime(ERA5dataFolder=ERA5dataFolder,calFileName=calFileName,TZ_offset=TZ,variable="Tair",lat=lat,long=long)
# 
# TairMu <- apply(X=datTairEns,MARGIN=2,FUN=mean)
# TairPrec <- 1/apply(X=datTairEns,MARGIN=2,FUN=var)
# dat2$TairMu <- TairMu 
# dat2$TairPrec <- TairPrec
# 
# TairMuDay <- apply(X=datTairEnsDay,MARGIN=2,FUN=mean)
# TairPrecDay <- 1/apply(X=datTairEnsDay,MARGIN=2,FUN=var)
# dat2$TairMuDay <- TairMuDay 
# dat2$TairPrecDay <- TairPrecDay

dayLengths <- numeric()

for(d in 1:length(days)){
  suntimes <- getSunlightTimes(date=days[d],
                               lat=lat,lon=long,keep=c("nauticalDawn","nauticalDusk"),
                               tz = "GMT") #GMT because I only care about difference
  dayLengths <- c(dayLengths,as.numeric(suntimes$nauticalDusk-suntimes$nauticalDawn))
}

dat2$D <- dayLengths
ICsdat <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(203,212),] 
dat2 <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(213,335),] #Starting July 1st

#nrowNum <- 365-212
# nrowNum <- 365-181
nrowNum <- 335-212
p <- matrix(nrow=nrowNum,ncol=0)
# TairMu <- matrix(nrow=nrowNum,ncol=0)
# TairMuDay <- matrix(nrow=nrowNum,ncol=0)
D <- matrix(nrow=nrowNum,ncol=0)
ICs <- matrix(nrow=10,ncol=0)
# TairPrec <- matrix(nrow=nrowNum,ncol=0)
# TairPrecDay <- matrix(nrow=nrowNum,ncol=0)
valNum <- 0
days2 <- matrix(nrow=nrowNum,ncol=0)

finalYrs <- numeric()
sofs <- numeric()
for(i in (lubridate::year(as.Date(dat2$dates[1]))):lubridate::year(as.Date(dat2$dates[length(dat2$dates)]))){
  subDat <- dat2[lubridate::year(as.Date(dat2$dates))==i,]
  #valNum <- valNum + 1
  valNum <- which(fittedDat[,'Year']==i)
  if(length(valNum)==0){
    Low <- NA
    High <- NA
  }else{
    Low <- fittedDat[valNum,'Low']
    High <- fittedDat[valNum,'High']
  }
  if(!is.na(Low)){
    newICs <- scales::rescale(ICsdat[lubridate::year(as.Date(ICsdat$dates))==i,]$p,from=c(Low,High))
    if(length(na.omit(newICs))>5){
      newCol <- scales::rescale(subDat$p,to=c(0,1),from=c(Low,High))
      p <- cbind(p,newCol)
      ICs <- cbind(ICs,newICs)
      days2 <- cbind(days2,as.Date(subDat$dates))
      finalYrs <- c(finalYrs,i)
      #sofs <- c(sofs,(fittedDat[valNum,'FallStartDay']-212))
      sofs <- c(sofs,(fittedDat[valNum,'FallStartDay']-212)) ######Change for start if needed
      # TairMu <- cbind(TairMu,subDat$TairMu)
      # TairMuDay <- cbind(TairMuDay,subDat$TairMuDay)
      D <- cbind(D,subDat$D)
      # TairPrec <- cbind(TairPrec,subDat$TairPrec)
      # TairPrecDay <- cbind(TairPrecDay,subDat$TairPrecDay)
    }
  }
}
p[p<0] <- 0
p[p>0.999] <- 0.999
ICs[ICs<0] <- 0
ICs[ICs>0.999] <- 0.999

dataFinal <- list(p=p,years=finalYrs,sofMean=mean(sofs))
dataFinal$n <- nrowNum
dataFinal$N <- ncol(dataFinal$p)


x1a <- numeric()
x1b <- numeric()
for(yr in 1:dataFinal$N){
  mu <- mean(ICs[,yr],na.rm=TRUE)
  vr <- var(ICs[,yr],na.rm = TRUE)
  x1a <- c(x1a,(mu**2-mu**3-mu*vr)/(vr))
  x1b <- c(x1b,(mu-2*mu**2+mu**3-vr+mu*vr)/(vr))
}

dataFinal$x1.a <- x1a
dataFinal$x1.b <- x1b

# dataFinal$TairMu <- TairMu
# dataFinal$TairPrec <- TairPrec
# dataFinal$TairMuDay <- TairMuDay
# dataFinal$TairPrecDay <- TairPrecDay
dataFinal$D <- D

dataFinal$TowerTair <- matrix(nrow=dataFinal$n,ncol=dataFinal$N)
dataFinal$co2 <- matrix(nrow=dataFinal$n,ncol=dataFinal$N)
dataFinal$par <- matrix(nrow=dataFinal$n,ncol=dataFinal$N)
dataFinal$RH <- matrix(nrow=dataFinal$n,ncol=dataFinal$N)

if(s==6){
  dailyMet <- read.csv("dailyMetPhenoCam.csv")
  for(yr in 1:dataFinal$N){
    dataFinal$TowerTair[,yr] <- dailyMet[lubridate::year(dailyMet$date)==dataFinal$years[yr],'airt.ac'][1:dataFinal$n]
    dataFinal$co2[,yr] <- dailyMet[lubridate::year(dailyMet$date)==dataFinal$years[yr],'co2'][1:dataFinal$n]
    dataFinal$par[,yr] <- dailyMet[lubridate::year(dailyMet$date)==dataFinal$years[yr],'parac'][1:dataFinal$n]
    dataFinal$RH[,yr] <- dailyMet[lubridate::year(dailyMet$date)==dataFinal$years[yr],'rh.ac'][1:dataFinal$n]
    
  }
}else if(s==17){
  dailyMet <- read.csv("dailyMetPhenoCam_harvardEMS.csv")
  for(yr in 1:dataFinal$N){
    dataFinal$TowerTair[,yr] <- dailyMet[lubridate::year(dailyMet$date)==dataFinal$years[yr],'ta.27.9m'][1:dataFinal$n]
    dataFinal$co2[,yr] <- dailyMet[lubridate::year(dailyMet$date)==dataFinal$years[yr],'co2'][1:dataFinal$n]
    dataFinal$par[,yr] <- dailyMet[lubridate::year(dailyMet$date)==dataFinal$years[yr],'par.29'][1:dataFinal$n]
    dataFinal$RH[,yr] <- dailyMet[lubridate::year(dailyMet$date)==dataFinal$years[yr],'rh.27.9m'][1:dataFinal$n]
    
  }
}

while(sum(is.na(dataFinal$TowerTair))>0){
  dataFinal$TowerTair[is.na(dataFinal$TowerTair)] <- dataFinal$TowerTair[which(is.na(dataFinal$TowerTair))-1]
}
while(sum(is.na(dataFinal$RH))>0){
  dataFinal$RH[is.na(dataFinal$RH)] <- dataFinal$RH[which(is.na(dataFinal$RH))-1]
}

while(sum(is.na(dataFinal$par))>0){
  dataFinal$par[is.na(dataFinal$par)] <- dataFinal$par[which(is.na(dataFinal$par))-1]
}
dataFinal$vpd <- get.vpd(dataFinal$RH,dataFinal$TowerTair)*0.1 #Convert to kPa
save(dataFinal,file=paste0(finalDataDirectory,siteName,"_dataFinal.RData"))



               