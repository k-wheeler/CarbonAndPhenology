#Create Data Object:
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(suncalc)
lat <- 42.5351
long <- -72.1744
dates <- seq(as.Date('2016-01-01'),as.Date('2020-12-31'),'day')
dates <- dates[format(dates,"%j")%in%seq(182,365)]
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
fileDates <- unique(subMet$date)
#subMet$date <- as.Date(subMet$date)
splitMet <- split(subMet,subMet$date)
for(d in 1:length(fileDates)){
  dateMet <- splitMet[d][[1]]
  datesInd <- which(dates==fileDates[d])
  dayMet <- dateMet[dateMet$datetime>sunRises[datesInd]&dateMet$datetime<sunSets[datesInd],]
  dailyMet$airt.ac[datesInd] <- mean(dayMet$airt.ac,na.rm=TRUE)
}
#plot(dates,dailyMet$airt.ac,pch=20)
load("Data/bbc1_dataFinal_includeJuly.RData")
dataFinal$TowerTair <- matrix(nrow=dataFinal$n,ncol=dataFinal$N)
splitDaily <- split(dailyMet,lubridate::year(dailyMet$date))
for(yr in 1:dataFinal$N){
  dataFinal$TowerTair[,yr] <- splitDaily[yr][[1]][,'airt.ac']
}
dataFinal$TowerTair[is.na(dataFinal$TowerTair)] <- dataFinal$TairMuDay[is.na(dataFinal$TowerTair)]

generalModel = "
model {
    ### Data Models for complete years
    for(yr in 1:(N)){
    for(i in 1:n){
    p[i,yr] ~ dnorm(x[i,yr],p.PC)
    }
    }
    
    #### Process Model
    for(yr in 1:(N)){
    for(i in 2:n){
    
    xmu[i,yr] <- max(min((x[(i-1),yr] + b4 * x[(i-1),yr]) + max(0,(b0 + b3 * TowerTair[i,yr] * D[i,yr])),x[1,yr]),0)
    x[i,yr] ~ dnorm(xmu[i,yr],p.proc)
    }
    }
    
    #### Priors
    for(yr in 1:N){ ##Initial Conditions
    x[1,yr] ~ dbeta(x1.a[yr],x1.b[yr]) I(0.001,0.999)
    }
    p.PC ~ dgamma(s1.PC,s2.PC)
    p.proc ~ dgamma(s1.proc,s2.proc)
    b0 ~ dunif(b0_lower,b0_upper)
    b3 ~ dunif(b3_lower,b3_upper)
    b4 ~ dunif(b4_lower,b4_upper)
    
  }
    "
variableNames <- c("x","p.proc","b0","b3","b4","p.PC")

nchain=5

outputFileName <- paste0("modelFits/bbc1_fullAutumn_dayTime_varBurn.RData")
partialFileName <- paste0("modelFits/bbc1_fullAutumn_dayTime_partial_varBurn.RData")
dataFinal$s1.PC <- 1.56
dataFinal$s2.PC <- 0.016
dataFinal$s1.proc <- 1.56
dataFinal$s2.proc <- 0.016
dataFinal$b0_lower <- -0.25
dataFinal$b0_upper <- 0
dataFinal$b3_lower <- 0
dataFinal$b3_upper <- 0.25
dataFinal$b4_lower <- -0.25
dataFinal$b4_upper <- 0

j.model <- try(jags.model(file = textConnection(generalModel),
                          data = dataFinal,
                          n.chains = nchain,
                          n.adapt = 3000))#Load Model Output 

out.burn <- try(runForecastIter(j.model=j.model,variableNames=variableNames,
                                baseNum = 15000,iterSize = 5000,effSize = 5000, maxIter=1000000,
                                partialFile = partialFileName))

##Thin the data:
if(typeof(out.burn)!=typeof(FALSE)){
  out.mat <- as.matrix(out.burn$params)
  thinAmount <- round(nrow(out.mat)/5000,digits=0)
  out.burn2 <- list()
  out.burn2$params <- window(out.burn$params,thin=thinAmount)
  out.burn2$predict <- window(out.burn$predict,thin=thinAmount)
  out.burn <- out.burn2
}
save(out.burn,file = outputFileName)
print(paste("saved:",outputFileName))
