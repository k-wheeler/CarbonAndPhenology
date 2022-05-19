#Organize Data
library(suncalc)
library('dplyr')

leafNames <- c("B1A","B1B","B1C","B1D","B1E","B1F",
               "B2A","B2B","B2C","B2D","B2E","B2F",
               "B3A","B3B","B3C","B3D","B3E","B3F",
               "O1A","O1B","O1C","O1D","O1E","O1F",
               "O2A","O2B","O2C","O2D","O2E","O2F",
               "O3A","O3B","O3C","O3D","O3E","O3F",
               "O4A","O4B","O4C","O4D","O4E","O4F")
trees <- c("B1","B2","B3","O1","O2","O3","O4")
letterSequence <- c("a","b","c","d","e","f")
dataDirectory <- "Data/CCI_Measurements/"
heightData <- read.csv('Data/leafHeights.csv',stringsAsFactors = FALSE)
heightData$Height <- heightData$Height*2.54 #Convert from inches to cm
leafNames <- heightData$LeafName[order(heightData$Height,decreasing = FALSE)]
output <- matrix(nrow=0,ncol=4)
for(t in 1:length(trees)){
  tree <- trees[t]
  #print(tree)
  files <- dir(path=dataDirectory,pattern=tree)
  for(f in 1:length(files)){
    dat <- read.csv(paste0(dataDirectory,files[f]),stringsAsFactors = FALSE,header=TRUE)
    dte <- strsplit(as.character(dat$Time.Date[2])," ")[[1]][2]
    # if(is.na(dte)){
    #   print(tree)
    #   print(f)
    # }
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
colfunc <- colorRampPalette(c("black", "white"))
cols <- c(colfunc(25)[1:18],colfunc(30)[1:24])
cols2 <- c(rep(cols[1],3),rep(cols[4],3),rep(cols[7],3),rep(cols[10],3),
           rep(cols[13],3),rep(cols[16],3),rep(cols[19],3),rep(cols[22],3),
           rep(cols[25],3),rep(cols[28],3),rep(cols[31],3),rep(cols[33],3),
           rep(cols[36],3),rep(cols[39],3))
#test <- output[order(output$Date),]
jpeg("leafCCI_data.jpeg",width=7,height=6.4,units = "in",res=1000)
par(mfrow=c(2,1))
par(mai=c(0.8,1,0.5,0.1))
subDat <- subset(output,Leaf==leafNames[19])
leafHeight <- heightData$Height[heightData$LeafName==leafNames[19]]
plot(subDat$Date,as.numeric(as.character(subDat$CCI_mean)),type = "l",
     ylim=c(0,25),xlim=range(output$Date),main="Oak Leaves (light is higher)",ylab="CCI",xlab="",bty="n",col=cols2[19])
for(i in 20:length(leafNames)){
  subDat <- subset(output,Leaf==leafNames[i])
  leafHeight <- heightData$Height[heightData$LeafName==leafNames[i]]
  lines(subDat$Date,as.numeric(as.character(subDat$CCI_mean)),col=cols2[i])
}

subDat <- subset(output,Leaf==leafNames[1])
leafHeight <- heightData$Height[heightData$LeafName==leafNames[1]]
plot(subDat$Date,as.numeric(as.character(subDat$CCI_mean)),type = "l",
     ylim=c(0,25),xlim=range(output$Date),
     xlab="Date", ylab="CCI",bty="n",main="Beech Leaves (light is higher)",col=cols2[1])
for(i in 2:18){
  subDat <- subset(output,Leaf==leafNames[i])
  leafHeight <- heightData$Height[heightData$LeafName==leafNames[19]]
  lines(subDat$Date,as.numeric(as.character(subDat$CCI_mean)),col=cols2[i])
}

dev.off()
dates <- seq(as.Date('2021-08-01'),as.Date('2021-11-30'),'day')
lat <- 42.5351
long <- -72.1744
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
#Read in met data
metHeights <- c(28,18,9,1)*100 #Heights in mm of met sensor locations

metHeightNames <- c('airt.ac','airt.mid','airt.low','airt.us1')
metDat <- read.csv('Data/hf282-01-hdwd-tower.csv',stringsAsFactors = FALSE)
#metDat$datetime <- as.POSIXlt(metDat$datetime)
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

# dailyMet$airt.ac<- sapply(split(subMet[,2],subMet$date),mean)
# dailyMet$airt.mid<- sapply(split(subMet[,3],subMet$date),mean)
# dailyMet$airt.low<- sapply(split(subMet[,4],subMet$date),mean)
# dailyMet$airt.ac[is.na(dailyMet$airt.ac)] <- dailyMet$airt.ac[which(is.na(dailyMet$airt.ac))-1]
# dailyMet$airt.mid[is.na(dailyMet$airt.mid)] <- dailyMet$airt.mid[which(is.na(dailyMet$airt.mid))-1]
# dailyMet$airt.low[is.na(dailyMet$airt.low)] <- dailyMet$airt.low[which(is.na(dailyMet$airt.low))-1]

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
# dailyMet$airt.us1<- sapply(split(subMet2[,2],subMet2$date),mean)
dailyMet$airt.us1[is.na(dailyMet$airt.us1)] <- dailyMet$airt.us1[which(is.na(dailyMet$airt.us1))-1]
dailyMet$airt.us1[is.na(dailyMet$airt.us1)] <- dailyMet$airt.us1[which(is.na(dailyMet$airt.us1))-1]
dailyMet$airt.us1[is.na(dailyMet$airt.us1)] <- dailyMet$airt.us1[which(is.na(dailyMet$airt.us1))-1]
dailyMet$airt.us1[is.na(dailyMet$airt.us1)] <- dailyMet$airt.us1[which(is.na(dailyMet$airt.us1))-1]

plot(dates,dailyMet$airt.ac,pch=20,ylim=c(-10,30),col="#0868ac",cex=2,xlab="Date",ylab="Tair (C)",main="Average Day Time Temperature")
points(dates,dailyMet$airt.mid,pch=20,col="#43a2ca",cex=2)
points(dates,dailyMet$airt.low,pch=20,col="#7bccc4",cex=2)
points(dates,dailyMet$airt.us1,pch=20,col="#bae4bc",cex=2)
legend('bottomleft',legend=metHeights,col=c("#0868ac","#43a2ca","#7bccc4","#bae4bc"),pch=rep(20,4),cex=1.5)

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
                                      as.numeric(as.character(finalMat$CCI_sd))),to=c(0,1))
  finalMat$CCI_sd[is.na(finalMat$CCI_sd)] <- mean(finalMat$CCI_sd,na.rm=TRUE)
  finalData <- list(leaf=lfName,CCI_means=as.numeric(as.character(finalMat$CCI_mean)),
                    CCI_precs=1/(as.numeric(as.character(finalMat$CCI_sd)))**2)
  
  finalData$n <- length(dates)
  finalData$dates <- dates
  finalData$D <- dayLengths
  finalData$height <- heightData[heightData$LeafName==lfName,2]
  lfMetName <- metHeightNames[which.min(abs(finalData$height - metHeights))]
  finalData$Tair <- dailyMet[,lfMetName]
  finalData$CCI_means[seq(which.min(finalData$CCI_means),finalData$n)] <- NA
  
  save(finalData,file=paste0('Data/finalData/',lfName,"_finalData.RData"))
  #plot(dates,finalData$CCI_means,pch=20)
}





               