#Find CCI measurements for each Licor File 
source('read_Licor.R')
source('fitA.R')
source('runLicorIter.R')
source('sharedVariables.R')
source('generalFunctions.R')
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(PEcAn.photosynthesis)#From: https://github.com/PecanProject/pecan/tree/develop/modules/photosynthesis/R

Allfiles <- paste0(LicorDataDirectory,dir(path=LicorDataDirectory,pattern="-"))
AllfilesClean <- dir(path=LicorDataDirectory,pattern="-")

CCIoutput <- matrix(nrow=0,ncol=4)
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
    CCIoutput <- rbind(CCIoutput,
                    cbind(paste0(tree,lfNums),rep(dte,length(means)),means,sds))
  }
}

beechRescaled <- as.data.frame(subset(CCIoutput,substr(CCIoutput[,1],1,1)=="B"))
colnames(beechRescaled) <- c("Leaf","Date","CCI_mean","CCI_sd")
beechRescaled[,3:4] <- scales::rescale(c(as.numeric(as.character(beechRescaled$CCI_mean)),
                                    as.numeric(as.character(beechRescaled$CCI_sd))),to=c(0.0001,0.9999))

oakRescaled <- as.data.frame(subset(CCIoutput,substr(CCIoutput[,1],1,1)=="O"))
colnames(oakRescaled) <- c("Leaf","Date","CCI_mean","CCI_sd")
oakRescaled[,3:4] <- scales::rescale(c(as.numeric(as.character(oakRescaled$CCI_mean)),
                                         as.numeric(as.character(oakRescaled$CCI_sd))),to=c(0.0001,0.9999))
allRescaled <- rbind(beechRescaled,oakRescaled)

allRescaled$Date <-   dte <- as.Date(allRescaled$Date ,format="%m/%d/%Y")
LicorMatchingValues <- as.data.frame(matrix(ncol=3,nrow=length(Allfiles)))
colnames(LicorMatchingValues) <- c("fname","CCI_mean","CCI_sd")
LicorMatchingValues$fname <- AllfilesClean
for(f in 1:length(Allfiles)){
  lfName <- strsplit(Allfiles[f],"-")[[1]][3]
  dt <- strsplit(Allfiles[f],"-")[[1]][2]
  if(substr(dt,3,3)=="S"){
    mth <- "09"
  }else if(substr(dt,3,3)=="A"){
    mth <- "08"
  }
  dy <- substr(dt,1,2)
  dte <- paste(mth,dy,"2021",sep="/")
  dte <- as.Date(dte,format="%m/%d/%Y")
  if(length(which(allRescaled$Leaf==lfName & allRescaled$Date==dte))>0){
    LicorMatchingValues$CCI_mean[f]<- allRescaled[allRescaled$Leaf==lfName & allRescaled$Date==dte,'CCI_mean']
    LicorMatchingValues$CCI_sd[f]<- allRescaled[allRescaled$Leaf==lfName & allRescaled$Date==dte,'CCI_sd']
  }else{
    subDat <- subset(allRescaled,allRescaled$Leaf==lfName)
    LicorMatchingValues$CCI_mean[f]<- subDat[which.min(abs(subDat$Date-dte)),'CCI_mean']
    LicorMatchingValues$CCI_sd[f]<- subDat[which.min(abs(subDat$Date-dte)),'CCI_sd']
  }
}
write.csv(LicorMatchingValues,quote = FALSE,row.names = FALSE,file="CCIforLicorMeasurements.csv")
