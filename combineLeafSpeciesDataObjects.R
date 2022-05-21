library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(doParallel)
library(abind)
source('generalModels.R')
source('sharedVariables.R')
source('generalFunctions.R')
source('metFunctions.R')
options(stringsAsFactors=FALSE)
species <- "oak"

if(species=="beech"){
  trees <- c("B1","B2","B3")
}else if(species=="oak"){
  trees <- c("O1","O2","O3","O4")
}
lNames <- c("A","B","E")

for(t in 1:length(trees)){
  tName <- trees[t]
  print(tName)
  for(l in 1:3){
    lfName <- paste0(tName,lNames[l])
    print(lfName)
    if(file.exists(paste0('Data/finalData/',lfName,"_finalData_withTran_withPhoto.RData"))){
    load(paste0('Data/finalData/',lfName,"_finalData_withTran_withPhoto.RData"))
    }else{
      load(paste0('Data/finalData/',lfName,"_finalData_withTran.RData"))
    }
    if(l==1){
      CCI_means <- finalData$CCI_means
      CCI_precs<- finalData$CCI_precs
      D <- finalData$D
      
      if(is.null(finalData$GPP)){
        Tair <- rep(NA,(finalData$n+1))
        tran <- NA
        NPP <- rep(NA,(finalData$n+1))
      }else{
        Tair <- finalData$Tair
        tran <- finalData$tran[4]
        NPP <- finalData$NPP
        GPP <- finalData$GPP
        R <- finalData$R
        aj <- finalData$aj
        ac <- finalData$ac
      }
    }else{
      CCI_means <- rbind(CCI_means,finalData$CCI_means)
      CCI_precs<- rbind(CCI_precs,finalData$CCI_precs)
      if(is.null(finalData$GPP)){
        Tair <- rbind(Tair,rep(NA,(finalData$n+1)))
        tran <- c(tran,NA)
        NPP <- rbind(NPP,rep(NA,(finalData$n)))
        GPP <- rbind(GPP,rep(NA,(finalData$n)))
        R <- rbind(R,rep(NA,(finalData$n)))
        aj <- rbind(aj,rep(NA,(finalData$n)))
        ac <- rbind(ac,rep(NA,(finalData$n)))
        
      }else{
        Tair <- rbind(Tair,finalData$Tair)
        tran <- c(tran,finalData$tran[4])
        NPP <- rbind(NPP,finalData$NPP)
        GPP <- rbind(GPP,finalData$GPP)
        R <- rbind(R,finalData$R)
        aj <- rbind(aj,finalData$aj)
        ac <- rbind(ac,finalData$ac)
      }
    }
    if(l==3 && t==1){
      allData <- list(CCI_means=CCI_means,
                      CCI_precs=CCI_precs,
                      NPP=NPP,
                      GPP=GPP,
                      aj=aj,
                      ac=ac,
                      Resp=R,
                      tran=tran,
                      n=finalData$n,
                      D=finalData$D,
                      Tair=Tair,
                      dates=finalData$dates,
                      treeN=length(trees))

    }else if(l==3){
      allData$CCI_means <- abind(allData$CCI_means,CCI_means,along=3)
      allData$CCI_precs <- abind(allData$CCI_precs,CCI_precs,along=3)
      allData$Tair <- abind(allData$Tair,Tair,along=3)
      allData$GPP <- abind(allData$GPP,GPP,along=3)
      allData$NPP <- abind(allData$NPP,NPP,along=3)
      allData$aj <- abind(allData$aj,aj,along=3)
      allData$ac <- abind(allData$ac,ac,along=3)
      allData$R <- abind(allData$R,R,along=3)
      allData$tran <- cbind(allData$tran,tran)
    }
  }
}

allData$N <- nrow(allData$Tair)
allData$CCI_precs[is.infinite(allData$CCI_precs)] <- getmode(allData$CCI_precs)
allData$CCI_precs[allData$CCI_precs>50000] <- getmode(allData$CCI_precs)
allData$CCI_means[allData$CCI_means==1] <- 0.999
allData$CCI_means[allData$CCI_means==0] <- 0.001
allData$n <- 122

save(file=paste0('Data/finalData/',"all",species,"LicorLeafData.RData"),allData)
