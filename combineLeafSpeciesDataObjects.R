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
species <- "beech"
photoOnly <- TRUE
feedback <- TRUE

if(species=="beech"){
  trees <- c("B1","B2","B3")
}else if(species=="oak"){
  trees <- c("O1","O2","O3","O4")
}
if(photoOnly){
  lNames <- c("A","B","E")
}else{
  if(species=="beech"){
    lNames <- c("B1A","B1B","B1C","B1D","B1E","B1F",
                "B2A","B2B","B2C","B2D","B2E","B2F",
                "B3A","B3B","B3C","B3D","B3E","B3F")
  }else if(species=="oak"){
    lNames <- c("O1A","O1B","O1C","O1D","O1E","O1F",
                "O2A","O2B","O2C","O2D","O2E","O2F",
                "O3A","O3B","O3C","O3D","O3E","O3F",
                "O4A","O4B","O4C","O4D","O4E","O4F")
  }
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
if(photoOnly){
  if(feedback){
    for(t in 1:length(trees)){
      tName <- trees[t]
      print(tName)
      for(l in 1:3){
        lfName <- paste0(tName,lNames[l])
        print(lfName)
        if(file.exists(paste0('Data/finalData/',lfName,"_finalData_withTran_withPhotoFeedback.RData"))){
          load(paste0('Data/finalData/',lfName,"_finalData_withTran_withPhotoFeedback.RData"))
        }else{
          load(paste0('Data/finalData/',lfName,"_finalData_withTran.RData"))
        }
        if(l==1){
          CCI_means <- finalData$CCI_means
          CCI_precs<- finalData$CCI_precs
          D <- finalData$D
          
          if(is.null(finalData$GPP)){
            Tair <- rep(NA,(finalData$n))
            NPP <- rep(NA,(finalData$n))
            vmax0 <- NA
            Jmax0 <- NA
            r0 <- NA
            
          }else{
            Tair <- finalData$Tair[1:122]
            NPP <- finalData$NPP[1:122,]
            GPP <- finalData$GPP[1:122,]
            R <- finalData$R[1:122,]
            aj <- finalData$aj[1:122,]
            ac <- finalData$ac[1:122,]
            vmax0 <- finalData$vmax0
            Jmax0 <- finalData$Jmax0
            r0 <- finalData$r0
          }
          tran <- finalData$tran[4]
          height <- finalData$height
        }else{
          CCI_means <- rbind(CCI_means,finalData$CCI_means[1:122])
          CCI_precs<- rbind(CCI_precs,finalData$CCI_precs[1:122])
          if(is.null(finalData$GPP)){
            Tair <- rbind(Tair,rep(NA,(finalData$n)))
            NPP <- abind(NPP,rep.col(rep(NA,(finalData$n)),100),along=3)
            GPP <- abind(GPP,rep.col(rep(NA,(finalData$n)),100),along=3)
            R <- abind(R,rep.col(rep(NA,(finalData$n)),100),along=3)
            aj <- abind(aj,rep.col(rep(NA,(finalData$n)),100),along=3)
            ac <- abind(ac,rep.col(rep(NA,(finalData$n)),100),along=3)
            vmax0 <- c(vmax0,NA)
            Jmax0 <- c(Jmax0,NA)
            r0 <- c(r0,NA)
            
          }else{
            Tair <- rbind(Tair,finalData$Tair[1:122])
            NPP <- abind(NPP,finalData$NPP[1:122,],along=3)
            GPP <- abind(GPP,finalData$GPP[1:122,],along=3)
            R <- abind(R,finalData$R[1:122,],along=3)
            aj <- abind(aj,finalData$aj[1:122,],along=3)
            ac <- abind(ac,finalData$ac[1:122,],along=3)
            vmax0 <- c(vmax0,finalData$vmax0)
            Jmax0 <- c(Jmax0,finalData$Jmax0)
            r0 <- c(r0,finalData$r0)
          }
          tran <- c(tran,finalData$tran[4])
          height <- c(height,finalData$height)
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
                          treeN=length(trees),
                          vmax0=vmax0,
                          Jmax0=Jmax0,
                          r0=r0,
                          height=height)
          
        }else if(l==3){
          allData$CCI_means <- abind(allData$CCI_means,CCI_means[,1:122],along=3)
          allData$CCI_precs <- abind(allData$CCI_precs,CCI_precs[,1:122],along=3)
          allData$Tair <- abind(allData$Tair,Tair[,1:122],along=3)
          allData$GPP <- abind(allData$GPP,GPP[1:122,,],along=4)
          allData$NPP <- abind(allData$NPP,NPP[1:122,,],along=4)
          allData$aj <- abind(allData$aj,aj[1:122,,],along=4)
          allData$ac <- abind(allData$ac,ac[1:122,,],along=4)
          allData$R <- abind(allData$R,R[1:122,,],along=4)
          allData$tran <- cbind(allData$tran,tran)
          allData$vmax0 <- cbind(allData$vmax0,vmax0)
          allData$Jmax0 <- cbind(allData$Jmax0,Jmax0)
          allData$r0 <- cbind(allData$r0,r0)
          allData$height <- cbind(allData$height,height)
        }
      }
    }
  }else{
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
          NPP <- rep(NA,(finalData$n+1))
          vmax0 <- NA
          Jmax0 <- NA
          r0 <- NA
          
        }else{
          Tair <- finalData$Tair
          NPP <- finalData$NPP
          GPP <- finalData$GPP
          R <- finalData$R
          aj <- finalData$aj
          ac <- finalData$ac
          vmax0 <- finalData$vmax0
          Jmax0 <- finalData$Jmax0
          r0 <- finalData$r0
        }
        tran <- finalData$tran[4]
        height <- finalData$height
      }else{
        CCI_means <- rbind(CCI_means,finalData$CCI_means)
        CCI_precs<- rbind(CCI_precs,finalData$CCI_precs)
        if(is.null(finalData$GPP)){
          Tair <- rbind(Tair,rep(NA,(finalData$n+1)))
          NPP <- rbind(NPP,rep(NA,(finalData$n)))
          GPP <- rbind(GPP,rep(NA,(finalData$n)))
          R <- rbind(R,rep(NA,(finalData$n)))
          aj <- rbind(aj,rep(NA,(finalData$n)))
          ac <- rbind(ac,rep(NA,(finalData$n)))
          vmax0 <- c(vmax0,NA)
          Jmax0 <- c(Jmax0,NA)
          r0 <- c(r0,NA)
          
        }else{
          Tair <- rbind(Tair,finalData$Tair)
          NPP <- rbind(NPP,finalData$NPP)
          GPP <- rbind(GPP,finalData$GPP)
          R <- rbind(R,finalData$R)
          aj <- rbind(aj,finalData$aj)
          ac <- rbind(ac,finalData$ac)
          vmax0 <- c(vmax0,finalData$vmax0)
          Jmax0 <- c(Jmax0,finalData$Jmax0)
          r0 <- c(r0,finalData$r0)
        }
        tran <- c(tran,finalData$tran[4])
        height <- c(height,finalData$height)
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
                        treeN=length(trees),
                        vmax0=vmax0,
                        Jmax0=Jmax0,
                        r0=r0,
                        height=height)
        
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
        allData$vmax0 <- cbind(allData$vmax0,vmax0)
        allData$Jmax0 <- cbind(allData$Jmax0,Jmax0)
        allData$r0 <- cbind(allData$r0,r0)
        allData$height <- cbind(allData$height,height)
      }
    }
  }
  }
  allData$N <- nrow(allData$Tair)
}else{
  for(l in 1:length(lNames)){
    lfName <- lNames[l]
    print(lfName)
    load(paste0('Data/finalData/',lfName,"_finalData_withTran.RData"))
    
    if(l==1){
      CCI_means <- finalData$CCI_means
      CCI_precs<- finalData$CCI_precs
      D <- finalData$D
      Tair <- finalData$Tair
      
      tran <- finalData$tran[4]
      height <- finalData$height
    }else{
      CCI_means <- cbind(CCI_means,finalData$CCI_means)
      CCI_precs<- cbind(CCI_precs,finalData$CCI_precs)
      Tair <- cbind(Tair,finalData$Tair)
      
      tran <- c(tran,finalData$tran[4])
      height <- c(height,finalData$height)
    }
    if(l==length(lNames)){
      allData <- list(CCI_means=CCI_means,
                      CCI_precs=CCI_precs,
                      tran=tran,
                      n=finalData$n,
                      D=finalData$D,
                      Tair=Tair,
                      dates=finalData$dates,
                      N=length(lNames),
                      height=height)
      
    }
  }
}


allData$CCI_precs[is.infinite(allData$CCI_precs)] <- getmode(allData$CCI_precs)
allData$CCI_precs[allData$CCI_precs>50000] <- getmode(allData$CCI_precs)
allData$CCI_means[allData$CCI_means==1] <- 0.999
allData$CCI_means[allData$CCI_means==0] <- 0.001
allData$n <- 122
if(photoOnly){
  if(feedback){
    save(file=paste0('Data/finalData/',"all",species,"LicorLeafDataFeedback.RData"),allData)
  }else{
  save(file=paste0('Data/finalData/',"all",species,"LicorLeafData.RData"),allData)
  }
}else{
  save(file=paste0('Data/finalData/',"all_noPhoto_",species,"LicorLeafData.RData"),allData)
}

