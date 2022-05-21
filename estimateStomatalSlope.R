#Estimate Stomatal Slope
source('generalModels.R')
source('sharedVariables.R')
source('generalFunctions.R')
library(PEcAn.photosynthesis)
#From: https://github.com/PecanProject/pecan/tree/develop/modules/photosynthesis/R

#Allfiles <- paste0(LicorDataDirectorytory,dir(path=LicorDataDirectorytory,pattern='-B'))
Allfiles <- paste0(LicorDataDirectorytory,intersect(dir(path=LicorDataDirectorytory,pattern='-O'),dir(path=LicorDataDirectorytory,pattern='AUG')))


for(i in 1:length(Allfiles)){
  print(i)
  files <- Allfiles[i]
  master = lapply(files, read_Licor)
  dat <- do.call("rbind", master)
  dat <- Licor_QAQC(files,dat)
  if(i==1){
    allDat <- dat
  }else{
    allDat <- rbind(allDat,dat)
  }
}

medlynX <- allDat$Photo/(allDat$CO2S*sqrt(allDat$VpdL))

plot(medlynX[medlynX>0.0015],allDat$Cond[medlynX>0.0015],pch=20)
abline(lm(allDat$Cond[medlynX>0.0015]~medlynX[medlynX>0.0015]),col="red")
sm <- summary(lm(allDat$Cond[medlynX>0.0015]~medlynX[medlynX>0.0015]))
sm$coefficients[2] #Beech Stomatal Slope
#3.737588 for stomatal slope for beech
#4.85734 for oak
sm$coefficients[1] #Beech Stomatal g0
#0.05215698 for stomatal g0 for beech
#0.1852357 for oak
