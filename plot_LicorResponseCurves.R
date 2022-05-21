library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(doParallel)
#From: https://github.com/PecanProject/pecan/tree/develop/modules/photosynthesis/R

n.cores <- 8

#register the cores.
registerDoParallel(cores=n.cores)

dataDirectory <- "Data/Licor_Measurements/"
leafNames <- c("B1A","B1B","B1C","B1D","B1E","B1F",
               "B2A","B2B","B2C","B2D","B2E","B2F",
               "B3A","B3B","B3C","B3D","B3E","B3F",
               "O1A","O1B","O1C","O1D","O1E","O1F",
               "O2A","O2B","O2C","O2D","O2E","O2F",
               "O3A","O3B","O3C","O3D","O3E","O3F",
               "O4A","O4B","O4C","O4D","O4E","O4F")

pdf(file="FieldWorkResponseCurvesConverged_Aug.pdf",height=5,width=12)
par(mfrow=c(1,2))
for(l in length(leafNames):1){
  print(leafNames[l])
  lfName <- leafNames[l]
  load(paste0('Data/finalData/',lfName,"_finalData_withTran.RData"))
  #files <- dir(path=dataDirectory,pattern=leafNames[l])
  files <- intersect(dir(path=dataDirectory,pattern=leafNames[l]),dir(path=dataDirectory,pattern="AUG"))
  if(length(files)==0){
    files <- intersect(dir(path=dataDirectory,pattern=leafNames[l]),dir(path=dataDirectory,pattern="SEP"))
    if(length(files)==0){
      next
    }
  }
  f=1
  fileName <- paste0(dataDirectory,files[f])
  print(fileName)
  
  dteChar <- strsplit(files[f],"-")[[1]][2]
  if(substr(dteChar,3,3)=="A"){
    dte <- as.Date(paste0("2021-08-",substr(dteChar,1,2)))
  }else if(substr(dteChar,3,3)=="S"){
    dte <- as.Date(paste0("2021-09-",substr(dteChar,1,2)))
  }
  print(dte)
  if(file.exists(paste0('LicorFits/',leafNames[l],"_",dte,"_Cleaned_LicorResponseCurve_varBurn.RData"))){
    dat <- read_Licor(filename=fileName)
    
    if(fileName%in% c("Data/Licor_Measurements/HarvardForest_Licor_Wheeler-06SEP21-B1A",
                      "Data/Licor_Measurements/HarvardForest_Licor_Wheeler-20AUG21-B2B")){
      dat <- dat[dat$PARi<1490,]
    }else if(fileName =="Data/Licor_Measurements/HarvardForest_Licor_Wheeler-10AUG21-B1A"){
      dat <- dat[-c(16,17),]
    }else if(fileName =="Data/Licor_Measurements/HarvardForest_Licor_Wheeler-17SEP21-B3B"){
      dat <- dat[-14,]
    }else if(fileName =="Data/Licor_Measurements/HarvardForest_Licor_Wheeler-20AUG21-B2E"){
      dat <- dat[-5,]
    }else if(fileName %in% c("Data/Licor_Measurements/HarvardForest_Licor_Wheeler-06SEP21-O4B",
                             "Data/Licor_Measurements/HarvardForest_Licor_Wheeler-12AUG21-O2B",
                             "Data/Licor_Measurements/HarvardForest_Licor_Wheeler-18SEP21-O4B")){
      next
    }else if(fileName =="Data/Licor_Measurements/HarvardForest_Licor_Wheeler-07SEP21-O4B"){
      dat <- dat[-13,]
    }else if(fileName =="Data/Licor_Measurements/HarvardForest_Licor_Wheeler-08SEP21-O2A"){
      dat <- dat[-13,]
    }else if(fileName =="Data/Licor_Measurements/HarvardForest_Licor_Wheeler-19SEP21-O1A"){ #Very high photosynthesis values...
      dat <- dat[-15,]
    }
    load(paste0('LicorFits/',leafNames[l],"_",dte,"_Cleaned_LicorResponseCurve_varBurn.RData"))
    plot_photo2(dat,out,byLeaf=FALSE,id=paste(leafNames[l],dte))
  }
}
dev.off()




# source('read_Licor.R')
# source('fitA.R')
# source('runLicorIter.R')
# source('plot_photo.R')
# library(PhenoForecast)
# library(PhenologyBayesModeling)
# library(rjags)
# library(runjags)
# library(doParallel)
# library(ecoforecastR)
# 
# n.cores <- 8
# 
# #register the cores.
# registerDoParallel(cores=n.cores)
# 
# dataDirectory <- "Data/Licor_Measurements/"
# leafNames <- c("B1A","B1B","B1C","B1D","B1E","B1F",
#                "B2A","B2B","B2C","B2D","B2E","B2F",
#                "B3A","B3B","B3C","B3D","B3E","B3F",
#                "O1A","O1B","O1C","O1D","O1E","O1F",
#                "O2A","O2B","O2C","O2D","O2E","O2F",
#                "O3A","O3B","O3C","O3D","O3E","O3F",
#                "O4A","O4B","O4C","O4D","O4E","O4F")
# allPhoto <- matrix(ncol=4,nrow=0) #LeafName, dte, Vcmax, Jmax
# allDates <- rep(Sys.Date(),51)
# #allPhoto <- rbind(allPhoto,c("Test",as.Date(Sys.Date()),NA,NA))
# j=1
# pdf(file="FieldWorkResponseCurvesConverged.pdf",height=5,width=6)
# #par(mfrow=c(1,2))
# for(l in 1:length(leafNames)){
#   print(leafNames[l])
#   files <- dir(path=dataDirectory,pattern=leafNames[l])
#   jmax <- numeric()
#   vcmax <- numeric()
#   dates <- rep(Sys.Date(),length(files))
#   if(length(files)>0){
#     for(f in 1:length(files)){
#       fileName <- paste0(dataDirectory,files[f])
#       print(fileName)
#       dteChar <- strsplit(files[f],"-")[[1]][2]
#       if(substr(dteChar,3,3)=="A"){
#         dte <- as.Date(paste0("2021-08-",substr(dteChar,1,2)))
#       }else if(substr(dteChar,3,3)=="S"){
#         dte <- as.Date(paste0("2021-09-",substr(dteChar,1,2)))
#       }
#       print(dte)
#       dat <- read_Licor(filename=fileName)
#       if(file.exists(paste0('LicorFits/',leafNames[l],"_",dte,"_LicorResponseCurve_varBurn.RData"))){
#       load(file=paste0('LicorFits/',leafNames[l],"_",dte,"_LicorResponseCurve_varBurn.RData"))
#       jmax <- c(jmax,mean(data.frame(as.matrix(out$params))$Jmax0))
#       vcmax <- c(vcmax,mean(data.frame(as.matrix(out$params))$vmax0))
#       dates[f] <- dte
#       allPhoto <- rbind(allPhoto,c(leafNames[l],dte,mean(data.frame(as.matrix(out$params))$Jmax0),mean(data.frame(as.matrix(out$params))$vmax0)))
#       allDates[j] <- dte
#       j <- j+1
#       plot_photo(dat,out,byLeaf=FALSE,id=paste0(leafNames[l],dte))
#       # plot(dat$Ci,dat$Photo,pch=20,main=paste(leafNames[l],dte))
#       # plot(dat$PARi,dat$Photo,pch=20,main=paste(leafNames[l],dte))
#       }else{
#         jmax <- c(jmax,NA)
#         vcmax <- c(vcmax,NA)
#       }
#     }
#     if(length(na.omit(vcmax)>0)){
#     plot(dates,vcmax,main=paste(leafNames[l],"Vcmax"),pch=20)
#     plot(dates,jmax,main=paste(leafNames[l],"Jmax"),pch=20)
#     }
#   }
# }
# dev.off()


##Note: leafNames changes to sorted by height
jpeg("allPhotoParametersAgainstTime.jpeg",width=7,height=7,units = "in",res=1000)
par(mfrow=c(2,1))
par(mai=c(0.8,1,0.5,0.1))
plot(allDates,as.numeric(allPhoto[,3]),pch=20,cex=0.2,xlab="Time",ylab="Jmax",ylim=c(20,280),main="Lighter is Higher Height (By Species)")
for(l in 1:length(leafNames)){
  subDat <- allPhoto[allPhoto[,1]==leafNames[l],]
  subDts <- allDates[allPhoto[,1]==leafNames[l]]
  lines(subDts[order(subDts)],as.numeric(subDat[,3])[order(subDts)],col=cols[l])
}

for(i in 1:nrow(allPhoto)){
  if(substr(allPhoto[i,1],1,1)=="B"){
    points(allDates[i],allPhoto[i,3],pch=20,col=cols[which(leafNames==allPhoto[i,1])])
  }else{
    points(allDates[i],allPhoto[i,3],pch=17,col=cols[which(leafNames==allPhoto[i,1])])
  }
}
legend('topright',c("Beech","Oak"),pch=c(20,17))

plot(allDates,as.numeric(allPhoto[,4]),pch=20,cex=0.2,xlab="Time",ylab="Vcmax",ylim=c(20,140),main="Lighter is Higher Height(By Species)")
for(l in 1:length(leafNames)){
  subDat <- allPhoto[allPhoto[,1]==leafNames[l],]
  subDts <- allDates[allPhoto[,1]==leafNames[l]]
  lines(subDts[order(subDts)],as.numeric(subDat[,4])[order(subDts)],col=cols[l])
}
for(i in 1:nrow(allPhoto)){
  if(substr(allPhoto[i,1],1,1)=="B"){
    points(allDates[i],allPhoto[i,4],pch=20,col=cols[which(leafNames==allPhoto[i,1])])
  }else{
    points(allDates[i],allPhoto[i,4],pch=17,col=cols[which(leafNames==allPhoto[i,1])])
  }
}
dev.off()


