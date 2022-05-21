
load("Data/bbc1_dataFinal_includeJuly.RData")
#dataFinal$TairMuDay[,dataFinal$N]

plot(dates,dailyMet$airt.ac,pch=20,ylim=c(-10,30),col="#0868ac",cex=2,xlab="Date",ylab="Tair (C)",main="Average Day Time Temperature")
points(dates,dailyMet$airt.mid,pch=20,col=scales::alpha("#43a2ca"),cex=2)
points(dates,dailyMet$airt.low,pch=20,col=scales::alpha("#7bccc4"),cex=2)
points(dates,dailyMet$airt.us1,pch=20,col=scales::alpha("#bae4bc"),cex=2)
#points(dates,dataFinal$TairMuDay[,dataFinal$N],pch=20,cex=2,col="red")
legend('bottomleft',legend=metHeights,col=c("#0868ac","#43a2ca","#7bccc4","#bae4bc"),pch=rep(20,4),cex=1.5)

plot(dates,dailyMet$parac,pch=20,col="#0868ac",cex=2,xlab="Date",ylab="Tair (C)",main="Average Day Time PAR")
points(dates,dailyMet$parmid,pch=20,col=scales::alpha("#43a2ca"),cex=2)
points(dates,dailyMet$parlow,pch=20,col=scales::alpha("#7bccc4"),cex=2)
points(dates,dailyMet$parus1,pch=20,col=scales::alpha("#bae4bc"),cex=2)
#points(dates,dataFinal$TairMuDay[,dataFinal$N],pch=20,cex=2,col="red")
legend('bottomleft',legend=metHeights,col=c("#0868ac","#43a2ca","#7bccc4","#bae4bc"),pch=rep(20,4),cex=1.5)

plot(dates,dailyMet$rh.ac,pch=20,col="#0868ac",cex=2,xlab="Date",ylab="Tair (C)",main="Average Day Time Temperature")
points(dates,dailyMet$rh.mid,pch=20,col=scales::alpha("#43a2ca"),cex=2)
points(dates,dailyMet$rh.low,pch=20,col=scales::alpha("#7bccc4"),cex=2)
#points(dates,dataFinal$TairMuDay[,dataFinal$N],pch=20,cex=2,col="red")
legend('bottomleft',legend=metHeights[1:3],col=c("#0868ac","#43a2ca","#7bccc4"),pch=rep(20,3),cex=1.5)

###
# plot(dates,dailyMet$airt.ac-dailyMet$airt.mid,pch=20)
# abline(h=0,col="gray",lty=2)
# 
# plot(dates,dailyMet$airt.mid-dailyMet$airt.low,pch=20)
# abline(h=0,col="gray",lty=2)
# 
# plot(dates,dailyMet$airt.low-dailyMet$airt.us1,pch=20)
# abline(h=0,col="gray",lty=2)
# 
# plot(dates,dailyMet$airt.ac-dailyMet$airt.us1,pch=20)
# abline(h=0,col="gray",lty=2)
# abline(v=as.Date("2019-01-01"),col="red")
# 
# allDat <- merge(subMet,subMet2,'datetime')
# plot(allDat$datetime,allDat$airt.ac-allDat$airt.us1,pch=20,cex=0.25,xlab="Time",ylab="Difference",
#      main="Difference in Temperature Between Levels")
# #abline(h=0,col="gray",lty=2,lwd=2)
# 
# plot(allDat$datetime,allDat$airt.ac,pch=20,cex=0.25,xlab="Time")
# 
# 
# points(allDat$datetime,allDat$airt.ac-allDat$airt.low,pch=20,cex=0.25,col=scales::alpha("cyan",0.5))
# abline(h=0,col="red",lty=2,lwd=3)
# legend('bottomleft',legend=c("28m - 1m","28m - 9m"),col=c("black","cyan"),pch=rep(20,2))


# colfunc <- colorRampPalette(c("black", "white"))
# cols <- c(colfunc(25)[1:18],colfunc(30)[1:24])
# cols2 <- c(rep(cols[1],3),rep(cols[4],3),rep(cols[7],3),rep(cols[10],3),
#            rep(cols[13],3),rep(cols[16],3),rep(cols[19],3),rep(cols[22],3),
#            rep(cols[25],3),rep(cols[28],3),rep(cols[31],3),rep(cols[33],3),
#            rep(cols[36],3),rep(cols[39],3))
# #test <- output[order(output$Date),]
# jpeg("leafCCI_data.jpeg",width=7,height=6.4,units = "in",res=1000)
# par(mfrow=c(2,1))
# par(mai=c(0.8,1,0.5,0.1))
# subDat <- subset(output,Leaf==leafNames[19])
# leafHeight <- heightData$Height[heightData$LeafName==leafNames[19]]
# plot(subDat$Date,as.numeric(as.character(subDat$CCI_mean)),type = "l",
#      ylim=c(0,25),xlim=range(output$Date),main="Oak Leaves (light is higher)",ylab="CCI",xlab="",bty="n",col=cols2[19])
# for(i in 20:length(leafNames)){
#   subDat <- subset(output,Leaf==leafNames[i])
#   leafHeight <- heightData$Height[heightData$LeafName==leafNames[i]]
#   lines(subDat$Date,as.numeric(as.character(subDat$CCI_mean)),col=cols2[i])
# }
# 
# subDat <- subset(output,Leaf==leafNames[1])
# leafHeight <- heightData$Height[heightData$LeafName==leafNames[1]]
# plot(subDat$Date,as.numeric(as.character(subDat$CCI_mean)),type = "l",
#      ylim=c(0,25),xlim=range(output$Date),
#      xlab="Date", ylab="CCI",bty="n",main="Beech Leaves (light is higher)",col=cols2[1])
# for(i in 2:18){
#   subDat <- subset(output,Leaf==leafNames[i])
#   leafHeight <- heightData$Height[heightData$LeafName==leafNames[19]]
#   lines(subDat$Date,as.numeric(as.character(subDat$CCI_mean)),col=cols2[i])
# }
# 
# dev.off()

outputFileName <- paste0("modelFits/","hiearchicalTree_rescaledTemp_allBeech_TairchlorophyllCycling_varBurn.RData")
outputFileName <- paste0("modelFits/","hiearchicalTree_rescaledTemp_allBeech_NPPchlorophyllCycling_varBurn.RData")
load(outputFileName)
pred.out <- as.matrix(out.burn$predict)
ci=apply(as.matrix(out.burn$predict),MARGIN=2,FUN=quantile,c(0.025,0.5,0.975))

pred.cols = intersect(grep('\\[1,', colnames(pred.out)),grep('1\\]', colnames(pred.out)))
plot(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[1,,1],pch=20)
ecoforecastR::ciEnvelope(seq(1,allData$n),ci[1,pred.cols],ci[3,pred.cols],col="lightblue")
points(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[1,,1],pch=20)

pred.cols = intersect(grep('\\[2,', colnames(pred.out)),grep('1\\]', colnames(pred.out)))
plot(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[2,,1],pch=20)
ecoforecastR::ciEnvelope(seq(1,allData$n),ci[1,pred.cols],ci[3,pred.cols],col="lightblue")
points(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[2,,1],pch=20)

pred.cols = intersect(grep('\\[3,', colnames(pred.out)),grep('1\\]', colnames(pred.out)))
plot(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[3,,1],pch=20)
ecoforecastR::ciEnvelope(seq(1,allData$n),ci[1,pred.cols],ci[3,pred.cols],col="lightblue")
points(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[3,,1],pch=20)

pred.cols = intersect(grep('\\[1,', colnames(pred.out)),grep('2\\]', colnames(pred.out)))
plot(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[1,,2],pch=20)
ecoforecastR::ciEnvelope(seq(1,allData$n),ci[1,pred.cols],ci[3,pred.cols],col="lightblue")
points(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[1,,2],pch=20)

pred.cols = intersect(grep('\\[2,', colnames(pred.out)),grep('2\\]', colnames(pred.out))) #********
plot(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[2,,2],pch=20,main="NPP Covariate Model Fit",
     xlab="Time",ylab="Rescaled CCI")
ecoforecastR::ciEnvelope(seq(1,allData$n),ci[1,pred.cols],ci[3,pred.cols],col="lightblue")
points(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[2,,2],pch=20)

pred.cols = intersect(grep('\\[1,', colnames(pred.out)),grep('3\\]', colnames(pred.out)))
plot(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[1,,3],pch=20)
ecoforecastR::ciEnvelope(seq(1,allData$n),ci[1,pred.cols],ci[3,pred.cols],col="lightblue")
points(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[1,,3],pch=20)

pred.cols = intersect(grep('\\[2,', colnames(pred.out)),grep('3\\]', colnames(pred.out)))
plot(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[2,,3],pch=20)
ecoforecastR::ciEnvelope(seq(1,allData$n),ci[1,pred.cols],ci[3,pred.cols],col="lightblue")
points(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[2,,3],pch=20)

pred.cols = intersect(grep('\\[3,', colnames(pred.out)),grep('3\\]', colnames(pred.out)))
plot(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[3,,3],pch=20)
ecoforecastR::ciEnvelope(seq(1,allData$n),ci[1,pred.cols],ci[3,pred.cols],col="lightblue")
points(seq(1,length(allData$CCI_means[1,,1])),allData$CCI_means[3,,3],pch=20)

source('read_Licor.R')
source('fitA.R')
source('runLicorIter.R')
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
allPhoto <- matrix(ncol=4,nrow=0) #LeafName, dte, Vcmax, Jmax
allDates <- rep(Sys.Date(),51)
#allPhoto <- rbind(allPhoto,c("Test",as.Date(Sys.Date()),NA,NA))
j=1
# pdf(file="FieldWorkResponseCurvesConverged.pdf",height=5,width=12)
# par(mfrow=c(1,2))
foreach(l=1:length(leafNames)) %dopar% {
  #for(l in length(leafNames):1){
  print(leafNames[l])
  #files <- dir(path=dataDirectory,pattern=leafNames[l])
  files <- dir(path=dataDirectory,pattern=leafNames[l])
  # jmax <- numeric()
  # vcmax <- numeric()
  # dates <- rep(Sys.Date(),length(files))
  if(length(files)>0){
    for(f in 1:length(files)){
      fileName <- paste0(dataDirectory,files[f])
      print(fileName)
      
      dteChar <- strsplit(files[f],"-")[[1]][2]
      if(substr(dteChar,3,3)=="A"){
        dte <- as.Date(paste0("2021-08-",substr(dteChar,1,2)))
      }else if(substr(dteChar,3,3)=="S"){
        dte <- as.Date(paste0("2021-09-",substr(dteChar,1,2)))
      }
      print(dte)
      if(!file.exists(paste0('LicorFits/',leafNames[l],"_",dte,"_Cleaned_LicorResponseCurve_varBurn.RData"))){
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
          #dat <- dat[-(1:nrow(dat)),]
          next
        }else if(fileName =="Data/Licor_Measurements/HarvardForest_Licor_Wheeler-07SEP21-O4B"){
          dat <- dat[-13,]
        }else if(fileName =="Data/Licor_Measurements/HarvardForest_Licor_Wheeler-08SEP21-O2A"){
          dat <- dat[-13,]
        }else if(fileName =="Data/Licor_Measurements/HarvardForest_Licor_Wheeler-19SEP21-O1A"){ #Very high photosynthesis values...
          dat <- dat[-15,]
        }
        
        out <- fitA(flux.data = dat,hiearchical = FALSE)
        save(out,file=paste0('LicorFits/',leafNames[l],"_",dte,"_Cleaned_LicorResponseCurve_varBurn.RData"))
        print(paste("Saved:",paste0('LicorFits/',leafNames[l],"_",dte,"_Cleaned_LicorResponseCurve_varBurn.RData")))
      }
      # jmax <- c(jmax,mean(data.frame(as.matrix(out$params))$Jmax0))
      # vcmax <- c(vcmax,mean(data.frame(as.matrix(out$params))$vmax0))
      # dates[f] <- dte
      # allPhoto <- rbind(allPhoto,c(leafNames[l],dte,mean(data.frame(as.matrix(out$params))$Jmax0),mean(data.frame(as.matrix(out$params))$vmax0)))
      # allDates[j] <- dte
      # j <- j+1
      # plot_photo2(dat,out,byLeaf=FALSE,id=paste(leafNames[l],dte))
      #plot(dat$Ci,dat$Photo,pch=20,main=paste(leafNames[l],dte))
      #plot(dat$PARi,dat$Photo,pch=20,main=paste(leafNames[l],dte))
    }
    # plot(dates,vcmax,main=paste(leafNames[l],"Vcmax"),pch=20)
    # plot(dates,jmax,main=paste(leafNames[l],"Jmax"),pch=20)
  }
}
# dev.off()
# ##Note: leafNames changes to sorted by height
# jpeg("allPhotoParametersAgainstTime.jpeg",width=7,height=7,units = "in",res=1000)
# par(mfrow=c(2,1))
# par(mai=c(0.8,1,0.5,0.1))
# plot(allDates,as.numeric(allPhoto[,3]),pch=20,cex=0.2,xlab="Time",ylab="Jmax",ylim=c(20,280),main="Lighter is Higher Height (By Species)")
# for(l in 1:length(leafNames)){
#   subDat <- allPhoto[allPhoto[,1]==leafNames[l],]
#   subDts <- allDates[allPhoto[,1]==leafNames[l]]
#   lines(subDts[order(subDts)],as.numeric(subDat[,3])[order(subDts)],col=cols[l])
# }
# 
# for(i in 1:nrow(allPhoto)){
#   if(substr(allPhoto[i,1],1,1)=="B"){
#     points(allDates[i],allPhoto[i,3],pch=20,col=cols[which(leafNames==allPhoto[i,1])])
#   }else{
#     points(allDates[i],allPhoto[i,3],pch=17,col=cols[which(leafNames==allPhoto[i,1])])
#   }
# }
# legend('topright',c("Beech","Oak"),pch=c(20,17))
# 
# plot(allDates,as.numeric(allPhoto[,4]),pch=20,cex=0.2,xlab="Time",ylab="Vcmax",ylim=c(20,140),main="Lighter is Higher Height(By Species)")
# for(l in 1:length(leafNames)){
#   subDat <- allPhoto[allPhoto[,1]==leafNames[l],]
#   subDts <- allDates[allPhoto[,1]==leafNames[l]]
#   lines(subDts[order(subDts)],as.numeric(subDat[,4])[order(subDts)],col=cols[l])
# }
# for(i in 1:nrow(allPhoto)){
#   if(substr(allPhoto[i,1],1,1)=="B"){
#     points(allDates[i],allPhoto[i,4],pch=20,col=cols[which(leafNames==allPhoto[i,1])])
#   }else{
#     points(allDates[i],allPhoto[i,4],pch=17,col=cols[which(leafNames==allPhoto[i,1])])
#   }
# }
# dev.off()
