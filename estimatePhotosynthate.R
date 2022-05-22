#Load in Libraries and Other Files ----
library(rjags)
library(runjags)
library(stats)
library(PEcAn.photosynthesis)
library(RColorBrewer)
source('generalModels.R')
source('sharedVariables.R')
source('generalFunctions.R')
options(stringsAsFactors = FALSE)

treeSpecies <- "oak"
licorReponseFit <- paste0('LicorFits/',treeSpecies,"_QAQC_August_allRE_LicorResponseCurve_varBurn.RData")
load(licorReponseFit)
out.mat <- as.data.frame(as.matrix(out$params))

#Load All Licor Data ----
Allfiles <- paste0(LicorDataDirectory,dir(path=LicorDataDirectory,pattern="-"))

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

#Farquhar Model Defined ----
#From: https://github.com/PecanProject/pecan/blob/develop/modules/photosynthesis/code/FBB_functions.R
farquhar = function(ci,Fparams,q,Temp,typ){
  vmax0 <- Fparams[1]
  Jmax0 <- Fparams[2]
  r0 <- Fparams[3]
  cp0 <- Fparams[4]
  alpha <- Fparams[5]
  
  R <- 8.3144621/1000 ## gas constant
  r.c <- 18.72
  r.H <- 46.39 
  Vc.c <- 26.35
  Vc.H <- 65.33
  Vo.c <- 22.98
  Vo.H <- 60.11
  cp.c <- 19.02
  cp.H <- 37.83
  cp.ref <- 42.75
  Kc.c <- 38.05
  Kc.H <- 79.43
  Kc.ref <- 404.9
  Ko.c <- 20.30
  Ko.H <- 36.38
  Ko.ref <- 278.4
  Kc <- 46
  Ko <- 22000#33000?
  po <- 21000
  
  ## Constants: June et al 2004, Funct Plant Bio
  Omega <- 18
  To <- 37+273.15    ## Representative value, would benefit from spp calibration!
  
  #Respiration
  r  <- r0 * exp(r.c - r.H/R/Temp) #Temp in Kelvin
  cp <- cp0 * exp(cp.c - cp.H/R/Temp)/cp.ref
  Kc.T <- Kc * exp(Kc.c - Kc.H/R/Temp)/Kc.ref
  Ko.T <- Ko * exp(Ko.c - Ko.H/R/Temp)/Ko.ref
  Jmax <- Jmax0 * exp(-(Temp-To)*(Temp-To)/(Omega*Omega))
  vmax <- vmax0 * exp(Vc.c - Vc.H/R/Temp)
  #q is PARi
  aj <- (alpha*q/(sqrt(1+(alpha*alpha*q*q)/(Jmax*Jmax))))*(ci-cp)/(4*ci+8*cp)    ## electron transport limited
  ac <- vmax*(ci-cp)/(ci+Kc.T*(1+po/Ko.T))       ## maximum rubisco limited without covariates
  if(typ=="NPP"){
    min(aj,ac) - r
  }else if(typ=="GPP"){
    min(aj,ac)
  }else if(typ=="R"){
    r
  }else if(typ=="aj"){
    aj
  }else if(typ=="ac"){
    ac
  }

}

#Medlyn Model Defined ----
medlyn = function(input,Mparams,Fparams,obs,typ){
  Temp <- obs[4] + 273.15 #Convert to Kelvin
  Ci <- obs[1] - input[1]/input[2]
  e1 <- farquhar(Ci,Fparams,obs[3],Temp=Temp,typ=typ) - input[1]
  e2 <- Mparams[1] + (1 + Mparams[2]/sqrt(obs[2])*input[1]/obs[1]-input[2])*100
  return(e1^2 + e2^2)
}

#Function to Solve Series-of-Equations ----
solve.model = function(vmax0,Jmax0,r0,cp0,alpha, m, g0,finalData,allDat,typ="NPP"){
  output = list()
  print(c(vmax0,Jmax0,r0,cp0,alpha))
  for(i in 1:finalData$n){					# loop over data points
    ic <- c(mean(allDat$Photo), mean(allDat$Cond)) 		# take initial conditions from actual data
    out <- optim(ic,			# solve simultaneously for An.pred and gs.pred
                 medlyn,
                 Mparams = c(g0,m),	        # Ballberry params 
                 Fparams = c(vmax0,Jmax0,r0,cp0,alpha),	# Farquhar params
                 obs = c(finalData$co2[i], finalData$vpd[i]+0.001, finalData$par[i],finalData$Tair[i]),# data
                 typ = typ)  			
    output$An.pred[i] = out$par[1]
    output$gs.pred[i] = out$par[2]/100 #Scale back to
  }	
  return(output$An.pred)
}

# ci <- 400 ###Change
#   
# out.mat.names <- names(out.mat)
# r0 <- mean(out.mat$r0)
# cp0 <- mean(out.mat$cp0)
# pdfName <- "estimatedPhotosynthate_oak.pdf"
# pdf(file=pdfName,height=8,width=10)
# par(mfcol=c(1,1))
# for(l in 1:length(licorLeafNames)){
#   lfName <- licorLeafNames[l]
#   print(lfName)
#   load(paste0('Data/finalData/',lfName,"_finalData.RData"))
#   cl <- paste0("Vleaf[",licorIndices[l],"]")
#   vmax0 <- mean(out.mat[,which(out.mat.names==cl)] + (out.mat$vmax0))
#   
#   cl <- paste0("Aleaf[",licorIndices[l],"]")
#   alpha <- mean(out.mat[,which(out.mat.names==cl)] + (out.mat$alpha0))
#   if(substr(lfName,1,1)=="B"){
#     m=3.737588
#     g0=0.05215698
#   }else if(substr(lfName,1,1)=="O"){
#     m=4.85734
#     g0=0.1852357
#   }
#   output <- solve.model(vmax0,Jmax0,r0,cp0,alpha, m, g0,finalData)
#   
#   #plot(finalData$dates,finalData$CCI_means[1:finalData$n],pch=20,main=lfName)
#   #plot(finalData$dates,finalData$Tair[1:finalData$n],type="l",main=lfName)
#   if(l==1){
#     plot(finalData$dates,output$An.pred,type="l",ylab="An",xlab="Time")
#   }else{
#     lines(finalData$dates,output$An.pred,col=l)
#   }
# 
#   #plot(finalData$dates,output$gs.pred,type="l",ylab="gs",xlab="Time",main=lfName)
#   print(quantile(finalData$Tair,c(0.025,0.5,0.975)))
# }
# dev.off()

#Calculate Photosynthate for PhenoCam (mean Oak)----
treeSpecies <- "oak"
licorReponseFit <- paste0('LicorFits/',treeSpecies,"_QAQC_August_allRE_LicorResponseCurve_varBurn.RData")
load(licorReponseFit)
#load(paste0(finalDataDirectory,"bbc1_dataFinal.RData"))
load(paste0(finalDataDirectory,"harvard_dataFinal.RData"))

out.mat <- as.data.frame(as.matrix(out$params))
vmax0 <- mean(out.mat$vmax0)
alpha <- mean(out.mat$alpha0)
r0 <- mean(out.mat$r0)
cp0 <- mean(out.mat$cp0)
Jmax0 <- mean(out.mat$Jmax0)
m=4.85734
g0=0.1852357
dataFinal$NPP <- matrix(nrow=dataFinal$n,ncol=dataFinal$N)
dataFinal$GPP <- matrix(nrow=dataFinal$n,ncol=dataFinal$N)
dataFinal$R <- matrix(nrow=dataFinal$n,ncol=dataFinal$N)
dataFinal$aj <- matrix(nrow=dataFinal$n,ncol=dataFinal$N)
dataFinal$ac <- matrix(nrow=dataFinal$n,ncol=dataFinal$N)
for(yr in 1:dataFinal$N){
  finalData <- list(co2=dataFinal$co2[,yr], 
                    vpd=dataFinal$vpd[,yr], 
                    par=dataFinal$par[,yr],
                    Tair=dataFinal$TowerTair[,yr],
                    n=dataFinal$n)
  dataFinal$NPP[,yr] <- solve.model(vmax0,Jmax0,r0,cp0,alpha, m, g0,finalData,allDat,typ="NPP")
  dataFinal$GPP[,yr] <- solve.model(vmax0,Jmax0,r0,cp0,alpha, m, g0,finalData,allDat,typ="GPP")
  dataFinal$R[,yr] <- solve.model(vmax0,Jmax0,r0,cp0,alpha, m, g0,finalData,allDat,typ="R")
  dataFinal$aj[,yr] <- solve.model(vmax0,Jmax0,r0,cp0,alpha, m, g0,finalData,allDat,typ="aj")
  dataFinal$ac[,yr] <- solve.model(vmax0,Jmax0,r0,cp0,alpha, m, g0,finalData,allDat,typ="ac")
}
#plot(dataFinal$aj[,1],pch=20)
#save(file=paste0(finalDataDirectory,"bbc1_dataFinal_withPhoto.RData"),dataFinal)
save(file=paste0(finalDataDirectory,"harvard_dataFinal_withPhoto.RData"),dataFinal)

#Calculate Photosynthate for Individual Leaves and Plot ----
# n <- length(licorLeafNames)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# cls = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# pdfName <- "estimatedPhotosynthate_all.pdf"
# pdf(file=pdfName,height=8,width=10)
# par(mfcol=c(1,1))
# for(typ in c("NPP","GPP","R")){
if(treeSpecies=="beech"){
  specificLicorLeafNames <- beechLicorAugLeafNames
}else if(treeSpecies=="oak"){
  specificLicorLeafNames <- oakLicorAugLeafNames
}

for(l in 1:length(specificLicorLeafNames)){
  lfName <- specificLicorLeafNames[l]
  print(lfName)
  load(paste0(finalDataDirectory,lfName,"_finalData_withTran.RData"))
  
  vmax0 <- mean(out.mat[,intersect(grep('Vleaf', colnames(out.mat)),grep(paste0(l,'\\]'), colnames(out.mat)))]+out.mat$vmax0)
  alpha <- mean(out.mat[,intersect(grep('Aleaf', colnames(out.mat)),grep(paste0(l,'\\]'), colnames(out.mat)))]+out.mat$alpha0)
  r0 <- mean(out.mat[,intersect(grep('rleaf', colnames(out.mat)),grep(paste0(l,'\\]'), colnames(out.mat)))]+out.mat$r0)
  cp0 <- mean(out.mat[,intersect(grep('cleaf', colnames(out.mat)),grep(paste0(l,'\\]'), colnames(out.mat)))]+out.mat$cp0)
  Jmax0 <- mean(out.mat[,intersect(grep('Jleaf', colnames(out.mat)),grep(paste0(l,'\\]'), colnames(out.mat)))]+out.mat$Jmax0)
  if(substr(lfName,1,1)=="B"){
    m=3.737588
    g0=0.05215698
  }else if(substr(lfName,1,1)=="O"){
    m=4.85734
    g0=0.1852357
  }
  finalData$NPP <- solve.model(vmax0,Jmax0,r0,cp0,alpha, m, g0,finalData,allDat,typ="NPP")
  finalData$GPP <- solve.model(vmax0,Jmax0,r0,cp0,alpha, m, g0,finalData,allDat,typ="GPP")
  finalData$R <- solve.model(vmax0,Jmax0,r0,cp0,alpha, m, g0,finalData,allDat,typ="R")
  finalData$aj <- solve.model(vmax0,Jmax0,r0,cp0,alpha, m, g0,finalData,allDat,typ="aj")
  finalData$ac <- solve.model(vmax0,Jmax0,r0,cp0,alpha, m, g0,finalData,allDat,typ="ac")
  # 
  # #plot(finalData$dates,finalData$CCI_means[1:finalData$n],pch=20,main=lfName)
  # #plot(finalData$dates,finalData$Tair[1:finalData$n],type="l",main=lfName)
  # if(typ=="NPP"){
  #   ylim=c(-1,25)
  #   finalData$NPP <- output$An.pred
  # }else if(typ=="GPP"){
  #   ylim=c(-1,30)
  #   finalData$GPP <- output$An.pred
  # }else if(typ=="R"){
  #   ylim=c(0,3)
  #   finalData$R <- output$An.pred
  # }
  # if(l==1){
  #   plot(finalData$dates,output$An.pred,type="l",ylab="An",xlab="Time",col=cls[l],main=typ,ylim=ylim)
  # }else{
  #   lines(finalData$dates,output$An.pred,col=cls[l])
  # }
  save(file=paste0('Data/finalData/',lfName,"_finalData_withTran_withPhoto.RData"),finalData)
  
  #plot(finalData$dates,output$gs.pred,type="l",ylab="gs",xlab="Time",main=lfName)
  #print(quantile(finalData$Tair,c(0.025,0.5,0.975)))
}
  #legend("topright",col=cls,lwd=rep(1,length(cls)),licorLeafNames)

#dev.off()


licorLeafNames <- c("B1A","B1B","B1E","B2A","B2B","B3A","B3B","B3E","O1A","O1B","O1E","O2A","O2B","O2E","O3A","O3B","O4A","O4B")
library(RColorBrewer)
n <- length(licorLeafNames)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
cls = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
transDates <- numeric()
for(l in 1:length(licorLeafNames)){
  lfName <- licorLeafNames[l]
  load(paste0('Data/finalData/',lfName,"_finalData_withTran.RData"))
  transDates <- c(transDates,finalData$tran[4])

}

colfunc <- colorRampPalette(c("black", "white"))
#colfunc <- colorRampPalette(c("red","yellow"))
#cols <- brewer.pal(61,"YlOrRd")
max(transDates)- min(transDates)
tranSeq <- seq(59,119)
cols <- colfunc(65)[1:61]


pdfName <- "estimatedPhotosynthate_all4.pdf"
pdf(file=pdfName,height=8,width=10)
par(mfcol=c(1,1))
for(l in 1:length(licorLeafNames)){
  lfName <- licorLeafNames[l]
  print(lfName)
  load(paste0('Data/finalData/',lfName,"_finalData_withTran.RData"))
  ylim=c(-1,25)
  cl <- cols[which(tranSeq==round(finalData$tran[4]))]
  if(l==1){
    plot(finalData$dates,finalData$NPP,type="l",ylab="NPP",xlab="Time",col=cl,main="NPP",ylim=ylim)
  }else{
    if(substr(lfName,1,1)=="B"){
      lines(finalData$dates,finalData$NPP,col=cl)
    }else{
      lines(finalData$dates,finalData$NPP,col=cl,lty=2)
    }

  }
  
}
legend("topright",lty=c(1,2),c("Beech","Oak"))
ylim=c(-1,30)
for(l in 1:length(licorLeafNames)){
  lfName <- licorLeafNames[l]
  print(lfName)
  load(paste0('Data/finalData/',lfName,"_finalData_withTran.RData"))
  cl <- cols[which(tranSeq==round(finalData$tran[4]))]
  if(l==1){
    plot(finalData$dates,finalData$GPP,type="l",ylab="GPP",xlab="Time",col=cl,main="GPP",ylim=ylim)
  }else{
    if(substr(lfName,1,1)=="B"){
      lines(finalData$dates,finalData$GPP,col=cl)
    }else{
      lines(finalData$dates,finalData$GPP,col=cl,lty=2)
    }
  }
}
legend("topright",lty=c(1,2),c("Beech","Oak"))
ylim=c(0,3)
for(l in 1:length(licorLeafNames)){
  lfName <- licorLeafNames[l]
  print(lfName)
  load(paste0('Data/finalData/',lfName,"_finalData_withTran.RData"))
  cl <- cols[which(tranSeq==round(finalData$tran[4]))]
  if(l==1){
    plot(finalData$dates,finalData$R,type="l",ylab="R",xlab="Time",col=cl,main="R",ylim=ylim)
  }else{
    if(substr(lfName,1,1)=="B"){
      lines(finalData$dates,finalData$R,col=cl)
    }else{
      lines(finalData$dates,finalData$R,col=cl,lty=2)
    }
  }
  
}
legend("topright",lty=c(1,2),c("Beech","Oak"))
  #legend("topright",col=cls,lwd=rep(1,length(cls)),licorLeafNames)

dev.off()

for(l in 1:length(licorLeafNames)){
  lfName <- licorLeafNames[l]
  print(lfName)
  load(paste0('Data/finalData/',lfName,"_finalData_withTran.RData"))
  cl <- cols[which(tranSeq==round(finalData$tran[4]))]
  if(l==1){
    plot(finalData$dates,finalData$co2[1:122],type="l",ylab="RH",xlab="Time",col=cl,main="CO2")#,ylim=c(0,110))
  }else{
    if(substr(lfName,1,1)=="B"){
      lines(finalData$dates,finalData$R[1:122],col=cl)
    }else{
      lines(finalData$dates,finalData$RH[1:122],col=cl,lty=2)
    }
  }
  
}







