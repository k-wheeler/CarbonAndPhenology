##Create Calibration Fits for mean air temperature no trigger model
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(doParallel)
dataDirectory <- "Data/finalData/"
leafNames <- c("B1A","B1B","B1C","B1D","B1E","B1F",
               "B2A","B2B","B2C","B2D","B2E","B2F",
               "B3A","B3B","B3C","B3D","B3E","B3F",
               "O1A","O1B","O1C","O1D","O1E","O1F",
               "O2A","O2B","O2C","O2D","O2E","O2F",
               "O3A","O3B","O3C","O3D","O3E","O3F",
               "O4A","O4B","O4C","O4D","O4E","O4F")

variableNames <- c("x","p.proc","b0","b3","b4")
nchain=5
generalModel = "
model {
    ### Data Models for complete years
    for(i in 1:n){
    CCI_means[i] ~ dnorm(x[i],CCI_precs[i])
    }
    
    #### Process Model
    for(i in 2:n){
    
    xmu[i] <- max(min((x[(i-1)] + -1*b4 * x[(i-1)]) + max(0,(-1*b0 + b3 * Tair[i] * D[i])),x[1]),0)
    x[i] ~ dnorm(xmu[i],p.proc)
    }
    
    #### Priors
    x[1] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    p.proc ~ dgamma(s1.proc,s2.proc)
    #p.PC ~ dgamma(s1.PC,s2.PC)
    b0 ~ dbeta(b0.a,b0.b) 
    b3 ~ dbeta(b3.a,b3.b) 
    b4 ~ dbeta(b4.a,b4.b) 
    
  }
    "

for(l in 1:length(leafNames)){
  lfName <- leafNames[l]
  print(lfName)
  load(paste0('Data/finalData/',lfName,"_finalData.RData"))
  
  outputFileName <- paste0("modelFits/",lfName,"_OGchlorophyllCycling_varBurn.RData")
  partialFileName <- paste0("modelFits/",lfName,"_OGchlorophyllCycling_partial_varBurn.RData")
  
  #load('bbc1_meanTemp_summer183_expBreak_slope_forecast_b3_calibration_varBurn.RData')
  load("modelFits/bbc1_fullAutumn_dayTime_varBurn.RData")
  out.mat <- as.data.frame(as.matrix(out.burn$params))
  mu <- -1*mean(out.mat$b0,na.rm=TRUE)
  vr <- var(out.mat$b0,na.rm = TRUE)
  finalData$b0.a <- (mu**2-mu**3-mu*vr)/(vr)
  finalData$b0.b <- (mu-2*mu**2+mu**3-vr+mu*vr)/(vr)
  mu <- mean(out.mat$b3,na.rm=TRUE)
  vr <- var(out.mat$b3,na.rm = TRUE)
  finalData$b3.a <- (mu**2-mu**3-mu*vr)/(vr)
  finalData$b3.b <- (mu-2*mu**2+mu**3-vr+mu*vr)/(vr)
  mu <- -1*mean(out.mat$b4,na.rm=TRUE)
  vr <- var(out.mat$b4,na.rm = TRUE)
  finalData$b4.a <- (mu**2-mu**3-mu*vr)/(vr)
  finalData$b4.b <- (mu-2*mu**2+mu**3-vr+mu*vr)/(vr)
  
  finalData$s1.PC <- 1.56
  finalData$s2.PC <- 0.016
  finalData$s1.proc <- 1.56
  finalData$s2.proc <- 0.016
  # finalData$b0_lower <- -0.25
  # finalData$b0_upper <- 0
  # finalData$b3_lower <- 0
  # finalData$b3_upper <- 0.25
  # finalData$b4_lower <- -0.25
  # finalData$b4_upper <- 0
  #Initial Conditions (uninformed and probably need to pick better priors)
  finalData$x1.a <- 10
  finalData$x1.b <- 2

  j.model   <- jags.model(file = textConnection(generalModel),
                          data = finalData,
                          n.chains = nchain,
                          n.adapt = 1000)
  out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,
                              baseNum = 50000,iterSize = 10000,effSize = 4000,
                              partialFile = partialFileName)
  ##Thin the data:
  out.mat <- as.matrix(out.burn$params)
  thinAmount <- round(nrow(out.mat)/5000,digits=0)
  out.burn2 <- list()
  out.burn2$params <- window(out.burn$params,thin=thinAmount)
  out.burn2$predict <- window(out.burn$predict,thin=thinAmount)
  out.burn <- out.burn2
  save(out.burn,file = outputFileName)
  summary(as.matrix(out.burn$params))

  ci=apply(as.matrix(out.burn$predict),MARGIN=2,FUN=quantile,c(0.025,0.5,0.975))
  plot(finalData$dates,finalData$CCI_means,pch=20)
  ecoforecastR::ciEnvelope(finalData$dates,ci[1,],ci[3,],col="lightblue")
  points(finalData$dates,finalData$CCI_means,pch=20)

}
load(paste0('Data/finalData/','B1A',"_finalData.RData"))
Temps <- finalData$Tair
forecastStep <- function(IC,b0,b1,b2,b3,b4,Q=0,n,NT,Tair,D){
  x <- matrix(NA,n,NT)
  if(length(IC)==1){
    Xprev <- rep(IC,n)
  }else{
    Xprev <- IC
  }
  
  for(t in 1:NT){
    bd <- Xprev + b4 * Xprev #b0 + (b1 * Xprev) + (b2 * Xprev **2) 
    syn <- b0 + b1 * Tair[t] + b2 * D[t] + b3 * Tair[t] * D[t] 
    if(length(syn)==1){
      syn <- rep(syn,n)
    }
    
    xNew <- numeric()
    mu <- rep(NA,n)
    for(i in 1:n){
      mu[i] <- bd[i] + max(syn[i],0)
      mu[i] <- max(0,min(mu[i],IC[min(i,length(IC))]))
      
      if(length(Q)>1){
        xNew <- c(xNew,rnorm(1,mu[i],Q[i]))
      }else{
        xNew <- c(xNew,rnorm(1,mu[i],Q))
      }
    }
    x[,t] <- xNew
    Xprev <- x[,t]
  }
  return(x)
}
load('bbc1_meanTemp_summer183_expBreak_slope_forecast_b3_calibration_varBurn.RData')
old.mat <- data.frame(as.matrix(out.burn$param))
load("modelFits/bbc1_fullAutumn_dayTime_varBurn.RData")
load("Data/bbc1_dataFinal_includeJuly.RData")
Nmc=1000
out.mat <- data.frame(as.matrix(out.burn$param))
prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)

b0 <- out.mat$b0[prow]
b4 <- out.mat$b4[prow]
b3 <- out.mat$b3[prow]
b1 <- rep(0,Nmc)
b2 <- rep(0,Nmc)
dates <- seq(as.Date('2021-08-01'),as.Date('2021-11-30'),'day')
NT=length(dates)
pdfName <- paste0("trialLeafChlorophyllCyclingModelFits_bbcReFit2.pdf")
pdf(file=pdfName,height=5,width=7)
for(l in 1:length(leafNames)){
  lfName <- leafNames[l]
  print(lfName)
  load(paste0('Data/finalData/',lfName,"_finalData.RData"))
  # ysPred <- forecastStep(IC=rbeta(Nmc,10,0.3),b0=b0,b1=b1,b2=b2,b3=b3,b4=b4,
  #                        n=Nmc,NT=NT,Tair=finalData$Tair,D=finalData$D)
  ysPred <- forecastStep(IC=rbeta(Nmc,10,0.3),b0=b0,b1=b1,b2=b2,b3=b3,b4=b4,
                         n=Nmc,NT=NT,Tair=Temps,D=finalData$D)
  ci = apply(ysPred,2,quantile,c(0.025,0.5,0.975))
  plot(finalData$dates,finalData$CCI_means,pch=20,main=lfName,ylim=c(0,1))
  ecoforecastR::ciEnvelope(dates,ci[1,],ci[3,],col="lightblue")
  points(finalData$dates,finalData$CCI_means,pch=20)
}
dev.off()