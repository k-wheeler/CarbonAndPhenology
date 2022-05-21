library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(doParallel)
library(scales)
source('generalModels.R')
source('sharedVariables.R')
source('generalFunctions.R')
#Fit Transitions ----

variables <- c("mS","mF","y[1]","k")
nchain=5
##' Create a Bayes Model for a deciduous broadleaf site
##'
##' @param yobs
##' @import rjags
##' @import runjags
##' @import PhenologyBayesModeling
##' @export
createChangepointModel_Fall <- function(yobs) {
  nchain = 5
  data <- list(yobs=yobs,x=seq(1,length(yobs)),n=length(yobs))
  data$obs.prec <- 100
  data$s1 <- 0.001
  data$s2 <- 0.00001
  data$c.alp <- 8
  data$c.bet <- 1.5
  DB_model <- "
  model{
  ##priors
  #prec ~ dgamma(s1,s2)
  y[1] ~ dbeta(c.alp,c.bet)
  mS ~ dunif(0,1)
  mF ~ dunif(0,1)
  k ~ dunif(0,183)
  
  for(i in 2:n){
  muS[i] <- y[(i-1)] - mS
  muF[i] <- y[(i-1)] - mF
  y[i] <- ifelse(x[i]>k,muF[i],muS[i])
  #y[i] ~ dnorm(mu[i],prec)
  yobs[i] ~ dnorm(y[i],obs.prec)
  }
  }
  "
  inits <- list()
  c("mS","mF","y[1]","k")
  for(i in 1:nchain){
    inits[[i]] <- list(mS = rnorm(1,0.0025,0.0001),
                       mF = rnorm(1,0.03,0.001),
                       k = rnorm(1,140,1))
  }
  
  j.model   <- jags.model(file = textConnection(DB_model),
                          data = data,
                          inits = inits,
                          n.chains = nchain,
                          n.adapt = 2000)
  return(j.model)
}
for(l in 1:length(leafNames)){
  lfName <- leafNames[l]
  print(lfName)
  load(paste0('Data/finalData/',lfName,"_finalData.RData"))
  
  outputFileName <- paste0("modelFits/",lfName,"_estimatedTransition_varBurn.RData")
  
  j.model <- createChangepointModel_Fall(yobs=finalData$CCI_means)

  var.burn <- runMCMC_Model(j.model = j.model,variableNames = variables, baseNum=10000,
                            iterSize = 5000,sampleCutoff = 5000)
  if(typeof(var.burn)!=typeof(FALSE)){
    out.mat <- as.matrix(var.burn)
    thinAmount <- round(nrow(out.mat)/5000,digits=0)
    var.burn <- window(var.burn,thin=thinAmount)
  }
  save(var.burn,file=outputFileName)
}

#Plot and Save Transitions with Data Objects ----
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(scales)
source('generalModels.R')
source('sharedVariables.R')
source('generalFunctions.R')

changepointModel <- function(y1,mS,mF,k,xseq){
  yseq <- y1
  for(x in xseq){
    if(x<=k){
      yseq <- c(yseq,yseq[length(yseq)]-mS)
    }else{
      yseq <- c(yseq,yseq[length(yseq)]-mF)
    }
  }
  return(yseq)
}

pdf(file="LeafTransitions.pdf",height=5,width=6)
for(l in 1:length(leafNames)){
  lfName <- leafNames[l]
  print(lfName)
  load(paste0(finalDataDirectory,lfName,"_finalData.RData"))
  
  outputFileName <- paste0("modelFits/",lfName,"_estimatedTransition_varBurn.RData")
  load(outputFileName)
  #plot(finalData$dates,finalData$CCI_means,pch=20,main=lfName,ylab="Rescaled CCI",xlab="Time")
  days <- seq(1,length(finalData$dates))
  var.mat <- data.frame(as.matrix(var.burn))
  ycred <- matrix(0,nrow=10000,ncol=length(days))
  
  rndNums <- sample(1:nrow(var.mat),10000,replace=T)
  mS <- var.mat$mS[rndNums]
  mF <- var.mat$mF[rndNums]
  y1 <- var.mat$y.1[rndNums]
  k <- var.mat$k[rndNums]
  for(g in 1:10000){
    ycred[g,] <- changepointModel(y1=y1[g],mS=mS[g],mF=mF[g],k=k[g],xseq=days[2:length(days)])
  }
  ci <- apply(ycred,2,quantile,c(0.025,0.5,0.975))
  tran.ci <- quantile(k,c(0.025,0.5,0.975))
  plot(days,finalData$CCI_means,pch=20,main=lfName,ylab="Rescaled CCI",xlab="Time")
  ecoforecastR::ciEnvelope(days,ci[1,],ci[3,],col=alpha("lightblue",0.5))
  points(finalData$dates,finalData$CCI_means,pch=20)
  abline(v=tran.ci[2],col="blue")
  abline(v=tran.ci[1],col="blue",lty=2)
  abline(v=tran.ci[3],col="blue",lty=2)
  finalData$tran <- c(tran.ci,mean(k))
  save(file=paste0(finalDataDirectory,lfName,"_finalData_withTran.RData"),finalData)
  
}
dev.off()
