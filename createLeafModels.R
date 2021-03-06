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
species="oak"
calculateDIC <- TRUE
includeAllLeaves <- FALSE
c=3

n.cores <- 8
registerDoParallel(cores=n.cores)
# combinations <- subset(combinations, !as.logical(combinations$missingYear))
# combinations <- subset(combinations, !as.logical(combinations$excludePostSOS))

#modelVersion <- "Tair_D"

foreach(c=1:40) %dopar% {
#for(c in 1:nrow(combinations)){
  missingLeaf <- combinations$missingYear[c]
  excludePostSOS <- combinations$excludePostSOS[c]
  modelVersion <- combinations$Cov[c]
  addition <- combinations$Addition[c]
  
  if(addition){
    if(species=="beech"){
      model <- generalModel_B_CovPlusTair_D
    }else if(species=="oak"){
      model <- generalModel_O_CovPlusTair_D
    }
  }else if(modelVersion=="Tair_D"){
    if(species=="beech"){
      model <- generalModel_B_Tair_D
    }else if(species=="oak"){
      model <- generalModel_O_Tair_D
    }
  }else{
    if(species=="beech"){
      model <- generalModel_B_Cov
    }else if(species=="oak"){
      model <- generalModel_O_Cov
    }
  }
  if(includeAllLeaves){
    model <- generalModel_Tair_D
  }
  
  if(missingLeaf){
    if(addition){
      outputFileName <- paste0("modelFits/",species,"_",modelVersion,"_missingLeaf_CovPlusTair_varBurn.RData")
      partialFileName <- paste0("modelFits/",species,"_",modelVersion,"_missingLeaf_CovPlusTair_partial_varBurn.RData")
    }else{
      outputFileName <- paste0("modelFits/",species,"_",modelVersion,"_missingLeaf_varBurn.RData")
      partialFileName <- paste0("modelFits/",species,"_",modelVersion,"_missingLeaf_partial_varBurn.RData")
    }
  }else{
    if(addition){
      outputFileName <- paste0("modelFits/",species,"_",modelVersion,"_full_CovPlusTair_varBurn.RData")
      partialFileName <- paste0("modelFits/",species,"_",modelVersion,"_full_CovPlusTair_partial_varBurn.RData")
    }else{
      outputFileName <- paste0("modelFits/",species,"_",modelVersion,"_full_varBurn.RData")
      partialFileName <- paste0("modelFits/",species,"_",modelVersion,"_full_partial_varBurn.RData")
    }
  }
  if(excludePostSOS){
    if(addition){
      outputFileName <- paste0("modelFits/",species,"_",modelVersion,"_excludePostSOS_CovPlusTair_varBurn.RData")
      partialFileName <- paste0("modelFits/",species,"_",modelVersion,"_excludePostSOS_CovPlusTair_partial_varBurn.RData")
    }else{
      outputFileName <- paste0("modelFits/",species,"_",modelVersion,"_excludePostSOS_varBurn.RData")
      partialFileName <- paste0("modelFits/",species,"_",modelVersion,"_excludePostSOS_partial_varBurn.RData")
    }
  }
  #outputFileName <- gsub("_varBurn","_noTreeRE_varBurn",outputFileName)
  if(includeAllLeaves){
    outputFileName <- gsub("_varBurn","_allLeaves_varBurn",outputFileName)
  }

  dicFileName <- gsub("varBurn","dic",outputFileName)
  if(addition){
    inputFileName <- paste0("modelFits/harvard_",modelVersion,"_full_CovPlusTair_varBurn.RData")
  }else{
    inputFileName <- paste0("modelFits/harvard_",modelVersion,"_full_varBurn.RData")
  }
  if(!file.exists(dicFileName)){
  if(file.exists(outputFileName)){ ###*********CHANGE FOR NOT DIC
    print(outputFileName)
    
    #variableNames <- c("x","p.proc","b0","b3","b4","b0_leaf","b3_leaf","b4_leaf","b0_tree","b3_tree","b4_tree")
    variableNames <- c("x","p.proc","b0","b3","b4","b0_leaf","b3_leaf","b4_leaf")
    
    if(addition){
      #variableNames <- c(variableNames,"b2","b2_leaf","b2_tree")
      variableNames <- c(variableNames,"b2","b2_leaf")
    }
    if(species=="beech"){
      load(file=paste0('Data/finalData/',"allbeechLicorLeafData.RData"))
    }else if(species=="oak"){
      load(file=paste0('Data/finalData/',"alloakLicorLeafData.RData"))
    }
    if(includeAllLeaves){
      load(file=paste0('Data/finalData/',"all_noPhoto_",species,"LicorLeafData.RData"))
    }
    
    load(inputFileName)
    if(modelVersion=="TairOnly"){
      allData$Cov <- allData$Tair
    }else if(modelVersion=="GPP"){
      allData$Cov <- allData$GPP
    }else if(modelVersion=="NPP"){
      allData$Cov <- allData$NPP
    }else if(modelVersion=="NPP_feedback"){
      load(file=paste0('Data/finalData/',"all",species,"LicorLeafDataFeedback.RData"))
      
      allData$Cov  <- allData$NPP
      allData$Cov  <- abind(allData$Cov,allData$Cov[,100,,],along=2)
      if(addition){
        if(species=="beech"){
          model <- generalModel_B_feedback_CovPlusTair_D
        }else if(species=="oak"){
          model <- generalModel_O_feedback_CovPlusTair_D
        }
      }else{
        if(species=="beech"){
          model <- generalModel_B_feedback
        }else if(species=="oak"){
          model <- generalModel_O_feedback
        }
      }
    }else if(modelVersion=="GPP_feedback"){
      load(file=paste0('Data/finalData/',"all",species,"LicorLeafDataFeedback.RData"))
      
      allData$Cov  <- allData$GPP
      allData$Cov  <- abind(allData$Cov,allData$Cov[,100,,],along=2)
      if(addition){
        if(species=="beech"){
          model <- generalModel_B_feedback_CovPlusTair_D
        }else if(species=="oak"){
          model <- generalModel_O_feedback_CovPlusTair_D
        }
      }else{
        if(species=="beech"){
          model <- generalModel_B_feedback
        }else if(species=="oak"){
          model <- generalModel_O_feedback
        }
      }
    }
    
    out.mat <- as.data.frame(as.matrix(out.burn$params))
    mu <- -1*mean(out.mat$b0,na.rm=TRUE)
    vr <- var(out.mat$b0,na.rm = TRUE)
    allData$b0.a <- (mu**2-mu**3-mu*vr)/(vr)
    allData$b0.b <- (mu-2*mu**2+mu**3-vr+mu*vr)/(vr)
    allData$b0_prec <- 1/vr
    mu <- mean(out.mat$b3,na.rm=TRUE)
    vr <- var(out.mat$b3,na.rm = TRUE)
    allData$b3.a <- (mu**2-mu**3-mu*vr)/(vr)
    allData$b3.b <- (mu-2*mu**2+mu**3-vr+mu*vr)/(vr)
    mu <- -1*mean(out.mat$b4,na.rm=TRUE)
    vr <- var(out.mat$b4,na.rm = TRUE)
    allData$b4.a <- (mu**2-mu**3-mu*vr)/(vr)
    allData$b4.b <- (mu-2*mu**2+mu**3-vr+mu*vr)/(vr)
    allData$s1.PC <- 1.56
    allData$s2.PC <- 0.016
    allData$s1.proc <- 1.56
    allData$s2.proc <- 0.016
    allData$x1.a <- 10
    allData$x1.b <- 2
    if(addition){
      mu <- mean(out.mat$b2,na.rm=TRUE)
      vr <- var(out.mat$b2,na.rm = TRUE)
      allData$b2.a <- (mu**2-mu**3-mu*vr)/(vr)
      allData$b2.b <- (mu-2*mu**2+mu**3-vr+mu*vr)/(vr)
    }
    if(missingLeaf){
      if(includeAllLeaves){
        allData$CCI_means[,2] <- NA
      }else{
        allData$CCI_means[2,,2] <- NA
      }

    }
    if(excludePostSOS){
      for(t in 1:allData$treeN){ #Probably a cleaner way to do this, but I don't want to think about it now
        for(l in 1:allData$N){
          allData$CCI_means[l,(round(allData$tran[l,t],digits=0):allData$n),t] <- NA
        }
      }
    }
    
    j.model   <- jags.model(file = textConnection(model),
                            data = allData,
                            n.chains = nchain,
                            n.adapt = 1000)
    print("Done creating model")
    if(calculateDIC){
      load(outputFileName)
      var.sum <- summary(out.burn$params)
      DIC <- dic.samples(j.model,n.iter = var.sum$end)
      save(DIC,file=dicFileName)
    }else{
      out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,
                                  baseNum = 25000,iterSize = 10000,effSize = 5000,
                                  partialFile = partialFileName)
      
      out.burn <- thinMCMCOutput(out.burn)  ##Thin the data
      print(paste("saved:",outputFileName))
      save(out.burn,file = outputFileName)
    }
  }
  }
}