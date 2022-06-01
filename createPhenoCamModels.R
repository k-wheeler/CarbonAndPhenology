#Create Data Object:
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(suncalc)
library(doParallel)
library(abind)
source('generalModels.R')
source('sharedVariables.R')
source('generalFunctions.R')
source('metFunctions.R')
options(stringsAsFactors=FALSE)

calculateDIC <- TRUE
c <- 37

n.cores <- 4
registerDoParallel(cores=n.cores)

#modelVersion <- "Tair_D"
#c(31,33,34,36,37,38,39,40)
#foreach(c=c(31,34)) %dopar% {
foreach(c=c(31,33,34,36,37,38,39,40)) %dopar% {
#for(c in 1:nrow(combinations)){
  
  missingYear <- combinations$missingYear[c]
  excludePostSOS <- combinations$excludePostSOS[c]
  #modelVersion <- "TairOnly"
  modelVersion <- combinations$Cov[c]
  addition <- combinations$Addition[c]
  interaction <- combinations$interaction[c]
  
  #dates <- seq(as.Date('2016-01-01'),as.Date('2020-12-31'),'day')
  dates <- seq(as.Date('2010-01-01'),as.Date('2020-12-31'),'day')
  dates <- dates[format(dates,"%j")%in%seq(213,335)]
  
  #load(paste0(finalDataDirectory,"bbc1_dataFinal_withPhoto.RData"))
  load(paste0(finalDataDirectory,"harvard_dataFinal_withPhoto.RData"))
  dataFinal$n <- 122
  variableNames <- c("x","p.proc","b0","b3","b4","p.PC")
  if(addition){
    variableNames <- c(variableNames,"b2")
  }
  print(variableNames)
  if(interaction){
    model <- generalModel_PC_Interaction 
    variableNames <- c(variableNames,"b1")
  }else if(addition){
    model <- generalModel_PC_CovPlusTair_D
  }else if(modelVersion=="Tair_D"){
    model <- generalModel_PC_Tair_D
  }else{
    model <- generalModel_PC_Cov
  }
  if(modelVersion=="TairOnly"){
    dataFinal$Cov <- dataFinal$TowerTair
  }else if(modelVersion=="GPP"){
    dataFinal$Cov <- dataFinal$GPP
  }else if(modelVersion=="NPP"){
    dataFinal$Cov <- dataFinal$NPP
  }else if(modelVersion=="GPP_aging"){
    dataFinal$Cov <- dataFinal$GPP_aging
  }else if(modelVersion=="NPP_aging"){
    dataFinal$Cov <- dataFinal$NPP_aging
  }else if(modelVersion=="NPP_feedback"){
    load(paste0(finalDataDirectory,"harvard_dataFinal_withPhotoFeedback.RData"))
    dataFinal$Cov <- dataFinal$NPP
    dataFinal$Cov <- abind(dataFinal$Cov,dataFinal$Cov[,,100],along=3)
    if(addition){
      model <- generalModel_PC_feedback_CovPlusTair_D
    }else{
      model <- generalModel_PC_feedback
    }

  }else if(modelVersion=="GPP_feedback"){
    load(paste0(finalDataDirectory,"harvard_dataFinal_withPhotoFeedback.RData"))
    dataFinal$Cov <- dataFinal$GPP
    dataFinal$Cov <- abind(dataFinal$Cov,dataFinal$Cov[,,100],along=3)
    if(addition){
      model <- generalModel_PC_feedback_CovPlusTair_D
    }else{
      model <- generalModel_PC_feedback
    }
  }

  if(missingYear){
    dataFinal$p[,2] <- NA 
    if(addition){
      if(interaction){
        outputFileName <- paste0("modelFits/harvard_",modelVersion,"_missingYear_interaction_varBurn.RData")
        partialFileName <- paste0("modelFits/harvard_",modelVersion,"_missingYear_interaction_partial_varBurn.RData")
      }else{
        outputFileName <- paste0("modelFits/harvard_",modelVersion,"_missingYear_CovPlusTair_varBurn.RData")
        partialFileName <- paste0("modelFits/harvard_",modelVersion,"_missingYear_CovPlusTair_partial_varBurn.RData")
      }
    }else{
      outputFileName <- paste0("modelFits/harvard_",modelVersion,"_missingYear_varBurn.RData")
      partialFileName <- paste0("modelFits/harvard_",modelVersion,"_missingYear_partial_varBurn.RData")
    }
  }else{
    if(addition){
      if(interaction){
        outputFileName <- paste0("modelFits/harvard_",modelVersion,"_full_interaction_varBurn.RData")
        partialFileName <- paste0("modelFits/harvard_",modelVersion,"_full_interaction_partial_varBurn.RData")
      }else{
        outputFileName <- paste0("modelFits/harvard_",modelVersion,"_full_CovPlusTair_varBurn.RData")
        partialFileName <- paste0("modelFits/harvard_",modelVersion,"_full_CovPlusTair_partial_varBurn.RData")
      }
    }else{
      outputFileName <- paste0("modelFits/harvard_",modelVersion,"_full_varBurn.RData")
      partialFileName <- paste0("modelFits/harvard_",modelVersion,"_full_partial_varBurn.RData")
    }
  }
  if(excludePostSOS){
    dataFinal$p[(dataFinal$sofMean-31):dataFinal$n,] <- NA
    if(addition){
      if(interaction){
        outputFileName <- paste0("modelFits/harvard_",modelVersion,"_excludePostSOS_interaction_varBurn.RData")
        partialFileName <- paste0("modelFits/harvard_",modelVersion,"_excludePostSOS_interaction_partial_varBurn.RData")
      }else{
        outputFileName <- paste0("modelFits/harvard_",modelVersion,"_excludePostSOS_interaction_varBurn.RData")
        partialFileName <- paste0("modelFits/harvard_",modelVersion,"_excludePostSOS_interaction_partial_varBurn.RData")
      }
    }else{
      outputFileName <- paste0("modelFits/harvard_",modelVersion,"_excludePostSOS_varBurn.RData")
      partialFileName <- paste0("modelFits/harvard_",modelVersion,"_excludePostSOS_partial_varBurn.RData")
    }
  }
  dicFileName <- gsub("varBurn","dic",outputFileName)
  
  print(outputFileName)
  print(dicFileName)
  if(file.exists(outputFileName) & !file.exists(dicFileName)){#******Change for non-DIC
  dataFinal$s1.PC <- 1.56
  dataFinal$s2.PC <- 0.016
  dataFinal$s1.proc <- 1.56
  dataFinal$s2.proc <- 0.016
  dataFinal$b0_lower <- -0.5
  dataFinal$b0_upper <- 0
  dataFinal$b3_lower <- 0
  dataFinal$b3_upper <- 0.5
  dataFinal$b4_lower <- -0.5
  dataFinal$b4_upper <- 0
  dataFinal$b1_upper <- 0.5
  dataFinal$b1_lower <- -0.5
  print(dim(dataFinal$Cov))
  j.model <- try(jags.model(file = textConnection(model),
                            data = dataFinal,
                            n.chains = nchain,
                            n.adapt = 1000))#Load Model Output 
  print("Done with Creating Model")
  if(calculateDIC){
    
    if(file.exists(outputFileName)){
      if(!file.exists(dicFileName)){
      print("Calculating DIC")
      load(outputFileName)
      var.sum <- summary(out.burn$params)
      DIC <- dic.samples(j.model,n.iter = var.sum$end)
      save(DIC,file=dicFileName)
      }
    }
  }else{
    print("Fitting Model")
    out.burn <- try(runForecastIter(j.model=j.model,variableNames=variableNames,
                                    baseNum = 50000,iterSize = 10000,effSize = 5000, maxIter=10000000,
                                    partialFile = partialFileName))

    ##Thin the data:
    out.burn <- thinMCMCOutput(out.burn)

    save(out.burn,file = outputFileName)
    print(paste("saved:",outputFileName))
  }
  }
}
