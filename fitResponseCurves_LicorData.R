source('read_Licor.R')
source('fitA.R')
source('runLicorIter.R')
source('sharedVariables.R')
source('generalFunctions.R')
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(PEcAn.photosynthesis)#From: https://github.com/PecanProject/pecan/tree/develop/modules/photosynthesis/R

hiearchical <- TRUE
fixedCCIeffect <- FALSE
treeSpecies <- "beech"
outputFileName <- paste0('LicorFits/',treeSpecies,"_QAQC_August_allRE_LicorResponseCurve_varBurn.RData")
if(treeSpecies=="oak"){
  Allfiles <- paste0(LicorDataDirectory,intersect(dir(path=LicorDataDirectory,pattern="-O"),
                                                  dir(path=LicorDataDirectory,pattern="AUG")))
  if(fixedCCIeffect){
    Allfiles <- paste0(LicorDataDirectory,dir(path=LicorDataDirectory,pattern="-O"))
    outputFileName <- paste0('LicorFits/',treeSpecies,"_QAQC_All_VJfixed_LicorResponseCurve_varBurn.RData")
  }
}else if(treeSpecies == "beech"){
  Allfiles <- paste0(LicorDataDirectory,intersect(dir(path=LicorDataDirectory,pattern="-B"),
                                                  dir(path=LicorDataDirectory,pattern="AUG")))
  if(fixedCCIeffect){
    Allfiles <- paste0(LicorDataDirectory,dir(path=LicorDataDirectory,pattern="-O"))
    outputFileName <- paste0('LicorFits/',treeSpecies,"_QAQC_All_VJfixed_LicorResponseCurve_varBurn.RData")
  }
}
print(outputFileName)

for(i in 1:length(Allfiles)){
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

if(hiearchical){
  A.model <- list(a.fixed = NULL, a.random = "leaf", 
                  V.fixed = NULL, V.random = "leaf",
                  J.fixed = NULL, J.random = "leaf",
                  r.fixed = NULL, r.random = "leaf",
                  c.fixed = NULL, c.random = "leaf",
                  n.iter = 5000, match = "fname")
  if(fixedCCIeffect){
    usedFiles <- unique(allDat$fname)
    cov.dat <- as.data.frame(matrix(nrow=length(usedFiles),ncol=1))
    cov.dat[,1] <- usedFiles
    colnames(cov.dat) <- "fname"
    
    CCIdat <- read.csv("CCIforLicorMeasurements.csv")
    cov.dat <- merge(cov.dat,CCIdat,'fname')[,1:2]
    
    colnames(cov.dat) <- c("fname","CCI_mean")
    A.model <- list(a.fixed = NULL, a.random = "leaf", 
                    V.fixed = "CCI_mean", V.random = NULL,
                    J.fixed = "CCI_mean", J.random = NULL,
                    r.fixed = NULL, r.random = "leaf",
                    c.fixed = NULL, c.random = "leaf",
                    n.iter = 5000, match = "fname")
    out <- fitA(flux.data=allDat,model = A.model,cov.data = cov.dat,hiearchical = hiearchical)
  }else{
    out <- fitA(allDat,model = A.model,hiearchical = hiearchical)
  }
}else{
  out <- fitA(allDat,model = A.model,hiearchical = hiearchical)
}

save(out,file=outputFileName)
print(paste("Saved:", outputFileName))

# load(paste0('LicorFits/',"beech_all_LicorResponseCurve_varBurn.RData"))
# pdf(file="FieldWorkResponseCurvesConvergedHiearchical_Beech.pdf",height=5,width=6)
# PEcAn.photosynthesis::plot_photo(dat,out,byLeaf=TRUE)
# dev.off()



