source('read_Licor.R')
source('fitA.R')
source('runLicorIter.R')
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)

library(PEcAn.photosynthesis)
#From: https://github.com/PecanProject/pecan/tree/develop/modules/photosynthesis/R

dataDirectory <- "Data/Licor_Measurements/"
leafNames <- c("B1A","B1B","B1C","B1D","B1E","B1F",
               "B2A","B2B","B2C","B2D","B2E","B2F",
               "B3A","B3B","B3C","B3D","B3E","B3F",
               "O1A","O1B","O1C","O1D","O1E","O1F",
               "O2A","O2B","O2C","O2D","O2E","O2F",
               "O3A","O3B","O3C","O3D","O3E","O3F",
               "O4A","O4B","O4C","O4D","O4E","O4F")
files <- paste0(dataDirectory,dir(path=dataDirectory,pattern="-B"))
master = lapply(files, read_Licor)
dat <- do.call("rbind", master)

A.model <- list(a.fixed = NULL, a.random = "leaf", 
                V.fixed = NULL, V.random = "leaf", 
                n.iter = 5000, match = "fname")

out <- fitA(dat,model = A.model)

save(out,file=paste0('LicorFits/',"oak_all_LicorResponseCurve_varBurn.RData"))
print(paste("Saved:",paste0('LicorFits/',"oak_all_LicorResponseCurve_varBurn.RData")))
      
load(paste0('LicorFits/',"beech_all_LicorResponseCurve_varBurn.RData"))
pdf(file="FieldWorkResponseCurvesConvergedHiearchical_Beech.pdf",height=5,width=6)
PEcAn.photosynthesis::plot_photo(dat,out,byLeaf=TRUE)
dev.off()


