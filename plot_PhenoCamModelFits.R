library('rjags')
library('runjags')
library('ecoforecastR')
library('RColorBrewer')
source('generalModels.R')
source('sharedVariables.R')
source('generalFunctions.R')
source('metFunctions.R')
Nmc <- 5000
load(paste0(finalDataDirectory,"harvard_dataFinal_withPhoto.RData"))
pdfName <- "PhenoCam_ModelFits.pdf"
pdf(file=pdfName,height=4,width=15)
par(mfrow=c(3,1))
par(mai=c(0.1,0.1,0.5,0.1))
plotFiles <- intersect(intersect(dir(path="modelFits/",pattern="harvard"),
                                 dir(path="modelFits/",pattern="varBurn")),
                       grep(list.files(path="modelFits/"),pattern="partial", invert=TRUE,value=TRUE))

for(f in 1:length(plotFiles)){
  load(paste0('modelFits/',plotFiles[f]))
  out.mat <- data.frame(as.matrix(out.burn$params))
  prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)
  
  b0 <- out.mat$b0[prow]
  b4 <- out.mat$b4[prow]
  b3 <- out.mat$b3[prow]
  pred.mat <- data.frame(as.matrix(out.burn$predict))
  ci = apply(pred.mat,2,quantile,c(0.025,0.5,0.975))
  
  plot(seq(1,length(dataFinal$p[1:122,])),dataFinal$p[1:122,],pch=20,ylim=c(0,1),main=plotFiles[f])
  ecoforecastR::ciEnvelope(seq(1,length(dataFinal$p[1:122,])),ci[1,],ci[3,],col=col.alpha("blue",0.5))
  points(seq(1,length(dataFinal$p[1:122,])),dataFinal$p[1:122,],pch=20,col="black")
  
}
dev.off()



# yearInt <- 2
#
# yearRemoved <- dataFinal$years[yearInt]
# fileName <- "modelFits/bbc1_GPP_missingYear_varBurn.RData"
# fileName <- "modelFits/bbc1_NPP_missingYear_varBurn.RData"
# print(fileName)
# load(fileName)
