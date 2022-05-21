#Estimating Chlorophyll cycling model priors from PhenoCam Fits (mostly for the random effect variances )
library(rjags)
library(runjags)
b0ALL <- numeric()
b3ALL <- numeric()
b4ALL <- numeric()
#load('bbc1_meanTemp_summer183_expBreak_slope_forecast_b3_calibration_varBurn.RData')

files <- paste0('previousPhenoCamFits/',dir(path="previousPhenoCamFits/",pattern="varBurn"))
for(f in 1:length(files)){
  load(files[f])
  out.mat <- as.data.frame(as.matrix(out.burn$params))
  rndNum <- sample.int(nrow(out.mat),1000,replace=TRUE)
  b0ALL <- c(b0ALL,mean(out.mat$b0[rndNum]))
  b3ALL <- c(b3ALL,mean(out.mat$b3[rndNum]))
  b4ALL <- c(b4ALL,mean(out.mat$b4[rndNum]))
}
mu <- 1/var(b0ALL)
var <- 500
alp <- mu**2/var #14400
bet <- mu/var #5

mu <- 1/var(b3ALL)
var <- mu*0.2
alp <- mu**2/var #700,000,000
bet <- mu/var #5

mu <- 1/var(b4ALL)
var <- mu*0.2
alp <- mu**2/var #70000
bet <- mu/var #5

