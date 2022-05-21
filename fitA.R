##' @name fitA
##' @title fitA
##' @author Mike Dietze
##' @author Xiaohui Feng
##' @author Kathryn Wheeler
##' @export
##' 
##' @param flux.data  data.frame of Licor data, concatenated by rows, and with a leading column 'fname' that is used to count the number of curves and match to covariates
##' @param cov.data   data.frame of covariate data. Column names used in formulas
##' @param model      list including at least 6 components: the fixed effects model for alpha (a.fixed) and Vcmax (V.fixed), the random effects for these (a.random, V.random), the variable used to match the gas-exchange and covariate data (match), and the number of MCMC interations (n.iter). Additional optional arguments: TPU = TRUE turns on TPU limitation; Temp == 'Bernacchi01' turns on the Bernacchi et al 2001 temperature correction. If this is turned on all parameters are estimated for 25C, otherwise no temperature correction is applied. Setting Temp = 'June2004' will turn on the June et al 2004 Funct Plant Biol temperature correction to Jmax. Note: these two corrections are not mutually exclusive, you can set Temp = c('June2004','Bernacchi2001')
##' @param hiearchical Boolean if you want hiearchical
##' 
##' Right now the fixed effects are specified as a string using the standard R lm formula syntax, but without the LHS variable (e.g. '~ SLA + chl + SLA:chl'). The tilde is optional. For random effects, the two options right now are just 'leaf' for leaf-level random effects and NULL. 'model' has a default that sets all effects to NULL (fit one curve to all data) and n.iter=1000.
##' 
fitA <- function(flux.data, cov.data = NULL, model = NULL,hiearchical) {
  #partialFile=paste0('LicorFits/',"beechTest_all_randomEffect_A_V_LicorResponseCurve_partial_varBurn.RData")
  ##  TO-DO: 
  ##  Random effects using design matrix
  ##  Model selection
  ##  output variable selection: Pred Loss, WAIC?
  ##  function to do: multiple response curves
  ##  specify priors in model object
  ##  integrate with meta-analysis
  
  library(rjags)
  
  if (is.null(model)) {
    model <- list(a.fixed = NULL, a.random = NULL, V.fixed = NULL, V.random = NULL, 
                  n.iter = 5000, match = "fname")
  }
  if(hiearchical){
    out.variables <- c("r0", "vmax0", "alpha0", "Jmax0", "cp0", "tau", "pmean","pA","Vleaf","Aleaf","Jleaf","rleaf",'cleaf')
  }else{
    out.variables <- c("r0", "Jmax0", "cp0", "tau","pA","vmax0","alpha0","pmean")
  }

  
  a.fixed <- model$a.fixed
  a.random <- model$a.random
  V.fixed <- model$V.fixed
  V.random <- model$V.random
  J.fixed <- model$J.fixed
  J.random <- model$J.random
  r.fixed <- model$r.fixed
  r.random <- model$r.random
  c.fixed <- model$c.fixed
  c.random <- model$c.random
  if (is.null(model$match)) {
    model$match <- "fname"
  }
  
  dat <- as.data.frame(flux.data)
  
  id         <- dat[, model$match]
  n.curves   <- length(unique(id))
  curve.id   <- as.numeric(as.factor(id))
  curve.code <- tapply(as.character(id), curve.id, unique)
  
  ## match between gas exchange data and covariates
  if (!is.null(cov.data)) {
    ord <- match(curve.code, as.character(cov.data[, model$match]))
    cov.data <- cov.data[ord, ]
  }
  
  ## Vcmax design matrix
  if (is.null(V.fixed)) {
    XV <- NULL
  } else {
    if (is.null(cov.data)) {
      print("Vcmax formula provided but covariate data is absent:", V.fixed)
    }
    if (length(grep("~", V.fixed)) == 0) {
      V.fixed <- paste("~", V.fixed)
    }
    XV      <- with(cov.data, model.matrix(formula(V.fixed)))
    XV.cols <- colnames(XV)
    XV.cols <- XV.cols[XV.cols != "(Intercept)"]
    XV      <- as.matrix(XV[, XV.cols])
    colnames(XV) <- XV.cols
    Vcenter <- apply(XV, 2, mean, na.rm = TRUE)
    XV      <- t(t(XV) - Vcenter)
  }
  
  ## Jmax design matrix
  if (is.null(J.fixed)) {
    XJ <- NULL
  } else {
    if (is.null(cov.data)) {
      print("Jmax formula provided but covariate data is absent:", J.fixed)
    }
    if (length(grep("~", J.fixed)) == 0) {
      J.fixed <- paste("~", J.fixed)
    }
    XJ      <- with(cov.data, model.matrix(formula(J.fixed)))
    XJ.cols <- colnames(XJ)
    XJ.cols <- XJ.cols[XJ.cols != "(Intercept)"]
    XJ      <- as.matrix(XJ[, XJ.cols])
    colnames(XJ) <- XJ.cols
    Jcenter <- apply(XJ, 2, mean, na.rm = TRUE)
    XJ      <- t(t(XJ) - Jcenter)
  }
  
  ## alpha design matrix
  if (is.null(a.fixed)) {
    Xa <- NULL
  } else {
    if (is.null(cov.data)) {
      print("alpha formula provided but covariate data is absent:", a.fixed)
    }
    a.fixed <- ifelse(length(grep("~", a.fixed)) == 0, paste("~", a.fixed), a.fixed)
    Xa      <- with(cov.data, model.matrix(formula(a.fixed)))
    Xa      <- as.matrix(Xa[, -which(colnames(Xa) == "(Intercept)")])
    acenter <- apply(Xa, 2, mean, na.rm = TRUE)
    Xa      <- t(t(Xa) - acenter)
  }
  ## cp0 design matrix
  if (is.null(c.fixed)) {
    Xc <- NULL
  } else {
    if (is.null(cov.data)) {
      print("cp0 formula provided but covariate data is absent:", c.fixed)
    }
    c.fixed <- ifelse(length(grep("~", c.fixed)) == 0, paste("~", c.fixed), c.fixed)
    Xc      <- with(cov.data, model.matrix(formula(c.fixed)))
    Xc      <- as.matrix(Xc[, -which(colnames(Xc) == "(Intercept)")])
    ccenter <- apply(Xc, 2, mean, na.rm = TRUE)
    Xc      <- t(t(Xc) - ccenter)
  }
  
  ## r0 design matrix
  if (is.null(r.fixed)) {
    Xr <- NULL
  } else {
    if (is.null(cov.data)) {
      print("r0 formula provided but covariate data is absent:", r.fixed)
    }
    r.fixed <- ifelse(length(grep("~", r.fixed)) == 0, paste("~", r.fixed), r.fixed)
    Xr      <- with(cov.data, model.matrix(formula(r.fixed)))
    Xr      <- as.matrix(Xr[, -which(colnames(Xr) == "(Intercept)")])
    rcenter <- apply(Xr, 2, mean, na.rm = TRUE)
    Xr      <- t(t(Xr) - rcenter)
  }
  
  ## Define JAGS model
  
  my.model <- "  
  model{
  ## Priors
  Jmax0 ~ dlnorm(4.7,2.7)             ## maximum electron transport rate prior
  alpha0~dnorm(0.25,100)             ##quantum yield  (mol electrons/mole photon) prior
  vmax0 ~dlnorm(4.6,2.7)             ## maximum rubisco capacity prior
  #Jmax ~ dweibull(2.0,260)          ## maximum electron transport rate prior Serbin 2012
  #alpha0 ~ dgamma(2.0,22.0)         ## quantum yield prior Serbin 2012
  #vmax0 ~ dweibull(1.7,80)          ## maximum rate of carboxylation prior Serbin 2012
  r0 ~ dlnorm(0.75,1.56)             ## leaf respiration prior
  #r ~ dweibull(2.0,6.0)             ## broad leaf respiration prior for trees
  cp0 ~ dlnorm(1.9,2.7)              ## CO2 compensation point prior
  tau ~ dgamma(0.1,0.1)
  #TPU  tpu~ dlnorm(3,2.8)             ##tpu
  ## Constants: Bernacchi et al 2001, PC&E, Table 1
  R <- 8.3144621 ## gas constant
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
  ## Constants: June et al 2004, Funct Plant Bio
  Omega <- 18
  To <- 35    ## Representative value, would benifit from spp calibration!
  ## Vcmax BETAS
  #RLEAF.V  tau.Vleaf~dgamma(0.1,0.1)          ## add random leaf effects
  #RLEAF.V  for(i in 1:nrep){                  
  #RLEAF.V   Vleaf[i]~dnorm(0,tau.Vleaf)
  #RLEAF.V  }
  ## Jmax BETAS
  #RLEAF.J  tau.Jleaf~dgamma(0.1,0.1)          ## add random leaf effects
  #RLEAF.J  for(i in 1:nrep){                  
  #RLEAF.J   Jleaf[i]~dnorm(0,tau.Jleaf)
  #RLEAF.J  }
  ## alpha BETAs
  #RLEAF.A  tau.Aleaf~dgamma(0.1,0.1)
  #RLEAF.A  for(i in 1:nrep){                  
  #RLEAF.A   Aleaf[i]~dnorm(0,tau.Aleaf)
  #RLEAF.A  }
  
  ## r0 BETAs
  #RLEAF.r  tau.rleaf~dgamma(0.1,0.1)
  #RLEAF.r  for(i in 1:nrep){                  
  #RLEAF.r   rleaf[i]~dnorm(0,tau.rleaf)
  #RLEAF.r  }
  
  ## cp0 BETAs
  #RLEAF.c  tau.cleaf~dgamma(0.1,0.1)
  #RLEAF.c  for(i in 1:nrep){                  
  #RLEAF.c   cleaf[i]~dnorm(0,tau.cleaf)
  #RLEAF.c  }
  
  for(i in 1:n) {
  r0.refT[i] <- r0 #rFORMULA
  r[i]  <- r0.refT[i] ##B01* exp(r.c - r.H/R/T[i])
  cp0.refT[i] <- cp0 #cFORMULA
  cp[i] <- cp0.refT[i] ##B01* exp(cp.c - cp.H/R/T[i])/cp.ref
  Kc.T[i] <- Kc ##B01* exp(Kc.c - Kc.H/R/T[i])/Kc.ref
  Ko.T[i] <- Ko ##B01* exp(Ko.c - Ko.H/R/T[i])/Ko.ref
  Jmax.refT[i] <- Jmax0 #JFORMULA
  Jmax[i] <- Jmax.refT[i] ##J04 * exp(-(T[i]-To)*(T[i]-To)/(Omega*Omega))
  alpha[i] <- alpha0 #AFORMULA
  al[i]<-(alpha[i]*q[i]/(sqrt(1+(alpha[i]*alpha[i]*q[i]*q[i])/(Jmax[i]*Jmax[i]))))*(pi[i]-cp[i])/(4*pi[i]+8*cp[i])    ## electron transport limited without covariates
  vmax.refT[i] <- vmax0 #VFORMULA
  vmax[i] <- vmax.refT[i] ##B01* exp(Vc.c - Vc.H/R/T[i])
  ae[i]<- vmax[i]*(pi[i]-cp[i])/(pi[i]+Kc.T[i]*(1+po/Ko.T[i]))                                                    ## maximum rubisco limited without covariates
  #TPU    ap[i]<-3*tpu                      ## phosphate limited
  pmean[i]<-min(al[i], ae[i]) - r[i]      ## predicted net photosynthesis
  an[i]~dnorm(pmean[i],tau)            ## likelihood
  pA[i] ~ dnorm(pmean[i],tau)          ## prediction
  }
  foo <- rep[1] + nrep + T[1]                ## prevent warnings
  }
  "
  
  
  ## prep data
  sel <- seq_len(nrow(dat))  #which(dat$spp == s)
  if("Tleaf" %in% names(dat)){
    if(max(dat$Tleaf) < 100){ # if Tleaf in C, convert to K
      dat$Tleaf <- dat$Tleaf + 273.15
    } 
  } else if (!"Tleaf" %in% names(dat)) {
    dat$Tleaf <- 25 + 273.15 ## if no Tleaf, assume 25C in Kelvin
    warning("No Leaf Temperature provided, setting to 25C\n",
            "To change add a column named Tleaf to flux.data data frame")
  }
  
  mydat <- list(an = dat$Photo[sel], 
                pi = dat$Ci[sel], 
                q = dat$PARi[sel],
                T = dat$Tleaf, 
                n = length(sel), Kc = 46, 
                Ko = 22000, 
                po = 21000, 
                rep = curve.id, 
                nrep = n.curves)
  #  Kc<-46                          ## Michaelis constant CO2 (Pa)
  #  Ko<-33000                       ## Michaelis constant O2  (Pa)
  #  po<-21000                       ## partial pressure of O2  (Pa)
  
  ## TPU Limitation
  if ("TPU" %in% names(model)) {
    if (model$TPU == TRUE) {
      my.model <- gsub(pattern = "#TPU", " ", my.model)
      out.variables <- c(out.variables, "tpu")
    }
  }
  
  ## Temperature scaling
  Vformula <- NULL
  if ("Temp" %in% names(model)) {
    if ("Bernacchi01" %in% model$Temp) {
      my.model <- gsub(pattern = "##B01", " ", my.model)
    }
    if ("June2004" %in% model$Temp) {
      my.model <- gsub(pattern = "##J04", " ", my.model)
    }
  }
  
  ## VCmax Formulas
  Vformula <- NULL
  if ("leaf" %in% V.random) {
    Vformula <- " + Vleaf[rep[i]]"
    my.model <- gsub(pattern = "#RLEAF.V", " ", my.model)
    out.variables <- c(out.variables, "tau.Vleaf")
  }
  
  if (!is.null(XV)) {
    Vnames <- gsub(" ", "_", colnames(XV))
    Vformula <- paste(Vformula,
                      paste0("+ betaV", Vnames, "*XV[rep[i],", seq_len(ncol(XV)), "]", collapse = " "))
    Vpriors <- paste0("     betaV", Vnames, "~dnorm(0,0.001)", collapse = "\n")
    my.model <- sub(pattern = "## Vcmax BETAS", Vpriors, my.model)
    mydat[["XV"]] <- XV
    out.variables <- c(out.variables, paste0("betaV", Vnames))
  }
  if (!is.null(Vformula)) {
    my.model <- sub(pattern = "#VFORMULA", Vformula, my.model)
  } 
  
  ## Jmax Formulas
  Jformula <- NULL
  if ("leaf" %in% J.random) {
    Jformula <- " + Jleaf[rep[i]]"
    my.model <- gsub(pattern = "#RLEAF.J", " ", my.model)
    out.variables <- c(out.variables, "tau.Jleaf")
  }
  
  if (!is.null(XJ)) {
    Jnames <- gsub(" ", "_", colnames(XJ))
    Jformula <- paste(Jformula,
                      paste0("+ betaJ", Jnames, "*XJ[rep[i],", seq_len(ncol(XJ)), "]", collapse = " "))
    Jpriors <- paste0("     betaJ", Jnames, "~dnorm(0,0.001)", collapse = "\n")
    my.model <- sub(pattern = "## Jmax BETAS", Jpriors, my.model)
    mydat[["XJ"]] <- XJ
    out.variables <- c(out.variables, paste0("betaJ", Jnames))
  }
  if (!is.null(Jformula)) {
    my.model <- sub(pattern = "#JFORMULA", Jformula, my.model)
  } 
  
  ## alpha Formulas
  Aformula <- NULL
  if ("leaf" %in% a.random) {
    Aformula <- " + Aleaf[rep[i]]"
    my.model <- gsub(pattern = "#RLEAF.A", "", my.model)
    out.variables <- c(out.variables, "tau.Aleaf")
  }
  
  if (!is.null(Xa)) {
    Anames <- gsub(" ", "_", colnames(Xa))
    Aformula <- paste(Aformula, paste0("+ betaA", Anames, "*Xa[rep[i],", 1:ncol(Xa), 
                                       "]", collapse = " "))
    apriors <- paste0("betaA", Anames, "~dnorm(0,0.001)", collapse = "\n")
    my.model <- sub(pattern = "## alpha BETAs", apriors, my.model)
    mydat[["Xa"]] <- Xa
    out.variables <- c(out.variables, paste0("betaA", Anames))
  }
  if (!is.null(Aformula)) {
    my.model <- sub(pattern = "#AFORMULA", Aformula, my.model)
  }
  
  ## cp0 Formulas
  cformula <- NULL
  if ("leaf" %in% c.random) {
    cformula <- " + cleaf[rep[i]]"
    my.model <- gsub(pattern = "#RLEAF.c", "", my.model)
    out.variables <- c(out.variables, "tau.cleaf")
  }
  
  if (!is.null(Xc)) {
    cnames <- gsub(" ", "_", colnames(Xc))
    cformula <- paste(cformula, paste0("+ betac", cnames, "*Xc[rep[i],", 1:ncol(Xc), 
                                       "]", collapse = " "))
    cpriors <- paste0("betac", cnames, "~dnorm(0,0.001)", collapse = "\n")
    my.model <- sub(pattern = "## cp0 BETAs", cpriors, my.model)
    mydat[["Xc"]] <- Xc
    out.variables <- c(out.variables, paste0("betac", cnames))
  }
  if (!is.null(cformula)) {
    my.model <- sub(pattern = "#cFORMULA", cformula, my.model)
  }
  
  ## r0 Formulas
  rformula <- NULL
  if ("leaf" %in% r.random) {
    rformula <- " + rleaf[rep[i]]"
    my.model <- gsub(pattern = "#RLEAF.r", "", my.model)
    out.variables <- c(out.variables, "tau.rleaf")
  }
  
  if (!is.null(Xr)) {
    rnames <- gsub(" ", "_", colnames(Xr))
    rformula <- paste(rformula, paste0("+ betar", rnames, "*Xr[rep[i],", 1:ncol(Xr), 
                                       "]", collapse = " "))
    rpriors <- paste0("betar", rnames, "~dnorm(0,0.001)", collapse = "\n")
    my.model <- sub(pattern = "## r0 BETAs", rpriors, my.model)
    mydat[["Xr"]] <- Xr
    out.variables <- c(out.variables, paste0("betar", rnames))
  }
  if (!is.null(rformula)) {
    my.model <- sub(pattern = "#rFORMULA", rformula, my.model)
  }
  
  ## Define initial conditions
  init <- list()
  init[[1]] <- list(r0 = 1.2, vmax0 = 39, alpha0 = 0.25, tau = 10, cp0 = 6, Jmax0 = 80)  ## tau.Vleaf=30,beta1=4, beta2=1,beta5=3,tau.Vmon=10,tpu=10,
  init[[2]] <- list(r0 = 1, vmax0 = 100, alpha0 = 0.2, tau = 20, cp0 = 4, Jmax0 = 150)  ##tau.Vleaf=20,beta1=1,beta2=1,beta5=-1,tau.Vmon=20,tpu=13,
  init[[3]] <- list(r0 = 2, vmax0 = 60, alpha0 = 0.28, tau = 20, cp0 = 5, Jmax0 = 60) 
  init[[4]] <- list(r0 = 2, vmax0 = 60, alpha0 = 0.28, tau = 20, cp0 = 5, Jmax0 = 60) 
  init[[5]] <- list(r0 = 1.2, vmax0 = 39, alpha0 = 0.25, tau = 10, cp0 = 6, Jmax0 = 80)
##tau.Vleaf=100,beta1=1,beta2=2,beta5=2,tau.Vmon=3,tpu=20,
  print(out.variables)
  print(my.model)
  mc3 <- jags.model(file = textConnection(my.model), data = mydat, inits = init, n.chains = 5,n.adapt = 3000)
  
  #mc3.out <- coda.samples(model = mc3, variable.names = out.variables, n.iter = model$n.iter)
  out <- runLicorIter(j.model=mc3,variableNames=out.variables,
                                  baseNum = 50000,iterSize = 50000,effSize = 5000)
  
  ## split output
  # out         <- list(params = NULL, predict = NULL, model = my.model)
  # mfit        <- as.matrix(mc3.out, chains = TRUE)
  # pred.cols   <- union(grep("pA", colnames(mfit)), grep("pmean", colnames(mfit)))
  # chain.col   <- which(colnames(mfit) == "CHAIN")
  # out$predict <- mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
  # out$params  <- mat2mcmc.list(mfit[, -pred.cols])
  out <- thinMCMCOutput(out)
  return(out)
} # fitA
