#PhenoCam Models ----
generalModel_PC_Tair_D = "
model {
    ### Data Models for complete years
    for(yr in 1:(N)){
    for(i in 1:n){
    p[i,yr] ~ dnorm(x[i,yr],p.PC)
    }
    }
    
    #### Process Model
    for(yr in 1:(N)){
    for(i in 2:n){
    
    xmu[i,yr] <- max(min((x[(i-1),yr] + b4 * x[(i-1),yr]) + max(0,(b0 + b3 * TowerTair[i,yr] * D[i,yr])),x[1,yr]),0)
    x[i,yr] ~ dnorm(xmu[i,yr],p.proc)
    }
    }
    
    #### Priors
    for(yr in 1:N){ ##Initial Conditions
    x[1,yr] ~ dbeta(x1.a[yr],x1.b[yr]) I(0.001,0.999)
    }
    p.PC ~ dgamma(s1.PC,s2.PC)
    p.proc ~ dgamma(s1.proc,s2.proc)
    b0 ~ dunif(b0_lower,b0_upper)
    b3 ~ dunif(b3_lower,b3_upper)
    b4 ~ dunif(b4_lower,b4_upper)
    
  }
    "

generalModel_PC_Cov = "
model {
    ### Data Models for complete years
    for(yr in 1:(N)){
    for(i in 1:n){
    p[i,yr] ~ dnorm(x[i,yr],p.PC)
    }
    }
    
    #### Process Model
    for(yr in 1:(N)){
    for(i in 2:n){
    
    xmu[i,yr] <- max(min((x[(i-1),yr] + b4 * x[(i-1),yr]) + max(0,(b0 + b3 * Cov[i,yr])),x[1,yr]),0)
    x[i,yr] ~ dnorm(xmu[i,yr],p.proc)
    }
    }
    
    #### Priors
    for(yr in 1:N){ ##Initial Conditions
    x[1,yr] ~ dbeta(x1.a[yr],x1.b[yr]) I(0.001,0.999)
    }
    p.PC ~ dgamma(s1.PC,s2.PC)
    p.proc ~ dgamma(s1.proc,s2.proc)
    b0 ~ dunif(b0_lower,b0_upper)
    b3 ~ dunif(b3_lower,b3_upper)
    b4 ~ dunif(b4_lower,b4_upper)
    
  }
    "

generalModel_PC_CovPlusTair_D = "
model {
    ### Data Models for complete years
    for(yr in 1:(N)){
    for(i in 1:n){
    p[i,yr] ~ dnorm(x[i,yr],p.PC)
    }
    }
    
    #### Process Model
    for(yr in 1:(N)){
    for(i in 2:n){
    
    xmu[i,yr] <- max(min((x[(i-1),yr] + b4 * x[(i-1),yr]) + max(0,(b0 + b2 * Cov[i,yr] + b3 * TowerTair[i,yr] * D[i,yr])),x[1,yr]),0)
    x[i,yr] ~ dnorm(xmu[i,yr],p.proc)
    }
    }
    
    #### Priors
    for(yr in 1:N){ ##Initial Conditions
    x[1,yr] ~ dbeta(x1.a[yr],x1.b[yr]) I(0.001,0.999)
    }
    p.PC ~ dgamma(s1.PC,s2.PC)
    p.proc ~ dgamma(s1.proc,s2.proc)
    b0 ~ dunif(b0_lower,b0_upper)
    b2 ~ dunif(b3_lower,b3_upper)
    b3 ~ dunif(b3_lower,b3_upper)
    b4 ~ dunif(b4_lower,b4_upper)
    
  }
    "

#Beech Leaves Models ----
generalModel_B_Tair_D = "
model {
    ### Data Models for complete years
    for(l in 1:3){ #loop over leaves
    for(i in 1:n){
    CCI_means[l,i,1] ~ dnorm(x[l,i,1],CCI_precs[l,i,1]) T(0.0001,0.9999)
    CCI_means[l,i,3] ~ dnorm(x[l,i,3],CCI_precs[l,i,3]) T(0.0001,0.9999)
    }
    }
    
    for(l in 1:2){ #loop over leaves
    for(i in 1:n){
    CCI_means[l,i,2] ~ dnorm(x[l,i,2],CCI_precs[l,i,2]) T(0.0001,0.9999)
    }
    }
    
    #### Process Model
    for(t in 1:treeN){
    b0_tree[t] ~ dnorm(b0,b0_prec) T(0.0001,0.9999)
    b3_tree[t] ~ dnorm(b3,b3_prec) T(0.0001,0.9999)
    b4_tree[t] ~ dnorm(b4,b4_prec) T(0.0001,0.9999)
    }
    
    for(t in c(1,3)){
    for(l in 1:3){ #loop over leaves
    b0_leaf[l,t] ~ dnorm(b0_tree[t],b0_prec_t) T(0.0001,0.9999)
    b3_leaf[l,t] ~ dnorm(b3_tree[t],b3_prec_t) T(0.0001,0.9999)
    b4_leaf[l,t] ~ dnorm(b4_tree[t],b4_prec_t) T(0.0001,0.9999)
    for(i in 2:n){
    xmu[l,i,t] <- max(min((x[l,(i-1),t] + -1*b4 * x[l,(i-1),t]) + max(0,(-1*b0_leaf[l,t] + b3_leaf[l,t] * Tair[l,i,t] * D[i])),x[l,1,t]),0)
    x[l,i,t] ~ dnorm(xmu[l,i,t],p.proc)
    }
    }
    }
    
    for(l in 1:2){ #loop over leaves
    b0_leaf[l,2] ~ dnorm(b0_tree[2],b0_prec_t) T(0.0001,0.9999)
    b3_leaf[l,2] ~ dnorm(b3_tree[2],b3_prec_t) T(0.0001,0.9999)
    b4_leaf[l,2] ~ dnorm(b4_tree[2],b4_prec_t) T(0.0001,0.9999)
    for(i in 2:n){
    xmu[l,i,2] <- max(min((x[l,(i-1),2] + -1*b4 * x[l,(i-1),2]) + max(0,(-1*b0_leaf[l,2] + b3_leaf[l,2] * Tair[l,i,2] * D[i])),x[l,1,2]),0)
    x[l,i,2] ~ dnorm(xmu[l,i,2],p.proc)
    }
    }
    
    #### Priors
    for(l in 1:3){ #loop over leaves
    x[l,1,1] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    x[l,1,3] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    }
    for(l in 1:2){ #loop over leaves
    x[l,1,2] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    }
  
    p.proc ~ dgamma(s1.proc,s2.proc)
    b0_prec ~ dgamma(14400,5) #Priors based on looking at variation of means in fitted PhenoCam data 
    b4_prec ~ dgamma(70000,5)
    b3_prec ~ dgamma(700000000,5)
    b0_prec_t <- b0_prec*2
    b3_prec_t <- b3_prec*2
    b4_prec_t <- b4_prec*2
    b0 ~ dbeta(b0.a,b0.b) 
    b3 ~ dbeta(b3.a,b3.b) 
    b4 ~ dbeta(b4.a,b4.b) 
  }
    "
generalModel_B_Cov = "
model {
    ### Data Models for complete years
    for(l in 1:3){ #loop over leaves
    for(i in 1:n){
    CCI_means[l,i,1] ~ dnorm(x[l,i,1],CCI_precs[l,i,1]) T(0.0001,0.9999)
    CCI_means[l,i,3] ~ dnorm(x[l,i,3],CCI_precs[l,i,3]) T(0.0001,0.9999)
    }
    }
    
    for(l in 1:2){ #loop over leaves
    for(i in 1:n){
    CCI_means[l,i,2] ~ dnorm(x[l,i,2],CCI_precs[l,i,2]) T(0.0001,0.9999)
    }
    }

    
    #### Process Model
    for(t in 1:treeN){
    b0_tree[t] ~ dnorm(b0,b0_prec) T(0.0001,0.9999)
    b3_tree[t] ~ dnorm(b3,b3_prec) T(0.0001,0.9999)
    b4_tree[t] ~ dnorm(b4,b4_prec) T(0.0001,0.9999)
    }
    
    for(t in c(1,3)){
    for(l in 1:3){ #loop over leaves
    b0_leaf[l,t] ~ dnorm(b0_tree[t],b0_prec_t) T(0.0001,0.9999)
    b3_leaf[l,t] ~ dnorm(b3_tree[t],b3_prec_t) T(0.0001,0.9999)
    b4_leaf[l,t] ~ dnorm(b4_tree[t],b4_prec_t) T(0.0001,0.9999)
    for(i in 2:n){
    xmu[l,i,t] <- max(min((x[l,(i-1),t] + -1*b4 * x[l,(i-1),t]) + max(0,(-1*b0_leaf[l,t] + b3_leaf[l,t] * Cov[l,i,t])),x[l,1,t]),0)
    x[l,i,t] ~ dnorm(xmu[l,i,t],p.proc)
    }
    }
    }
    
    for(l in 1:2){ #loop over leaves
    b0_leaf[l,2] ~ dnorm(b0_tree[2],b0_prec_t) T(0.0001,0.9999)
    b3_leaf[l,2] ~ dnorm(b3_tree[2],b3_prec_t) T(0.0001,0.9999)
    b4_leaf[l,2] ~ dnorm(b4_tree[2],b4_prec_t) T(0.0001,0.9999)
    for(i in 2:n){
    xmu[l,i,2] <- max(min((x[l,(i-1),2] + -1*b4 * x[l,(i-1),2]) + max(0,(-1*b0_leaf[l,2] + b3_leaf[l,2] * Cov[l,i,2])),x[l,1,2]),0)
    x[l,i,2] ~ dnorm(xmu[l,i,2],p.proc)
    }
    }
    
    #### Priors
    for(l in 1:3){ #loop over leaves
    x[l,1,1] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    x[l,1,3] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    }
    for(l in 1:2){ #loop over leaves
    x[l,1,2] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    }
  
    p.proc ~ dgamma(s1.proc,s2.proc)
    b0_prec ~ dgamma(14400,5) #Priors based on looking at variation of means in fitted PhenoCam data 
    b4_prec ~ dgamma(70000,5)
    b3_prec ~ dgamma(700000000,5)
    b0_prec_t <- b0_prec*2
    b3_prec_t <- b3_prec*2
    b4_prec_t <- b4_prec*2
    b0 ~ dbeta(b0.a,b0.b) 
    b3 ~ dbeta(b3.a,b3.b) 
    b4 ~ dbeta(b4.a,b4.b) 
  }
    "
generalModel_B_CovPlusTair_D = "
model {
    ### Data Models for complete years
    for(l in 1:3){ #loop over leaves
    for(i in 1:n){
    CCI_means[l,i,1] ~ dnorm(x[l,i,1],CCI_precs[l,i,1]) T(0.0001,0.9999)
    CCI_means[l,i,3] ~ dnorm(x[l,i,3],CCI_precs[l,i,3]) T(0.0001,0.9999)
    }
    }
    
    for(l in 1:2){ #loop over leaves
    for(i in 1:n){
    CCI_means[l,i,2] ~ dnorm(x[l,i,2],CCI_precs[l,i,2]) T(0.0001,0.9999)
    }
    }
    
    #### Process Model
    for(t in 1:treeN){
    b0_tree[t] ~ dnorm(b0,b0_prec) T(0.0001,0.9999)
    b3_tree[t] ~ dnorm(b3,b3_prec) T(0.0001,0.9999)
    b2_tree[t] ~ dnorm(b2,b2_prec) T(0.0001,0.9999)
    b4_tree[t] ~ dnorm(b4,b4_prec) T(0.0001,0.9999)
    }
    
    for(t in c(1,3)){
    for(l in 1:3){ #loop over leaves
    b0_leaf[l,t] ~ dnorm(b0_tree[t],b0_prec_t) T(0.0001,0.9999)
    b2_leaf[l,t] ~ dnorm(b2_tree[t],b2_prec_t) T(0.0001,0.9999)
    b3_leaf[l,t] ~ dnorm(b3_tree[t],b3_prec_t) T(0.0001,0.9999)
    b4_leaf[l,t] ~ dnorm(b4_tree[t],b4_prec_t) T(0.0001,0.9999)
    for(i in 2:n){
    xmu[l,i,t] <- max(min((x[l,(i-1),t] + -1*b4 * x[l,(i-1),t]) + max(0,(-1*b0_leaf[l,t] + b2_leaf[l,t] * Cov[l,i,t] +b3_leaf[l,t] * Tair[l,i,t] * D[i])),x[l,1,t]),0)
    x[l,i,t] ~ dnorm(xmu[l,i,t],p.proc)
    }
    }
    }
    
    for(l in 1:2){ #loop over leaves
    b0_leaf[l,2] ~ dnorm(b0_tree[2],b0_prec_t) T(0.0001,0.9999)
    b2_leaf[l,2] ~ dnorm(b2_tree[2],b2_prec_t) T(0.0001,0.9999)
    b3_leaf[l,2] ~ dnorm(b3_tree[2],b3_prec_t) T(0.0001,0.9999)
    b4_leaf[l,2] ~ dnorm(b4_tree[2],b4_prec_t) T(0.0001,0.9999)
    for(i in 2:n){
    xmu[l,i,2] <- max(min((x[l,(i-1),2] + -1*b4 * x[l,(i-1),2]) + max(0,(-1*b0_leaf[l,2] + b2_leaf[l,2] * Cov[l,i,2] + b3_leaf[l,2] * Tair[l,i,2] * D[i])),x[l,1,2]),0)
    x[l,i,2] ~ dnorm(xmu[l,i,2],p.proc)
    }
    }
    
    #### Priors
    for(l in 1:3){ #loop over leaves
    x[l,1,1] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    x[l,1,3] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    }
    for(l in 1:2){ #loop over leaves
    x[l,1,2] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    }
  
    p.proc ~ dgamma(s1.proc,s2.proc)
    b0_prec ~ dgamma(14400,5) #Priors based on looking at variation of means in fitted PhenoCam data 
    b4_prec ~ dgamma(70000,5)
    b3_prec ~ dgamma(700000000,5)
    b2_prec ~ dgamma(70000,5)
    b0_prec_t <- b0_prec*2
    b3_prec_t <- b3_prec*2
    b4_prec_t <- b4_prec*2
    b2_prec_t <- b2_prec*2
    b0 ~ dbeta(b0.a,b0.b) 
    b3 ~ dbeta(b3.a,b3.b) 
    b4 ~ dbeta(b4.a,b4.b) 
    b2 ~ dbeta(b2.a,b2.b)
  }
    "

#Oak Leaves Models ----
generalModel_O_Tair_D = "
model {
    ### Data Models for complete years
    for(l in 1:3){ #loop over leaves
    for(i in 1:n){
    CCI_means[l,i,1] ~ dnorm(x[l,i,1],CCI_precs[l,i,1]) T(0.0001,0.9999)
    CCI_means[l,i,2] ~ dnorm(x[l,i,2],CCI_precs[l,i,2]) T(0.0001,0.9999)
    }
    }
    
    for(l in 1:2){ #loop over leaves
    for(i in 1:n){
    CCI_means[l,i,3] ~ dnorm(x[l,i,3],CCI_precs[l,i,3]) T(0.0001,0.9999)
    CCI_means[l,i,4] ~ dnorm(x[l,i,4],CCI_precs[l,i,4]) T(0.0001,0.9999)
    }
    }

    
    #### Process Model
    for(t in 1:treeN){
    b0_tree[t] ~ dnorm(b0,b0_prec) T(0.0001,0.9999)
    b3_tree[t] ~ dnorm(b3,b3_prec) T(0.0001,0.9999)
    b4_tree[t] ~ dnorm(b4,b4_prec) T(0.0001,0.9999)
    }
    
    for(t in 1:2){
    for(l in 1:3){ #loop over leaves
    b0_leaf[l,t] ~ dnorm(b0_tree[t],b0_prec_t) T(0.0001,0.9999)
    b3_leaf[l,t] ~ dnorm(b3_tree[t],b3_prec_t) T(0.0001,0.9999)
    b4_leaf[l,t] ~ dnorm(b4_tree[t],b4_prec_t) T(0.0001,0.9999)
    for(i in 2:n){
    xmu[l,i,t] <- max(min((x[l,(i-1),t] + -1*b4 * x[l,(i-1),t]) + max(0,(-1*b0_leaf[l,t] + b3_leaf[l,t] * Tair[l,i,t] * D[i])),x[l,1,t]),0)
    x[l,i,t] ~ dnorm(xmu[l,i,t],p.proc)
    }
    }
    }
    
    for(t in 3:4){
    for(l in 1:2){ #loop over leaves
    b0_leaf[l,t] ~ dnorm(b0_tree[t],b0_prec_t) T(0.0001,0.9999)
    b3_leaf[l,t] ~ dnorm(b3_tree[t],b3_prec_t) T(0.0001,0.9999)
    b4_leaf[l,t] ~ dnorm(b4_tree[t],b4_prec_t) T(0.0001,0.9999)
    for(i in 2:n){
    xmu[l,i,t] <- max(min((x[l,(i-1),t] + -1*b4 * x[l,(i-1),t]) + max(0,(-1*b0_leaf[l,t] + b3_leaf[l,t] * Tair[l,i,t] * D[i])),x[l,1,t]),0)
    x[l,i,t] ~ dnorm(xmu[l,i,t],p.proc)
    }
    }
    }
    
    #### Priors
    for(l in 1:3){ #loop over leaves
    x[l,1,1] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    x[l,1,2] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    }
    for(l in 1:2){ #loop over leaves
    x[l,1,3] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    x[l,1,4] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    }
  
    p.proc ~ dgamma(s1.proc,s2.proc)
    b0_prec ~ dgamma(14400,5) #Priors based on looking at variation of means in fitted PhenoCam data 
    b4_prec ~ dgamma(70000,5)
    b3_prec ~ dgamma(700000000,5)
    b0_prec_t <- b0_prec*2
    b3_prec_t <- b3_prec*2
    b4_prec_t <- b4_prec*2
    b0 ~ dbeta(b0.a,b0.b) 
    b3 ~ dbeta(b3.a,b3.b) 
    b4 ~ dbeta(b4.a,b4.b) 
  }
    "

generalModel_O_Cov = "
model {
    ### Data Models for complete years
    for(l in 1:3){ #loop over leaves
    for(i in 1:n){
    CCI_means[l,i,1] ~ dnorm(x[l,i,1],CCI_precs[l,i,1]) T(0.0001,0.9999)
    CCI_means[l,i,2] ~ dnorm(x[l,i,2],CCI_precs[l,i,2]) T(0.0001,0.9999)
    }
    }
    
    for(l in 1:2){ #loop over leaves
    for(i in 1:n){
    CCI_means[l,i,3] ~ dnorm(x[l,i,3],CCI_precs[l,i,3]) T(0.0001,0.9999)
    CCI_means[l,i,4] ~ dnorm(x[l,i,4],CCI_precs[l,i,4]) T(0.0001,0.9999)
    }
    }

    
    #### Process Model
    for(t in 1:treeN){
    b0_tree[t] ~ dnorm(b0,b0_prec) T(0.0001,0.9999)
    b3_tree[t] ~ dnorm(b3,b3_prec) T(0.0001,0.9999)
    b4_tree[t] ~ dnorm(b4,b4_prec) T(0.0001,0.9999)
    }
    
    for(t in 1:2){
    for(l in 1:3){ #loop over leaves
    b0_leaf[l,t] ~ dnorm(b0_tree[t],b0_prec_t) T(0.0001,0.9999)
    b3_leaf[l,t] ~ dnorm(b3_tree[t],b3_prec_t) T(0.0001,0.9999)
    b4_leaf[l,t] ~ dnorm(b4_tree[t],b4_prec_t) T(0.0001,0.9999)
    for(i in 2:n){
    xmu[l,i,t] <- max(min((x[l,(i-1),t] + -1*b4 * x[l,(i-1),t]) + max(0,(-1*b0_leaf[l,t] + b3_leaf[l,t] * Cov[l,i,t])),x[l,1,t]),0)
    x[l,i,t] ~ dnorm(xmu[l,i,t],p.proc)
    }
    }
    }
    
    for(t in 3:4){
    for(l in 1:2){ #loop over leaves
    b0_leaf[l,t] ~ dnorm(b0_tree[t],b0_prec_t) T(0.0001,0.9999)
    b3_leaf[l,t] ~ dnorm(b3_tree[t],b3_prec_t) T(0.0001,0.9999)
    b4_leaf[l,t] ~ dnorm(b4_tree[t],b4_prec_t) T(0.0001,0.9999)
    for(i in 2:n){
    xmu[l,i,t] <- max(min((x[l,(i-1),t] + -1*b4 * x[l,(i-1),t]) + max(0,(-1*b0_leaf[l,t] + b3_leaf[l,t] * Cov[l,i,t])),x[l,1,t]),0)
    x[l,i,t] ~ dnorm(xmu[l,i,t],p.proc)
    }
    }
    }
    
    #### Priors
    for(l in 1:3){ #loop over leaves
    x[l,1,1] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    x[l,1,2] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    }
    for(l in 1:2){ #loop over leaves
    x[l,1,3] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    x[l,1,4] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    }
  
    p.proc ~ dgamma(s1.proc,s2.proc)
    b0_prec ~ dgamma(14400,5) #Priors based on looking at variation of means in fitted PhenoCam data 
    b4_prec ~ dgamma(70000,5)
    b3_prec ~ dgamma(700000000,5)
    b0_prec_t <- b0_prec*2
    b3_prec_t <- b3_prec*2
    b4_prec_t <- b4_prec*2
    b0 ~ dbeta(b0.a,b0.b) 
    b3 ~ dbeta(b3.a,b3.b) 
    b4 ~ dbeta(b4.a,b4.b) 
  }
    "

generalModel_O_CovPlusTair_D = "
model {
    ### Data Models for complete years
    for(l in 1:3){ #loop over leaves
    for(i in 1:n){
    CCI_means[l,i,1] ~ dnorm(x[l,i,1],CCI_precs[l,i,1]) T(0.0001,0.9999)
    CCI_means[l,i,2] ~ dnorm(x[l,i,2],CCI_precs[l,i,2]) T(0.0001,0.9999)
    }
    }
    
    for(l in 1:2){ #loop over leaves
    for(i in 1:n){
    CCI_means[l,i,3] ~ dnorm(x[l,i,3],CCI_precs[l,i,3]) T(0.0001,0.9999)
    CCI_means[l,i,4] ~ dnorm(x[l,i,4],CCI_precs[l,i,4]) T(0.0001,0.9999)
    }
    }

    
    #### Process Model
    for(t in 1:treeN){
    b0_tree[t] ~ dnorm(b0,b0_prec) T(0.0001,0.9999)
    b2_tree[t] ~ dnorm(b2,b2_prec) T(0.0001,0.9999)
    b3_tree[t] ~ dnorm(b3,b3_prec) T(0.0001,0.9999)
    b4_tree[t] ~ dnorm(b4,b4_prec) T(0.0001,0.9999)
    }
    
    for(t in 1:2){
    for(l in 1:3){ #loop over leaves
    b0_leaf[l,t] ~ dnorm(b0_tree[t],b0_prec_t) T(0.0001,0.9999)
    b2_leaf[l,t] ~ dnorm(b2_tree[t],b2_prec_t) T(0.0001,0.9999)
    b3_leaf[l,t] ~ dnorm(b3_tree[t],b3_prec_t) T(0.0001,0.9999)
    b4_leaf[l,t] ~ dnorm(b4_tree[t],b4_prec_t) T(0.0001,0.9999)
    for(i in 2:n){
    xmu[l,i,t] <- max(min((x[l,(i-1),t] + -1*b4 * x[l,(i-1),t]) + max(0,(-1*b0_leaf[l,t] + b2_leaf[l,t] * Cov[l,i,t] +b3_leaf[l,t] * Tair[l,i,t] * D[i])),x[l,1,t]),0)
    x[l,i,t] ~ dnorm(xmu[l,i,t],p.proc)
    }
    }
    }
    
    for(t in 3:4){
    for(l in 1:2){ #loop over leaves
    b0_leaf[l,t] ~ dnorm(b0_tree[t],b0_prec_t) T(0.0001,0.9999)
    b2_leaf[l,t] ~ dnorm(b2_tree[t],b2_prec_t) T(0.0001,0.9999)
    b3_leaf[l,t] ~ dnorm(b3_tree[t],b3_prec_t) T(0.0001,0.9999)
    b4_leaf[l,t] ~ dnorm(b4_tree[t],b4_prec_t) T(0.0001,0.9999)
    for(i in 2:n){
    xmu[l,i,t] <- max(min((x[l,(i-1),t] + -1*b4 * x[l,(i-1),t]) + max(0,(-1*b0_leaf[l,t] + b2_leaf[l,t] * Cov[l,i,t]+ b3_leaf[l,t] * Tair[l,i,t] * D[i])),x[l,1,t]),0)
    x[l,i,t] ~ dnorm(xmu[l,i,t],p.proc)
    }
    }
    }
    
    #### Priors
    for(l in 1:3){ #loop over leaves
    x[l,1,1] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    x[l,1,2] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    }
    for(l in 1:2){ #loop over leaves
    x[l,1,3] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    x[l,1,4] ~ dbeta(x1.a,x1.b) I(0.001,0.999)
    }
  
    p.proc ~ dgamma(s1.proc,s2.proc)
    b0_prec ~ dgamma(14400,5) #Priors based on looking at variation of means in fitted PhenoCam data 
    b2_prec ~ dgamma(70000,5)
    b4_prec ~ dgamma(70000,5)
    b3_prec ~ dgamma(700000000,5)
    b0_prec_t <- b0_prec*2
    b3_prec_t <- b3_prec*2
    b2_prec_t <- b2_prec*2
    b4_prec_t <- b4_prec*2
    b0 ~ dbeta(b0.a,b0.b) 
    b3 ~ dbeta(b3.a,b3.b) 
    b4 ~ dbeta(b4.a,b4.b) 
    b2 ~ dbeta(b2.a,b2.b)

  }
    "

















