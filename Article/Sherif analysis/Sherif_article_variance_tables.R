
library(MMP)
library(tidyverse)
library(ggplot2)

## Data preparation
data("sherifdat")
# excluding last time point (individual measurement)
sherifdat <- subset(sherifdat, time <= 2)
# recode time to start from 0
sherifdat$time <- sherifdat$time + 1



# HomCEM
CEI2.h <- ce(y ~ 1+time, 
             ~ 1 | person, 
             ~ 1 | group, 
             emergence = ~ -1 + time, 
             time = "time",
             method = "CEI2", 
             method.team = "OU.homeostasis",
             data = sherifdat)
Covs <- get.Cov(CEI2.h$covariances,CEI2.h$object)

var.Err <- round(diag(Covs$SigmaE),2)
var.Indv <- round(diag(Covs$SigmaI),2)
var.Team <- round(diag(Covs$SigmaT),2)
var.Tot <- round(diag(Covs$Sigma),2)


Var.abs <- cbind(var.Err,var.Indv,var.Team,var.Tot)

var.Err.rel <- round(diag(Covs$SigmaE)/diag(Covs$Sigma),2)
var.Indv.rel <- round(diag(Covs$SigmaI)/diag(Covs$Sigma),2)
var.Team.rel <- round(diag(Covs$SigmaT)/diag(Covs$Sigma),2)
cat('V[I_t]/V[I_0] = ',var.Err.rel,'\n')
cat('V[I_t]/V[Y_t] = ',var.Indv.rel,'\n')
cat('V[G_t]/V[Y_t] = ',var.Team.rel,'\n')
rel <- cbind(var.Err.rel,var.Indv.rel,var.Team.rel)

