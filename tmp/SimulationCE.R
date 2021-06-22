### Simulation
# 2021-06-22

library(parallel)


## function for generating a data set and run all models
ce_sim <- function(l3n,l2n,l1n,
                   mu0,
                   mu1, # fixed effects
                   sg00=NULL, # group intercept variance
                   sg11=NULL, # group slope variance
                   sg01=NULL, # group intercept/slope cov
                   sp0, # individual prior variance
                   sp, # individual variance
                   se, # error variance
                   sigma=NULL, # standard deviation GP
                   corr=NULL, # correlation largest distance between time points
                   delta1,
                   delta2=NULL,
                   datatype) # "HeCEM","HoCEM","HeCEM+GP","HoCEM+GP" 
{
  
  ## Generate data         
  if (datatype == "HeCEM") {
    ## Consensus emergence type I (CEM2 type)
    data <- gendat_HeCE(l3n,l2n,l1n,
                        mu0,mu1,
                        sg00,sg11,sg01, 
                        sp0,
                        sp,
                        se,
                        delta1)
  } 
  if (datatype == "HoCEM") {
    ## Consensus emergence type II (CEI2 type)
    data <- gendat_HoCE(l3n,l2n,l1n,
                        mu0,mu1,
                        sg00,sg11,sg01, 
                        sp0,
                        sp,
                        se,
                        delta1)
  } 
  
  if (datatype == "HeCEM+GP") {
    ## Consensus emergence type I (CEM2 type) + group process
    data <- gendat_HeCEGP(l3n,l2n,l1n,
                          mu0,mu1, 
                          sp0,
                          sp,
                          se,
                          sigma,
                          corr,
                          delta1,
                          delta2)
  } 
  if (datatype == "HoCEM+GP") {
    ## Consensus emergence type II (CEI2 type) + group process
    data <- gendat_HoCEGP(l3n,l2n,l1n,
                          mu0,mu1, 
                          sp0,
                          sp,
                          se,
                          sigma,
                          corr,
                          delta1,
                          delta2)
  } 
  
  
  ## Estimate models
  ## HeCEM
  # HeCEM regular
  HeCEM <- try(ce(y ~ 1+time, 
                  ~ 1 | person, 
                  ~ 1 + time | group, 
                  emergence = ~ 1 + time, 
                  method = "CEM2", 
                  data = data))
  
  
  # HeCEM + group process
  HeCEMGP <- try(ce(y ~ 1+time, 
                    ~ 1 | person, 
                    ~ 1 | group, 
                    emergence = ~ 1 + time, 
                    method = "CEM2", 
                    time = "time",
                    method.team = "OU.homeostasis",
                    data = data))
  
  
  ## HoCEM
  # HoCEM regular
  HoCEM <- try(ce(y ~ 1+time, 
                  ~ 1 | person, 
                  ~ 1 + time | group, 
                  emergence = ~ -1 + time, 
                  method = "CEI2", 
                  data = data))
  
  
  # HoCEM + group process
  HoCEMGP <- try(ce(y ~ 1+time, 
                    ~ 1 | person, 
                    ~ 1 | group, 
                    emergence = ~ -1 + time, 
                    method = "CEI2", 
                    time = "time",
                    method.team = "OU.homeostasis",
                    data = data))
  
  ## GP
  # GP regular
  gp <- try(ce(y ~ 1+time, 
               ~ 1 | person, 
               ~ 1 + time | group, 
               emergence = ~ 1, 
               method = "GP",
               time = "time",
               data = data))
  
  # GP  + group process
  gpGP <- try(ce(y ~ 1+time, 
                 ~ 1 | person, 
                 ~ 1 | group, 
                 emergence = ~ 1, 
                 method = "GP",
                 time = "time",
                 method.team = "OU.homeostasis",
                 data = data))
  
  # what output do we want to save? 
  # delta1, convergence, seed
  # N = number of reps = 1000
  
  
  if(!inherits(HeCEM, 'try-error')){
    d1 <- unlist(HeCEM$covariances$indv[[2]])[2]
  } else {
    d1 <- NA
  }
  
  if(!inherits(HeCEMGP, 'try-error')){
    d2 <- unlist(HeCEMGP$covariances$indv[[2]])[2]
  } else {
    d2 <- NA
  }
  
  if(!inherits(HoCEM, 'try-error')){
    d3 <- unlist(HoCEM$covariances$indv[[2]])[1]
  } else {
    d3 <- NA
  }
  
  if(!inherits(HoCEMGP, 'try-error')){
    d4 <- unlist(HoCEMGP$covariances$indv[[2]])[1]
  } else {
    d4 <- NA
  }
  
  if(!inherits(gp, 'try-error')){
    d5 <- unlist(gp$covariances$indv[[2]])[1]
  } else {
    d5 <- NA
  }
  
  if(!inherits(gpGP, 'try-error')){
    d6 <- unlist(gpGP$covariances$indv[[2]])[1]
  } else {
    d6 <- NA
  }
  
  out <- c(HeCEM=d1,HoCEM=d3,GP=d5,HeCEMGP=d2,HoCEMGP=d4,GPGP=d6)
  
  ## return delta1 for each model
  return(out)
}



# # of reps
N <- 30

system.time(res <- sapply(1:N, function(i,...) {
  set.seed(123+i); 
  c(ce_sim(l3n=3,l2n=3,l1n=3,
           mu0=3,mu1=0.01,
           sg00=0.05,sg11=0.005,sg01=-0.005,
           sp0=0.01,
           sp=0.2,
           se=0.01,
           delta1=-0.8,
           datatype = "HeCEM"),
    seed =123+i)}  ))

simresults <- cbind(delta1 = rowMeans(res[1:6,],na.rm=T), NAs=rowMeans(is.na(res[1:6,])))


# setup
no_cores <- detectCores() - 1 
cl<-makeCluster(no_cores)
clusterExport(cl,ls())
clusterEvalQ(cl, library("MMP"))


system.time(res <- parSapply(cl,1:N, function(i,...) {
    set.seed(123+i); 
    c(ce_sim(l3n=3,l2n=3,l1n=3,
             mu0=3,mu1=0.01,
             sg00=0.05,sg11=0.005,sg01=-0.005,
             sp0=0.01,
             sp=0.2,
             se=0.01,
             #sigma = 0.01,
             #corr = 0.05,
             delta1=-0.8,
             #delta2 = 0.25,
             datatype = "HeCEM"),
      seed =123+i)}  ))

simresults <- cbind(delta1 = rowMeans(res[1:6,],na.rm=T), NAs=rowMeans(is.na(res[1:6,])))

stopCluster(cl)

