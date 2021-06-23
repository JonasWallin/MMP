## data generating functions for simulation

## Simulate type I consensus (CEM2 type)

gendat_HeCE <- function(l3n,l2n,l1n,
                        mu0,mu1, # fixed effects
                        sg00,sg11,sg01, # group var-cov
                        sp0, # individual prior variance
                        sp, # individual variance
                        se, # error variance
                        delta1){
  dat=expand.grid(time = 0:(l1n-1),
                  person = 1:l2n,
                  group = 1:l3n)
  g <- MASS::mvrnorm(l3n, c(mu0,mu1), matrix( c(sg00,sg01,sg01,sg11), 2) )
  dat$g1<-g[,1][dat[,3]]
  dat$g2<-g[,2][dat[,3]]
  dat$g <- dat$g1+dat$g2*dat$time
  dat$p0 <- rep(rnorm(l3n*l2n,0,sd=sqrt(sp0)), each=l1n)
  dat$p<- rnorm(l3n*l2n*l1n,dat$g+dat$p0,sd=sqrt(sp*exp(2*delta1*dat$time)))
  dat$y<- rnorm(l3n*l2n*l1n,dat$p,sd=sqrt(se))
  
  # return
  return(dat)
}

## Simulate type I consensus (CEM2 type)
## with group GP

gendat_HeCEGP <- function(l3n,l2n,l1n,
                          mu0,mu1, # fixed effects
                          sp0, # individual prior variance
                          sp, # individual variance
                          se, # error variance
                          sigma, # standard deviation GP
                          corr, # correlation largest distance between time points
                          delta1, # consensus
                          delta2 # group delta
){
  dat=expand.grid(time = 0:(l1n-1),
                  person = 1:l2n,
                  group = 1:l3n)
  
  time <- 0:(l1n-1)
  Dist <- as.matrix(dist(time))
  
  ## räkna ut kappa om vi vill ha correlation x mellan 
  ## längsta avståndet i data
  x <- 0.05 # correlation
  k <- max(Dist)/(-log(x))
  
  param <- c(sigma, -log(k)) # exp(param[2] ) = 1/k 
  Sigma <- OUcov(Dist, param) 
  
  # mean vector
  mean <- mu0 + mu1*time
  
  # var-cov
  covs <- diag(exp(delta2*time))%*%Sigma%*%diag(exp(delta2*time)) # sigma i mvrnorm
  
  g <- MASS::mvrnorm(l3n, mean, covs)
  gpdat <- data.frame(cbind("group"=rep(1:l3n,each=l1n),"time"=rep(time,l3n),"g"=as.vector(t(g))))
  dat <- dplyr::left_join(dat,gpdat,by=c("group","time"))
  
  dat$p0 <- rep(rnorm(l3n*l2n,0,sd=sqrt(sp0)), each=l1n)
  dat$p<- rnorm(l3n*l2n*l1n,dat$g+dat$p0,sd=sqrt(sp*exp(2*delta1*dat$time)))
  dat$y<- rnorm(l3n*l2n*l1n,dat$p,sd=sqrt(se))
  
  # return
  return(dat)
}

## generates ce type II (CEI2)
gendat_HoCE <- function(l3n,l2n,l1n,
                        mu0,mu1, # fixed effects
                        sg00,sg11,sg01, # group var-cov
                        sp0, # individual prior variance
                        sp, # individual variance
                        se, # error variance
                        delta1){
  dat=expand.grid(time = 0:(l1n-1),
                  person = 1:l2n,
                  group = 1:l3n)
  g <- MASS::mvrnorm(l3n, c(mu0,mu1), matrix( c(sg00,sg01,sg01,sg11), 2) )
  dat$g1<-g[,1][dat[,3]]
  dat$g2<-g[,2][dat[,3]]
  dat$g <- dat$g1+dat$g2*dat$time
  dat$p0 <- rep(rnorm(l3n*l2n,0,sd=sqrt(sp0)), each=l1n)
  dat$ep<- rep(rnorm(l3n*l2n,0,sd=sqrt(sp)),each=l1n)
  dat$p <- dat$g+dat$p0+dat$ep*exp(delta1*dat$time)
  dat$y<- rnorm(l3n*l2n*l1n,dat$p,sd=sqrt(se))
  
  # return
  return(dat)
}

## Simulate type II consensus (CEI2 type)
## with group GP

gendat_HoCEGP <- function(l3n,l2n,l1n,
                          mu0,mu1, # fixed effects
                          sp0, # individual prior variance
                          sp, # individual variance
                          se, # error variance
                          sigma, # standard deviation GP
                          corr, # correlation largest distance between time points
                          delta1, # consensus
                          delta2 # group delta
){
  dat=expand.grid(time = 0:(l1n-1),
                  person = 1:l2n,
                  group = 1:l3n)
  
  time <- 0:(l1n-1)
  Dist <- as.matrix(dist(time))
  
  ## räkna ut kappa om vi vill ha correlation x mellan 
  ## längsta avståndet i data
  x <- 0.05 # correlation
  k <- max(Dist)/(-log(x))
  
  param <- c(sigma, -log(k)) # exp(param[2] ) = 1/k 
  Sigma <- OUcov(Dist, param) 
  
  # mean vector
  mean <- mu0 + mu1*time
  
  # var-cov
  covs <- diag(exp(delta2*time))%*%Sigma%*%diag(exp(delta2*time)) # sigma i mvrnorm
  
  g <- MASS::mvrnorm(l3n, mean, covs)
  gpdat <- data.frame(cbind("group"=rep(1:l3n,each=l1n),"time"=rep(time,l3n),"g"=as.vector(t(g))))
  dat <- dplyr::left_join(dat,gpdat,by=c("group","time"))
  
  dat$p0 <- rep(rnorm(l3n*l2n,0,sd=sqrt(sp0)), each=l1n)
  dat$ep<- rep(rnorm(l3n*l2n,0,sd=sqrt(sp)),each=l1n)
  dat$p <- dat$g+dat$p0+dat$ep*exp(delta1*dat$time)
  dat$y<- rnorm(l3n*l2n*l1n,dat$p,sd=sqrt(se))
  
  # return
  return(dat)
}





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
  
  # Lang et al CEM
  CEM <-try(nlme::lme(y ~ time, random = list(group=pdSymm(~time),
                                             person=pdIdent(~1)),data=data,
                     weights=varExp( form = ~ time),
                     control=lmeControl(maxIter=15000,msMaxIter=15000)))
  
  
  
  # what output do we want to save? 
  # delta1, sp, convergence, seed
  # N = number of reps = 1000
  
  time <- unique(data$time)
  
  if(!inherits(HeCEM, 'try-error')){
    #d1 <- unlist(HeCEM$covariances$indv[[2]])[2]
    HeCEM_param <- c(d1 = unlist(HeCEM$covariances$indv[[2]])[2],
                     sp =exp(unlist(HeCEM$covariances$indv[[2]])[1]),
                     r_tmax = rt(HeCEM, time)[max(time+1),2])
  } else {
    #d1 <- NA
    HeCEM_param <- c(d1 = NA, sp = NA, r_tmax = NA)
  }
  
  if(!inherits(HeCEMGP, 'try-error')){
    #d2 <- unlist(HeCEMGP$covariances$indv[[2]])[2]
    HeCEMGP_param <- c(d1 = unlist(HeCEMGP$covariances$indv[[2]])[2],
                     sp =exp(unlist(HeCEMGP$covariances$indv[[2]])[1]),
                     r_tmax = rt(HeCEMGP, time)[max(time+1),2])
  } else {
    #d2 <- NA
    HeCEMGP_param <- c(d1 = NA, sp = NA, r_tmax = NA)
  }
  
  if(!inherits(HoCEM, 'try-error')){
    #d3 <- unlist(HoCEM$covariances$indv[[2]])[1]
    HoCEM_param <- c(d1 = unlist(HoCEM$covariances$indv[[2]])[1],
                      sp =exp(2*unlist(HoCEM$covariances$indv[[2]])[2]),
                     r_tmax = rt(HoCEM, time)[max(time+1),2])
  } else {
    #d3 <- NA
    HoCEM_param <- c(d1 = NA, sp = NA, r_tmax = NA)
  }
  
  if(!inherits(HoCEMGP, 'try-error')){
    #d4 <- unlist(HoCEMGP$covariances$indv[[2]])[1]
    HoCEMGP_param <- c(d1 = unlist(HoCEMGP$covariances$indv[[2]])[1],
                       sp =exp(2*unlist(HoCEMGP$covariances$indv[[2]])[2]),
                       r_tmax = rt(HoCEMGP, time)[max(time+1),2])
  } else {
    #d4 <- NA
    HoCEMGP_paramc(d1 = NA, sp = NA, r_tmax = NA)
  }
  
  if(!inherits(gp, 'try-error')){
    #d5 <- unlist(gp$covariances$indv[[2]])[1]
    GP_param <- c(d1 = unlist(gp$covariances$indv[[2]])[1],
                  sp =((exp(unlist(gp$covariances$indv[[2]])[2]))^2)/(2*exp(unlist(gp$covariances$indv[[2]])[3])),
                  r_tmax = rt(gp, time)[max(time+1),2])
  } else {
    #d5 <- NA
    GP_param <- c(d1 = NA, sp = NA, r_tmax = NA)
  }
  
  if(!inherits(gpGP, 'try-error')){
    #d6 <- unlist(gpGP$covariances$indv[[2]])[1]
    gpGP_param <- c(d1 = unlist(gpGP$covariances$indv[[2]])[1],
                  sp =((exp(unlist(gp$covariances$indv[[2]])[2]))^2)/(2*exp(unlist(gp$covariances$indv[[2]])[3])),
                  r_tmax = rt(gpGP, time)[max(time+1),2])
  } else {
    #d6 <- NA
    gpGP_param  <- c(d1 = NA, sp = NA, r_tmax = NA)
  }
  
  if((!inherits(CEM, 'try-error'))) {
    CEM_param <- c(d1 =CEM$modelStruct$varStruct,
                   se=as.numeric(VarCorr(CEM)[6,1]),
                   r_tmax = NA)
  } else {
    CEM_param <- c(d1 = NA, se = NA, r_tmax = NA)
  }
  
  out <- list(HeCEM=HeCEM_param,
              HoCEM=HoCEM_param,
              GP=GP_param,
              HeCEMGP=HeCEMGP_param,
              HoCEMGP=HoCEMGP_param,
              GPGP=gpGP_param,
              CEM=CEM_param)
  
  ## return delta1 for each model
  return(out)
}


## OLD
## Simulate type I consensus (CEM type)

gendat_I <- function(l3n,l2n,l1n,
                     mu0,mu1, # fixed effects
                     sg00,sg11,sg01, # group var-cov
                     sp, # individual variance
                     se, # error variance
                     delta1){
  dat=expand.grid(time = 0:(l1n-1),
                  person = 1:l2n,
                  group = 1:l3n)
  g <- MASS::mvrnorm(l3n, c(mu0,mu1), matrix( c(sg00,sg01,sg01,sg11), 2) )
  dat$g1<-g[,1][dat[,3]]
  dat$g2<-g[,2][dat[,3]]
  dat$g <- dat$g1+dat$g2*dat$time
  dat$p0 <- rep(rnorm(l3n*l2n,0,sd=sqrt(sp)), each=l1n)
  dat$p<- rnorm(l3n*l2n*l1n,dat$p0+dat$g,sd=sqrt(se*exp(2*delta1*dat$time)))
  
  # return
  return(dat)
}


## Simulate type I consensus (CEM type)
## with group GP

gendat_IGP <- function(l3n,l2n,l1n,
                       mu0,mu1, # fixed effects
                       sp, # individual variance
                       se, # error variance
                       sigma, # standard deviation GP
                       corr, # correlation largest distance between time points
                       delta1, # consensus
                       delta2 # group delta
){
  dat=expand.grid(time = 0:(l1n-1),
                  person = 1:l2n,
                  group = 1:l3n)
  
  time <- 0:(l1n-1)
  Dist <- as.matrix(dist(time))
  
  ## räkna ut kappa om vi vill ha correlation x mellan 
  ## längsta avståndet i data
  x <- 0.05 # correlation
  k <- max(Dist)/(-log(x))
  
  param <- c(sigma, -log(k)) # exp(param[2] ) = 1/k 
  Sigma <- OUcov(Dist, param) 
  
  # mean vector
  mean <- mu0 + mu1*time
  
  # var-cov
  covs <- diag(exp(delta2*time))%*%Sigma%*%diag(exp(delta2*time)) # sigma i mvrnorm
  
  g <- MASS::mvrnorm(l3n, mean, covs)
  gpdat <- data.frame(cbind("group"=rep(1:l3n,each=l1n),"time"=rep(time,l3n),"g"=as.vector(t(g))))
  dat <- dplyr::left_join(dat,gpdat,by=c("group","time"))
  
  dat$p0 <- rep(rnorm(l3n*l2n,0,sd=sqrt(sp)), each=l1n)
  dat$p<- rnorm(l3n*l2n*l1n,dat$p0+dat$g,sd=sqrt(se*exp(2*delta1*dat$time)))
  
  # return
  return(dat)
}


## generates ce type II (CEI)
gendat_II <- function(l3n,l2n,l1n,
                      mu0,mu1, # fixed effects
                      sg00,sg11,sg01, # group var-cov
                      sp, # individual variance
                      se, # error variance
                      delta1){
  dat=expand.grid(time = 0:(l1n-1),
                  person = 1:l2n,
                  group = 1:l3n)
  g <- MASS::mvrnorm(l3n, c(mu0,mu1), matrix( c(sg00,sg01,sg01,sg11), 2) )
  dat$g1<-g[,1][dat[,3]]
  dat$g2<-g[,2][dat[,3]]
  dat$g <- dat$g1+dat$g2*dat$time
  dat$ep<- rep(rnorm(l3n*l2n,0,sd=sqrt(sp)),each=l1n)
  dat$p <- dat$g+dat$ep*exp(delta1*dat$time)
  dat$y<- rnorm(l3n*l2n*l1n,dat$p,sd=sqrt(se))
  
  # return
  return(dat)
}

## Simulate type II consensus (CEI type)
## with group GP

gendat_IIGP <- function(l3n,l2n,l1n,
                        mu0,mu1, # fixed effects
                        sp, # individual variance
                        se, # error variance
                        sigma, # standard deviation GP
                        corr, # correlation largest distance between time points
                        delta1, # consensus
                        delta2 # group delta
){
  dat=expand.grid(time = 0:(l1n-1),
                  person = 1:l2n,
                  group = 1:l3n)
  
  time <- 0:(l1n-1)
  Dist <- as.matrix(dist(time))
  
  ## räkna ut kappa om vi vill ha correlation x mellan 
  ## längsta avståndet i data
  x <- 0.05 # correlation
  k <- max(Dist)/(-log(x))
  
  param <- c(sigma, -log(k)) # exp(param[2] ) = 1/k 
  Sigma <- OUcov(Dist, param) 
  
  # mean vector
  mean <- mu0 + mu1*time
  
  # var-cov
  covs <- diag(exp(delta2*time))%*%Sigma%*%diag(exp(delta2*time)) # sigma i mvrnorm
  
  g <- MASS::mvrnorm(l3n, mean, covs)
  gpdat <- data.frame(cbind("group"=rep(1:l3n,each=l1n),"time"=rep(time,l3n),"g"=as.vector(t(g))))
  dat <- dplyr::left_join(dat,gpdat,by=c("group","time"))
  
  dat$ep<- rep(rnorm(l3n*l2n,0,sd=sqrt(sp)),each=l1n)
  dat$p <- dat$g+dat$ep*exp(delta1*dat$time)
  dat$y<- rnorm(l3n*l2n*l1n,dat$p,sd=sqrt(se))
  
  # return
  return(dat)
}

