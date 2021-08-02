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
                          sgp, # GP variance at time 0
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
  k <- max(Dist)/(-log(corr))
  
  param <- c(log(sqrt((2 * 1/k) *sgp)), -log(k)) # exp(param[2] ) = 1/k 
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
                          sgp, # variance of gp at time 0
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
  k <- max(Dist)/(-log(corr))
  
  param <- c(log(sqrt((2 * 1/k) *sgp)), -log(k)) # exp(param[2] ) = 1/k 
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
                   sgp=NULL, # variance of GP at time 0
                   corr=NULL, # correlation largest distance between time points
                   delta1,
                   delta2=NULL,
                   datatype, # "HeCEM","HoCEM","HeCEM+GP","HoCEM+GP"
                   REML = F,
                   index,
                   path=NULL)  
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
                          sgp,
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
                          sgp,
                          corr,
                          delta1,
                          delta2)
  } 
  
  
  ## Estimate models
  ## HeCEM
  # HeCEM regular
  HeCEM <- tryCatch(ce(y ~ 1+time,
                       ~ 1 | person, 
                       ~ 1 + time | group, 
                       emergence = ~ 1 + time, 
                       method = "CEM2", 
                       data = data,
                       REML = REML),
                    error = function(e){NA})
  
  
  # HeCEM + group process
  HeCEMGP <- tryCatch(ce(y ~ 1+time,
                         ~ 1 | person, 
                         ~ 1 | group, 
                         emergence = ~ 1 + time, 
                         method = "CEM2", 
                         time = "time",
                         method.team = "OU.homeostasis",
                         data = data,
                         REML = REML), 
                      error = function(e){NA})
  
  
  ## HoCEM
  # HoCEM regular
  HoCEM <- tryCatch(ce(y ~ 1+time,
                       ~ 1 | person,
                       ~ 1 + time | group,
                       emergence = ~ -1 + time,
                       method = "CEI2",
                       data = data,
                       REML = REML),
                    error = function(e){NA})
  
  
  # HoCEM + group process
  HoCEMGP <- tryCatch(ce(y ~ 1+time,
                         ~ 1 | person, 
                         ~ 1 | group, 
                         emergence = ~ -1 + time, 
                         method = "CEI2", 
                         time = "time",
                         method.team = "OU.homeostasis",
                         data = data,
                         REML = REML),
                      error = function(e){NA})
  
  ## GP
  # GP regular
  gp <- tryCatch(ce(y ~ 1+time,
                    ~ 1 | person, 
                    ~ 1 + time | group, 
                    emergence = ~ 1, 
                    method = "GP",
                    time = "time",
                    data = data,
                    REML = REML),
                 error = function(e){NA})
  
  # GP  + group process
  gpGP <- tryCatch(ce(y ~ 1+time,
                      ~ 1 | person,
                      ~ 1 | group, 
                      emergence = ~ 1, 
                      method = "GP",
                      time = "time",
                      method.team = "OU.homeostasis",
                      data = data,
                      REML = REML),
                   error = function(e){NA})
  
  # Lang et al CEM
  CEM <-tryCatch(nlme::lme(y ~ time, random = list(group=pdSymm(~time),
                                             person=pdIdent(~1)),data=data,
                     weights=varExp( form = ~ time),
                     control=lmeControl(maxIter=15000,msMaxIter=15000)),
                 error = function(e){NA})
  
  ## results
  
  time <- unique(data$time)
  
  HeCEM_param <- tryCatch(c(d1 = unlist(HeCEM$covariances$indv[[2]])[2]/2,
                            sp =exp(unlist(HeCEM$covariances$indv[[2]])[1]),
                            r_tmax = rt(HeCEM, time)[max(time+1),2],
                            loglik = HeCEM$loglik),
                          error = function(e){
                            c(d1 = NA, sp = NA, r_tmax = NA, loglik=NA)})
  
  HeCEMGP_param <- tryCatch(c(d1 = unlist(HeCEMGP$covariances$indv[[2]])[2]/2,
                              sp =exp(unlist(HeCEMGP$covariances$indv[[2]])[1]),
                              r_tmax = rt(HeCEMGP, time)[max(time+1),2],
                              loglik = HeCEMGP$loglik),
                            error = function(e){
                              c(d1 = NA, sp = NA, r_tmax = NA, loglik=NA)})
  

  HoCEM_param <- tryCatch(c(d1 = unlist(HoCEM$covariances$indv[[2]])[1],
                            sp =exp(2*unlist(HoCEM$covariances$indv[[2]])[2]),
                            r_tmax = rt(HoCEM, time)[max(time+1),2],
                            loglik = HoCEM$loglik),
                          error = function(e){
                            c(d1 = NA, sp = NA, r_tmax = NA, loglik=NA)}) 
  

  HoCEMGP_param <- tryCatch(c(d1 = unlist(HoCEMGP$covariances$indv[[2]])[1],
                              sp =exp(2*unlist(HoCEMGP$covariances$indv[[2]])[2]),
                              r_tmax = rt(HoCEMGP, time)[max(time+1),2],
                              loglik = HoCEMGP$loglik),
                            error = function(e){
                              c(d1 = NA, sp = NA, r_tmax = NA, loglik=NA)})

  gp_param <- tryCatch(c(d1 = unlist(gp$covariances$indv[[2]])[1],
                         sp =((exp(unlist(gp$covariances$indv[[2]])[2]))^2)/(2*exp(unlist(gp$covariances$indv[[2]])[3])),
                         r_tmax = rt(gp, time)[max(time+1),2],
                         loglik = gp$loglik),
                       error = function(e){
                         c(d1 = NA, sp = NA, r_tmax = NA, loglik=NA)})

  gpGP_param <- tryCatch(c(d1 = unlist(gpGP$covariances$indv[[2]])[1],
                           sp =((exp(unlist(gpGP$covariances$indv[[2]])[2]))^2)/(2*exp(unlist(gpGP$covariances$indv[[2]])[3])),
                           r_tmax = rt(gpGP, time)[max(time+1),2],
                           loglik = gpGP$loglik),
                         error = function(e){
                           c(d1 = NA, sp = NA, r_tmax = NA, loglik=NA)})

  CEM_param <- tryCatch(c(d1 =CEM$modelStruct$varStruct,
                          se=as.numeric(VarCorr(CEM)[6,1]),
                          r_tmax = sqrt((as.numeric(VarCorr(CEM)[5,1]) + as.numeric(VarCorr(CEM)[6,1])*exp(2*CEM$modelStruct$varStruct*time))/
                                          (as.numeric(VarCorr(CEM)[5,1]) + as.numeric(VarCorr(CEM)[6,1])*exp(2*CEM$modelStruct$varStruct*time))[1])[max(time+1)],
                          loglik= as.numeric(CEM$logLik)),
                        error = function(e){
                          c(d1 = NA, sp = NA, r_tmax = NA, loglik=NA)})

  
  out <- list(HeCEM=HeCEM_param,
              HoCEM=HoCEM_param,
              GP=gp_param,
              HeCEMGP=HeCEMGP_param,
              HoCEMGP=HoCEMGP_param,
              GPGP=gpGP_param,
              CEM=CEM_param)
  
  ## save stuff
  # want to save seed, index, data, models and output
  all_output <- list(seed_info = c(seed = 123+index, index=index),
                     data = data, 
                     fitted_models = list(HeCEM=HeCEM,
                                          HoCEM=HoCEM,
                                          GP=gp,
                                          HeCEMGP=HeCEMGP,
                                          HoCEMGP=HoCEMGP,
                                          GPGP=gpGP,
                                          CEM=CEM),
                     output = out)
  
  if (!is.null(path)){
    saveRDS(all_output,file =paste0(path,  delta1,  "_",index, ".Rdata"))
  }
  
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

