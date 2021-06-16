## data generating functions for simulation

## Simulate type I consensus (CEM2 type)

gendat_CEM2 <- function(l3n,l2n,l1n,
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

gendat_CEM2GP <- function(l3n,l2n,l1n,
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
gendat_CEI2 <- function(l3n,l2n,l1n,
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

gendat_CEI2GP <- function(l3n,l2n,l1n,
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

