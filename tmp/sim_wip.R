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

## generates ce type II
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

set.seed(123)
# data type I ce
## same values as used in Lang et al 2019 except for sp and se,
## changed so total variance should be the same (?)
data_I <- gendat_I(l3n=5,l2n=5,l1n=10,
                   mu0=3,mu1=0.01,
                   sg00=0.05,sg11=0.005,sg01=-0.005, 
                   sp=0.1,
                   se=0.2,
                   delta1=-0.25)

# data type II ce
data_II <- gendat_II(l3n=5,l2n=5,l1n=10,
                   mu0=3,mu1=0.01,
                   sg00=0.05,sg11=0.005,sg01=-0.005,
                   sp=0.2,
                   se=0.1,
                   delta1=-0.25)

# data type I ce + GP
data_IGP <- gendat_IGP(l3n=5,l2n=5,l1n=4,
                     mu0=3,mu1=0.01,
                     sp=0.01,
                     se=0.2,
                     sigma=0.1, # standard deviation GP
                     corr=0.05, # correlation largest distance between time points
                     delta1=-0.25,
                     delta2=-0.15)

# data type II ce + GP
data_IIGP <- gendat_IIGP(l3n=5,l2n=5,l1n=4,
                     mu0=3,mu1=0.01,
                     sp=0.2,
                     se=0.01,
                     sigma=0.1, # standard deviation GP
                     corr=0.05, # correlation largest distance between time points
                     delta1=-0.25,
                     delta2=-0.15)

## ce type I, cei model
cei_I <- ce(p ~ 1+time, 
          ~ 1 | person, 
          ~ 1 + time | group, 
          emergence = ~ -1 + time, 
          method = "CEI", 
          data = data_I)

summary.ce(cei_I)

## ce type I + GP, cem model 
cem_IGP <- ce(p ~ 1+time, 
          ~ 1 | person, 
          ~ 1 + time | group, 
          emergence = ~ 1 + time, 
          method = "CEM", 
          data = data_IGP)

# error in solve.default(t(L), y_i) : 
# system is computationally singular: reciprocal condition number = 2.48091e-23 
summary.ce(cem_IGP)

# same analysis with nlme (works fine)
library(nlme)

cem_lme_I <- lme(p~time, 
    random = list(group=pdLogChol(~time),
                  person=pdIdent(~1)),
    data=data_IGP, weights=varExp(form =~time))

summary(cem_lme_I)

## ce type II, cei model
cei_IIGP <- ce(y ~ 1+time, 
            ~ 1 | person, 
            ~ 1 + time | group, 
            emergence = ~ -1 + time, 
            method = "CEI", 
            data = data_IIGP)
summary.ce(cei_IIGP)

cei_IIGPh <- ce(y ~ 1+time, 
               ~ 1 | person, 
               ~ 1| group, 
               emergence = ~ -1 + time, 
               method = "CEI", 
               method.team = "OU.homeostasis",
               data = data_IIGP)
summary.ce(cei_IIGP)


cem_II <- ce(y ~ 1+time, 
            ~ 1 | person, 
            ~ 1 + time | group, 
            emergence = ~ 1 + time, 
            method = "CEM", 
            data = data_II)
summary.ce(cem_II)

cem_lme_II <- lme(y~time, 
                 random = list(group=pdLogChol(~time),
                               person=pdIdent(~1)),
                 data=data_II,
                 weights=varExp(form =~time))

summary(cem_lme_II)
## gives similar results