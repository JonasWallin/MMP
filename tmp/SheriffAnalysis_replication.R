### Replication of SheriffAnalysis.R but with functions


library(MMP)
data("sherifdat")


## Adjusted CEM (CEI)

##
# Setting up the model
# T ~ N(\mu0, \sigma_T)
# I ~ N(T, \exp(-delta * time)*\sigma_I)
# y ~ N(I, \sigma_y)
##

data <- getData(y.centered ~ time, 
                ~ 1 | person, 
                ~ 1 | group, 
                emergence = ~ -1 + time, 
                method = "CEI", 
                data = sherifdat)

object <- dataToObject(data)


param2 <- param0(object)

res2 <-optim(param2, function(x){-loglik(x, object) })
COVARIANCE_BETA2 <- getCovbeta(res2$par, object)

betas <-getbeta(res2$par, object)

paramList2 <- paramToList(res2$par, object)

m1 <- ce(y.centered ~ time, 
         ~ 1 | person, 
         ~ 1 | group, 
         emergence = ~ -1 + time, 
         method = "CEI", 
         data = sherifdat)

## results - all are the same!
betas == m1$betas
COVARIANCE_BETA2 == m1$cov_beta
paramList2[1]
m1$covariances[1]
paramList2[2]
m1$covariances[2]
paramList2[3]
m1$covariances[3]



## CEM model

##
# Setting up the model 2
# T ~ N(\mu0, \sigma_T)
# I ~ N(T,\sigma_I)
# y ~ N(I,  \exp(-delta * time)*\sigma_y)
##

data3 <- getData(y.centered ~ 1, 
                ~ 1 | person, 
                ~ 1 | group, 
                emergence = ~ 1 + time, 
                method = "CEM", 
                data = sherifdat)



object3 <- dataToObject(data3)

param3 <- param0(object3)

res3 <-optim(param3, function(x){-loglik(x, object3) })

## gives covariances for fixed effect betas
COVARIANCE_BETA3 <- getCovbeta(res3$par, object3)

## gives fixed effect betas
betas <-getbeta(res3$par, object3)

## covariances:
paramList3 <- paramToList(res3$par, object3)

# sigma squared:
exp(0.438)


# delta:
exp(-2.003/2) # (or /2? or some other alteration?)
-2.003/2 # prbably this, value/2

# var for indv intercept
nlme::pdLogChol(-4.19)

# cov for random stuff on group level
nlme::pdLogChol(c(0.00037,-1.77742,0.08971))


m2 <- ce(y.centered ~ 1, 
         ~ 1 | person, 
         ~ 1 | group, 
         emergence = ~ 1 + time, 
         method = "CEM", 
         data = sherifdat)

## results - all are the same!
betas == m2$betas
COVARIANCE_BETA3 == m2$cov_beta
paramList3[1]
m2$covariances[1]
paramList3[2]
m2$covariances[2]
paramList3[3]
m2$covariances[3]

## GP model

##
# Setting up the model 3
# T ~ N(\mu0, \sigma_T)
# I ~ T + GP(\theta)
# y ~ N(I, \sigma_y)
##

data4 <- getData(y.centered ~ 1, 
                 ~ 1 | person, 
                 ~ 1 | group, 
                 emergence = ~ 1, 
                 method = "GP",
                 time = "time",
                 data = sherifdat)

object4 <- dataToObject(data4)

param4 <- param0(object4)

res4 <-optim(param4, function(x){-loglik(x, object4) })
res4 <-optim(res4$par, function(x){-loglik(x, object4) }) # should this be run twice? yes
COVARIANCE_BETA4 <- getCovbeta(res4$par, object4)

paramList4 <- paramToList(res4$par, object4)
betas <-getbeta(res4$par, object4)

m3 <- ce(y.centered ~ 1, 
         ~ 1 | person, 
         ~ 1 | group, 
         emergence = ~ 1, 
         method = "GP",
         time = "time",
         data = sherifdat)


## results - none of them match. if running optim only once they match, which is what ce does
## so all are the same!
betas == m3$betas
COVARIANCE_BETA4 == m3$cov_beta 
paramList4[1]
m3$covariances[1]
paramList4[2]
m3$covariances[2]
paramList4[3]
m3$covariances[3]


## plots

###
#'
#' computing r for empirical and the three models
#'
###
if(save.fig){
  pdf('r_est.pdf')  
}else{
  x11()
}
r <- r.emperical(sherifdat$y.centered, sherifdat$group, sherifdat$time)
plot(r$time, r$r, ylim=c(0,1),cex=2, pch=19)
#
time <- seq(min(r$time), max(r$time), length.out = 100)

# CEM
r_2 <- exp(paramList3$error[[1]][1] + paramList3$error[[1]][2]*time) +exp(2* paramList3$indv[[1]])
r_2 <- sqrt(r_2/r_2[1])
lines(time, r_2, col='red')

# CEI
r_1 <- exp(-2*paramList2$indv[[1]][1]*time + 2*paramList2$indv[[1]][2]) + exp(paramList2$error[[1]])
r_1 <- sqrt(r_1/r_1[1])
lines(time, r_1, col='blue')

# GP
Cov <- rep(0, length(time))
for(i in 1:length(time)){
  indv <- list(D=as.matrix(1),time=time[i])
  Cov[i] <- object4$indvCovs[[1]]$get_Cov(paramList4$indv[[1]], indv)  + exp(paramList2$error[[1]]) # should the last really be this? Something is wrong here
}
r_3 <- sqrt(Cov/Cov[1])
lines(time, r_3, col='green',cex=2)

if(save.fig)
  dev.off()

betas <-getbeta(res$par, TeamObj)
Smooth <- smoothIndivual(res$par, TeamObj)
Smooth2 <- smoothIndivual(res2$par, TeamObjv2)
Smooth3 <- smoothIndivual(res3$par, TeamObjv3)


if(save.fig){
  pdf('smooth.pdf')  
}else{
  x11()
}
par(mfrow=c(2,3))
time <- c(0,1,2)
for(j in 1:2){
  for(i in 1:3){
    plot(time, as.vector(Matrix::t(TeamObj$teams[[j]]$indv[[i]]$A)%*%TeamObj$teams[[j]]$data$Y), ylab='obs',ylim=c(-1.5,1.5) ,cex=2, pch=19)
    m <- Smooth$teams[[j]]$indv[[i]]$mean
    s <-  sqrt(Smooth$teams[[j]]$indv[[i]]$var  + exp(paramList$error[[1]]))
    lines(time,m + 2*s, col='blue' )
    lines(time,m - 2*s, col='blue' )
    m <- Smooth2$teams[[j]]$indv[[i]]$mean
    s <-  sqrt(Smooth2$teams[[j]]$indv[[i]]$var  +  exp(paramList2$error[[1]][1] + paramList2$error[[1]][2]*time))
    lines(time,m + 2*s, col='red' )
    lines(time,m - 2*s, col='red' )
    
    m <- Smooth3$teams[[j]]$indv[[i]]$mean
    s <-  sqrt(Smooth3$teams[[j]]$indv[[i]]$var  +  exp(paramList3$error[[1]][1] ))
    lines(time,m + 2*s, col='green' )
    lines(time,m - 2*s, col='green' )
  }
}
if(save.fig)
  dev.off()


