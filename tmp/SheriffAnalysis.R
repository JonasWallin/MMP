###
# Building the analysis of the classical Sheriff data
#
#
#
# D: 2020-03-16
##
library(MMP)
library(ggplot2)
library(dplyr)
save.fig=T
data("sherifdat")
graphics.off()

centered = T
##
# Setting up the model
# T ~ N(\mu0, \sigma_T)
# I ~ N(T, \exp(-delta * time)*\sigma_I)
# y ~ N(I, \sigma_y)
##
if(centered){
  y <-sherifdat$y.centered
}else{
  y <-sherifdat$y
}
TeamObj <- dataToObject(y = y,
                        team = sherifdat$group,
                        indv = sherifdat$person,
                        Xf   = rep(1, length(sherifdat$y)),
                        XI   = rep(1, length(sherifdat$y)),
                        XT   = rep(1,  length(sherifdat$y)),
                        WLI = F,
                        wEI  = sherifdat$time)
param <- param0(TeamObj)

res<-optim(param, function(x){-loglik(x, TeamObj) })
paramList <- paramToList(res$par, TeamObj)


##
# Setting up the model 2
# T ~ N(\mu0, \sigma_T)
# I ~ N(T,\sigma_I)
# y ~ N(I,  \exp(-delta * time)*\sigma_y)
##

TeamObjv2 <- dataToObject(y = y,
                        team = sherifdat$group,
                        indv = sherifdat$person,
                        Xf   = rep(1, length(sherifdat$y)),
                        XI   = rep(1, length(sherifdat$y)),
                        XT   = rep(1,  length(sherifdat$y)),
                        wNOISE = cbind(rep(1,  length(sherifdat$y)),sherifdat$time))
param <- param0(TeamObjv2)
res2<-optim(param, function(x){-loglik(x, TeamObjv2) })
paramList2 <- paramToList(res2$par, TeamObjv2)

##
# Setting up the model 3
# T ~ N(\mu0, \sigma_T)
# I ~ T + GP(\theta)
# y ~ N(I, \sigma_y)
##
TeamObjv3 <- dataToObject(y = y,
                        team = sherifdat$group,
                        indv = sherifdat$person,
                        Xf   = rep(1, length(sherifdat$y)),
                        XT   = rep(1,  length(sherifdat$y)),
                        TI = T,
                        time  = sherifdat$time)

param <- param0(TeamObjv3)
res3<-optim(param, function(x){-loglik(x, TeamObjv3) })
res3<-optim(res3$par, function(x){-loglik(x, TeamObjv3) })
paramList3 <- paramToList(res3$par, TeamObjv3)
Y <- y
pl <- ggplot(data=sherifdat, aes(x=time,y=Y, colour = person)) + geom_point( size=2) +
  xlab("time") + ylab("inches")
pl <- pl +  facet_wrap(~group)
x11()
print(pl)


###
#'
#' computing r for emperical and the three models
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
r_2 <- exp(paramList2$error[[1]][1] + paramList2$error[[1]][2]*time) +exp(2* paramList2$indv[[1]])
r_2 <- sqrt(r_2/r_2[1])
lines(time, r_2, col='red')
r_1 <- exp(-2*paramList$indv[[1]][1]*time + 2*paramList$indv[[1]][2]) + exp(paramList$error[[1]])
r_1 <- sqrt(r_1/r_1[1])
lines(time, r_1, col='blue')


Cov <- rep(0, length(time))
for(i in 1:length(time)){
  indv <- list(time=time[i])
  Cov[i] <- TeamObjv3$indvCovs[[1]]$get_Cov(paramList3$indv[[1]], indv) + + exp(paramList$error[[1]])
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
