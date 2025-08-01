###
# example using predicit for the sheriff data
#
#
###

library(MMP)
library(tidyverse)
## excluding individual measurement occasions
data("sherifdat")
sherifdat <- subset(sherifdat, time <= 2)
#sherifdat_lang <- subset(sherifdat, time <= 2 & time >=0)
sherifdat$time <- sherifdat$time + 1
#sherifdat$y <- sherifdat$y.centered



null <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 + time | group, 
           emergence = ~ 1, 
           method = "CEM2", 
           data = sherifdat,
           REML=T)

CEM2 <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 + time | group, 
           emergence = ~ 1 + time, 
           method = "CEM2", 
           data = sherifdat,
           REML=T)
CEI2 <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 + time | group, 
           emergence = ~ -1 + time, # to only extend with P_ij^0, still keep -1
           method = "CEI2", 
           data = sherifdat,
           REML = T)

GP <- ce(y ~ 1+time, 
         ~ 1 | person, 
         ~ 1 + time | group, 
         emergence = ~ 1, 
         method = "GP",
         time = "time",
         data = sherifdat,
         REML=T)
CEI2.h <- ce(y ~ 1+time, 
             ~ 1 | person, 
             ~ 1 | group, 
             emergence = ~ -1 + time, # 1 inkluderar "indivudal baseline variance"
             time = "time",
             method = "CEI2", 
             method.team = "OU.homeostasis",
             data = sherifdat,
             REML=T)
#suppose now I want to do preidiction of two indivuals I create a data.frame
#containing all relevant covariates
n.time = 100
sherifdat.new <- data.frame(person = c(rep(1,n.time), rep(2,n.time)),
                            group  = c(rep(1,n.time), rep(1,n.time)),
                            time   = c(seq(0,4.5,length.out = n.time),seq(0,4.5,length.out = n.time) ))
pred.GP.new   <- predict.ce(GP, sherifdat.new)
pred.CEI2.new <- predict.ce(CEI2, sherifdat.new)
CEI2.h.new <- predict.ce(CEI2.h, sherifdat.new)
index.1 <- sherifdat.new$person==1 & sherifdat.new$group==1
index.0 <- sherifdat$person==1 & sherifdat$group==1
plot(sherifdat.new$time[index.1], pred.GP.new[index.1,1],col='red',type='l',ylim=c(-3,9))
points(sherifdat$time[index.0],sherifdat$y[index.0])
lines(sherifdat.new$time[index.1],  pred.GP.new[index.1,1] + 2*sqrt(pred.GP.new[index.1,2]),col='red',lty=2)
lines(sherifdat.new$time[index.1],  pred.GP.new[index.1,1] - 2*sqrt(pred.GP.new[index.1,2]),col='red',lty=2)
lines(sherifdat.new$time[index.1], pred.CEI2.new[index.1,1],col='blue')
lines(sherifdat.new$time[index.1],  pred.CEI2.new[index.1,1] + 2*sqrt(pred.CEI2.new[index.1,2]),col='blue',lty=2)
lines(sherifdat.new$time[index.1],  pred.CEI2.new[index.1,1] - 2*sqrt(pred.CEI2.new[index.1,2]),col='blue',lty=2)
lines(sherifdat.new$time[index.1], CEI2.h.new[index.1,1],col='black')
lines(sherifdat.new$time[index.1],  CEI2.h.new[index.1,1] + 2*sqrt(CEI2.h.new[index.1,2]),col='black',lty=2)
lines(sherifdat.new$time[index.1],  CEI2.h.new[index.1,1] - 2*sqrt(CEI2.h.new[index.1,2]),col='black',lty=2)