library(MMP)
library(ggplot2)
library(dplyr)

data("sherifdat")
sherifdat <- subset(sherifdat, time >= 0 & time <= 2)

## CEM model

##
# Setting up the model 2
# T ~ N(\mu0, \sigma_T)
# I ~ N(T,\sigma_I)
# y ~ N(I,  \exp(-delta * time)*\sigma_y)
##

CEM <- ce(y ~ 1, 
         ~ 1 | person, 
         ~ 1 | group, 
         emergence = ~ 1 + time, 
         method = "CEM", 
         data = sherifdat)



## Adjusted CEM (CEI)

##
# Setting up the model
# T ~ N(\mu0, \sigma_T)
# I ~ N(T, \exp(-delta * time)*\sigma_I)
# y ~ N(I, \sigma_y)
##

CEI <- ce(y ~ 1, 
         ~ 1 | person, 
         ~ 1 | group, 
         emergence = ~ -1 + time, 
         method = "CEI", 
         data = sherifdat)

summary.ce(CEI)

## GP model

##
# Setting up the model 3
# T ~ N(\mu0, \sigma_T)
# I ~ T + GP(\theta)
# y ~ N(I, \sigma_y)
##


GP <- ce(y ~ 1, 
         ~ -1 | person, 
         ~ 1 | group, 
         emergence = ~ 1, 
         method = "GP",
         time = "time",
         data = sherifdat)

summary.ce(GP)

r.plot(CEM, CEI, GP, sherifdat$y,sherifdat$group, sherifdat$time)

smooth.plot(CEM, CEI, GP, sherifdat$y, sherifdat$time)

