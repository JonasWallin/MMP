### Replication of SheriffAnalysis.R but with functions


library(MMP)
data("sherifdat")


## Adjusted CEM

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

paramList2 <- paramToList(res2$par, object)

## gets the same results. BUT! needs to have WLI = F, WLT = T

## CEM model

data3 <- getData(y.centered ~ 1, 
                ~ 1 | person, 
                ~ 1 | group, 
                emergence = ~ 1 + time, 
                method = "CEM", 
                data = sherifdat)


object3 <- dataToObject(data3)

param3 <- param0(object3)

res3 <-optim(param3, function(x){-loglik(x, object3) })
COVARIANCE_BETA3 <- getCovbeta(res3$par, object3)

paramList3 <- paramToList(res3$par, object3)

## Works but needs WLI = T, WLT = T.


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
res4 <-optim(res4$par, function(x){-loglik(x, object4) }) # should this be run twice?
COVARIANCE_BETA4 <- getCovbeta(res4$par, object4)

paramList4 <- paramToList(res4$par, object4)

# Works but needs to have WLI = F, WLT = T