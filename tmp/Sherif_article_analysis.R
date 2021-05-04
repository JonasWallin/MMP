### Sherif data analysis for article
### 2021-02-22

library(MMP)
library(tidyverse)
## excluding individual measurement occasions
data("sherifdat")
sherifdat <- subset(sherifdat, time <= 2)
#sherifdat_lang <- subset(sherifdat, time <= 2 & time >=0)
sherifdat$time <- sherifdat$time + 1
#sherifdat$y <- sherifdat$y.centered

# recreate CEM from Lang et al bookchapter
CEM <- ce(y ~ 1+time, 
          ~ 1 | person, 
          ~ 1 + time | group, 
          emergence = ~ 1 + time, 
          method = "CEM", 
          data = sherifdat)

summary.ce(CEM)
CEM$res$convergence
CEM$res

CEM.null <- ce(y ~ 1+time, 
          ~ 1 | person, 
          ~ 1 + time | group, 
          emergence = ~ 1, 
          method = "CEM", 
          data = sherifdat)

summary.ce(CEM.null)
CEM.null$res

CEM.bridge <- ce(y ~ 1+time, 
          ~ 1 | person, 
          ~ 1 | group, 
          emergence = ~ 1 + time, 
          time = "time",
          method = "CEM",
          method.team = "OU",
          data = sherifdat)

summary.ce(CEM.bridge)
CEM.bridge$res
CEM.homeostasis <- ce(y ~ 1+time, 
                 ~ 1 | person, 
                 ~ 1 | group, 
                 emergence = ~ 1 + time, 
                 time = "time",
                 method = "CEM",
                 method.team = "OU.homeostasis",
                 data = sherifdat)

summary.ce(CEM.homeostasis)


# bridge null model?
bridge.null <- ce(y ~ 1+time, 
                 ~ 1 | person, 
                 ~ 1 | group, 
                 emergence = ~ 1, 
                 time = "time",
                 method = "CEM",
                 method.team = "OU",
                 data = sherifdat)

summary.ce(bridge.null)
bridge.null$res$convergence

## Adjusted CEM
CEI.bridge <- ce(y ~ 1+time, 
          ~ 1 | person, 
          ~ 1 | group, 
          emergence = ~ -1 + time, # 1 inkluderar "indivudal baseline variance"
          time = "time",
          method = "CEI", 
          method.team = "OU",
          data = sherifdat)

summary.ce(CEI.bridge)
CEI.bridge$res$convergence

## Adjusted CEM
CEI <- ce(y ~ 1+time, 
          ~ 1 | person, 
          ~ 1 + time| group, 
          emergence = ~ -1 + time, # 1 inkluderar "indivudal baseline variance"
          method = "CEI", 
          data = sherifdat)

summary.ce(CEI)
CEI$res

## GP

GP <- ce(y ~ 1+time, 
         ~ 1 | person, 
         ~ 1 + time | group, 
         emergence = ~ 1, 
         method = "GP",
         time = "time",
         data = sherifdat)

summary.ce(GP)
GP$res

# GP bridge
GP.bridge <- ce(y ~ 1+time, 
         ~ 1 | person, 
         ~ 1 | group, 
         emergence = ~ 1, 
         method = "GP",
         method.team = "OU",
         time = "time",
         data = sherifdat)

summary.ce(GP.bridge)
GP.bridge$res

# r plot
r.plot(CEM, CEI, GP, sherifdat$y,sherifdat$group, sherifdat$time)
r.plot(CEM.bridge, CEI.bridge, GP.bridge, sherifdat$y,sherifdat$group, sherifdat$time)

smooth.plot(CEM, CEI, GP, sherifdat$y, sherifdat$time, groups.to.plot = c(1,8))
smooth.plot(CEM.bridge, CEI.bridge, GP.bridge, sherifdat$y, sherifdat$time, groups.to.plot = c(1,8))
