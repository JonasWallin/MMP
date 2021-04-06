### Sherif data analysis for article
### 2021-02-22

library(MMP)
library(tidyverse)
## excluding individual measurement occasions
data("sherifdat")
sherifdat <- subset(sherifdat, time <= 2)
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

CEM.null <- ce(y ~ 1+time, 
          ~ 1 | person, 
          ~ 1 + time | group, 
          emergence = ~ 1, 
          method = "CEM", 
          data = sherifdat)

summary.ce(CEM.null)

CEM.bridge <- ce(y ~ 1+time, 
          ~ 1 | person, 
          ~ 1 | group, 
          emergence = ~ 1 + time, 
          time = "time",
          method = "CEM",
          method.team = "OU",
          data = sherifdat)

summary.ce(CEM.bridge)

## Adjusted CEM
CEI.bridge <- ce(y ~ 1+time, 
          ~ 1 | person, 
          ~ 1  | group, 
          emergence = ~ -1 + time, # 1 inkluderar "indivudal baseline variance",
                                   # hur tolkar vi den?
                                  # Kontrollerar för individuella grundskillnader? 
                                  # ex mer positiv till att börja med
          time = "time",
          method = "CEI", 
          method.team = "OU",
          data = sherifdat)

summary.ce(CEI.bridge)
## Adjusted CEM
CEI <- ce(y ~ 1+time, 
          ~ 1 | person, 
          ~ 1 + time| group, 
          emergence = ~ -1 + time, # 1 inkluderar "indivudal baseline variance",
          # hur tolkar vi den?
          # Kontrollerar för individuella grundskillnader? 
          # ex mer positiv till att börja med
          method = "CEI", 
          data = sherifdat)

summary.ce(CEI)




## GP

GP <- ce(y ~ 1+time, 
         ~ 1 | person, 
         ~ 1 + time | group, 
         emergence = ~ 1, 
         method = "GP",
         time = "time",
         data = sherifdat)

summary.ce(GP)

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



# r plot
r.plot(CEM, CEI, GP, sherifdat$y,sherifdat$group, sherifdat$time)
r.plot(CEM.bridge, CEI.bridge, GP.bridge, sherifdat$y,sherifdat$group, sherifdat$time)

smooth.plot(CEM, CEI, GP, sherifdat$y, sherifdat$time, groups.to.plot = c(1,3))
smooth.plot(CEM.bridge, CEI.bridge, GP.bridge, sherifdat$y, sherifdat$time, groups.to.plot = c(1,3))
smooth.plot(CEM, CEI.bridge, GP, sherifdat$y, sherifdat$time, groups.to.plot = c(1,8))

smooth.plot(CEM, CEI, GP, sherifdat$y, sherifdat$time, groups.to.plot = c(6,8))
smooth.plot(CEM.bridge, CEI.bridge, GP.bridge, sherifdat$y, sherifdat$time, groups.to.plot = c(6,8))
