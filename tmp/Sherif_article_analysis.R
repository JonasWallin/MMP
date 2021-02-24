### Sherif data analysis for article
### 2021-02-22

library(MMP)

## excluding individual measurement occasions
data("sherifdat")
sherifdat <- subset(sherifdat, time <= 2)

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


## Adjusted CEM
CEI <- ce(y ~ 1+time, 
          ~ 1 | person, 
          ~ 1 + time | group, 
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

# r plot
r.plot(CEM, CEI, GP, sherifdat$y,sherifdat$group, sherifdat$time)

smooth.plot(CEM, CEI, GP, sherifdat$y, sherifdat$time)
