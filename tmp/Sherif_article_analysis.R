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



null <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 + time | group, 
           emergence = ~ 1, 
           method = "CEM2", 
           data = sherifdat)

CEM2 <- ce(y ~ 1+time, 
               ~ 1 | person, 
               ~ 1 + time | group, 
               emergence = ~ 1 + time, 
               method = "CEM2", 
               data = sherifdat)
CEI2 <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 + time | group, 
           emergence = ~ -1 + time, # to only extend with P_ij^0, still keep -1
           method = "CEI2", 
           data = sherifdat,
           REML = F)

GP <- ce(y ~ 1+time, 
         ~ 1 | person, 
         ~ 1 + time | group, 
         emergence = ~ 1, 
         method = "GP",
         time = "time",
         data = sherifdat)

CEM2.h <- ce(y ~ 1+time, 
             ~ 1 | person, 
             ~ 1 | group, 
             emergence = ~ 1 + time, 
             time = "time",
             method = "CEM2",
             method.team = "OU.homeostasis",
             data = sherifdat)


CEI2.h <- ce(y ~ 1+time, 
             ~ 1 | person, 
             ~ 1 | group, 
             emergence = ~ -1 + time, # 1 inkluderar "indivudal baseline variance"
             time = "time",
             method = "CEI2", 
             method.team = "OU.homeostasis",
             data = sherifdat)

GP.h <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 | group, 
           emergence = ~ 1, 
           method = "GP",
           method.team = "OU.homeostasis",
           time = "time",
           data = sherifdat)

weigths <- akaike.weight(list(CEM2,CEI2,GP,CEM2.h,CEI2.h,GP.h), c("CEM2","CEI2","GP","CEM2.h","CEI2.h","GP.h"))
round(weigths[,5],2)

# r plot 
r.plot(list(GP,CEI2, CEM2,GP.h,CEI2.h, CEM2.h),y, group, time, sherifdat, names = c("GP","HoCEM","HeCEM","GP+GP","HoCEM+GP","HeCEM+GP"))
r.plot(list(GP,CEI2, CEM2),y, group, time, sherifdat, names = c("GP","HoCEM","HeCEM"))

r.plot_old2(CEM, CEI, GP, sherifdat$y,sherifdat$group, sherifdat$time)
r.plot_old2(CEM.homeostasis, CEI.h, GP.h, sherifdat$y,sherifdat$group, sherifdat$time)

# smoothing plots TODO: ADJUST TO NEW METHODS
smooth.plot(CEM2, CEI2, GP, sherifdat$y, sherifdat$time, groups.to.plot = c(1,8))
smooth.plot(CEM.homeostasis, CEI.h, GP.h, sherifdat$y, sherifdat$time, groups.to.plot = c(1,8))



### same analysis on centered data


null <- ce(y.centered ~ 1+time, 
           ~ 1 | person, 
           ~ 1 + time | group, 
           emergence = ~ 1, 
           method = "CEM2", 
           data = sherifdat)

CEM2 <- ce(y.centered ~ 1+time, 
           ~ 1 | person, 
           ~ 1 + time | group, 
           emergence = ~ 1 + time, 
           method = "CEM2", 
           data = sherifdat)
CEI2 <- ce(y.centered ~ 1+time, 
           ~ 1 | person, 
           ~ 1 + time | group, 
           emergence = ~ -1 + time, # to only extend with P_ij^0, still keep -1
           method = "CEI2", 
           data = sherifdat)

GP <- ce(y.centered ~ 1+time, 
         ~ 1 | person, 
         ~ 1 + time | group, 
         emergence = ~ 1, 
         method = "GP",
         time = "time",
         data = sherifdat)

CEM2.h <- ce(y.centered ~ 1+time, 
             ~ 1 | person, 
             ~ 1 | group, 
             emergence = ~ 1 + time, 
             time = "time",
             method = "CEM2",
             method.team = "OU.homeostasis",
             data = sherifdat)


CEI2.h <- ce(y.centered ~ 1+time, 
             ~ 1 | person, 
             ~ 1 | group, 
             emergence = ~ -1 + time, # 1 inkluderar "indivudal baseline variance"
             time = "time",
             method = "CEI2", 
             method.team = "OU.homeostasis",
             data = sherifdat)

GP.h <- ce(y.centered ~ 1+time, 
           ~ 1 | person, 
           ~ 1 | group, 
           emergence = ~ 1, 
           method = "GP",
           method.team = "OU.homeostasis",
           time = "time",
           data = sherifdat)

weigths <- akaike.weight(list(CEM2,CEI2,GP,CEM2.h,CEI2.h,GP.h), c("CEM2","CEI2","GP","CEM2.h","CEI2.h","GP.h"))
round(weigths[,5],2)

# r plot 
r.plot(list(GP,CEI2, CEM2,GP.h,CEI2.h, CEM2.h),y, group, time, sherifdat, names = c("GP","HoCEM","HeCEM","GP+GP","HoCEM+GP","HeCEM+GP"))

r.plot_old2(CEM, CEI, GP, sherifdat$y,sherifdat$group, sherifdat$time)
r.plot_old2(CEM.homeostasis, CEI.h, GP.h, sherifdat$y,sherifdat$group, sherifdat$time)
