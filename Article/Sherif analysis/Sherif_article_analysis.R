### Sherif data analysis for article
### date: 2021-02-22

## Things to rename:
# CEM2 = HetCEM
# CEI2 = HomCEM
# CEM2.h = HetCEM with nonlinear team dynamics
# CEI2.h = HomCEM with nonlinear team dynamics
# GP.h = GP model with nonlinear team dynamics
# method.team name "OU.homeostasis" = "nonlinear" (or something similar)

# Note: use summary.ce(model_name) to obtain a model summary,
# e.g. summary.ce(CEM2)

library(MMP)
library(tidyverse)

## Data preparation
data("sherifdat")
# excluding last time point (individual measurement)
sherifdat <- subset(sherifdat, time <= 2)
# recode time to start from 0
sherifdat$time <- sherifdat$time + 1

## Linear team dynamics
# null model 
null <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 + time | group, 
           emergence = ~ 1, 
           method = "CEM2", 
           data = sherifdat)

# HetCEM 
CEM2 <- ce(y ~ 1+time, 
               ~ 1 | person, 
               ~ 1 + time | group, 
               emergence = ~ 1 + time, 
               method = "CEM2", 
               data = sherifdat)

# HomCEM 
CEI2 <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 + time | group, 
           emergence = ~ -1 + time, 
           method = "CEI2", 
           data = sherifdat,
           REML = F)

# GP 
GP <- ce(y ~ 1+time, 
         ~ 1 | person, 
         ~ 1 + time | group, 
         emergence = ~ 1, 
         method = "GP",
         time = "time",
         data = sherifdat)

## Nonlinear team dynamics
# HetCEM
CEM2.h <- ce(y ~ 1+time, 
             ~ 1 | person, 
             ~ 1 | group, 
             emergence = ~ 1 + time, 
             time = "time",
             method = "CEM2",
             method.team = "OU.homeostasis",
             data = sherifdat)

# HomCEM
CEI2.h <- ce(y ~ 1+time, 
             ~ 1 | person, 
             ~ 1 | group, 
             emergence = ~ -1 + time, 
             time = "time",
             method = "CEI2", 
             method.team = "OU.homeostasis",
             data = sherifdat)

# GP
GP.h <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 | group, 
           emergence = ~ 1, 
           method = "GP",
           method.team = "OU.homeostasis",
           time = "time",
           data = sherifdat)

# Akaike weights
weigths <- akaike.weight(list(CEM2,CEI2,GP,CEM2.h,CEI2.h,GP.h), 
                         c("CEM2","CEI2","GP","CEM2.h","CEI2.h","GP.h"))
round(weigths[,5],2)

# r plot 
r.plot(list(GP,CEI2, CEM2,GP.h,CEI2.h, CEM2.h),
       sherifdat$y,sherifdat$group, sherifdat$time, 
       names = c("GP","HomCEM","HetCEM","GP GP","HomCEM GP","HetCEM GP"))

## smoothing plots 
# Linear team dynamics
smooth.plot(models=list(CEM2,CEI2,GP),
             group = sherifdat$group,
             person = sherifdat$person, 
             time = sherifdat$time, 
             y = sherifdat$y,
             groups.to.plot = c(1,8), 
             names = c("HetCEM", "HomCEM", "GP"))

# Nonlinear team dynamics
smooth.plot(models=list(CEM2.h,CEI2.h,GP.h),
            group = sherifdat$group,
            person = sherifdat$person, 
            time = sherifdat$time, 
            y = sherifdat$y,
            groups.to.plot = c(1,8), 
            names = c("HetCEM GP", "HomCEM GP", "GP GP"))

# separate plots for each model with vs without nonlinear team dynamics
# (not included in article)
smooth.plot(models=list(CEM2,CEM2.h),
            group = sherifdat$group,
            person = sherifdat$person, 
            time = sherifdat$time, 
            y = sherifdat$y,
            groups.to.plot = c(1,8), 
            names = c("HetCEM","HetCEM GP"))

smooth.plot(models=list(CEI2,CEI2.h),
            group = sherifdat$group,
            person = sherifdat$person, 
            time = sherifdat$time, 
            y = sherifdat$y,
            groups.to.plot = c(1,8), 
            names = c("HomCEM","HomCEM GP"))

smooth.plot(models=list(GP,GP.h),
            group = sherifdat$group,
            person = sherifdat$person, 
            time = sherifdat$time, 
            y = sherifdat$y,
            groups.to.plot = c(1,8), 
            names = c("GP","GP GP"))


