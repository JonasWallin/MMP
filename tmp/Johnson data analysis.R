# Johnson data JoM (2014)
# 2023-09-01
# see coding info in:
# /Dropbox/buisness/Projekt 2 Team Emergence YB JW FD/Johnson et al 2014 JOM
# /Johnson et al - 2014 - JoM - special issue Bayes - team performance with ineq constrained hypo
# /1. Data/1. Raw/PCE coding 01Jan2012.xlsx
library(ggplot2)
library(readxl)
library(dplyr)
PCEraw <- read_excel("~/Library/CloudStorage/Dropbox/buisness/Projekt 2 Team Emergence YB JW FD/Johnson et al 2014 JOM/Johnson et al - 2014 - JoM - special issue Bayes - team performance with ineq constrained hypo/1. Data/1. Raw/PCE dataset 08Jan2012.xlsx", 
                     sheet = "student", na = ".")

PCEraw$person <- PCEraw$ids-6*(PCEraw$grp-1)

# exclude groups with only first two time points
PCEsub <- subset(PCEraw, 
                 grp!=6 & grp!=7 & grp!=11 & grp!=12
                 & grp!=15 & grp!=16 & grp!=24 & grp!=26
                 & grp!=27 & grp!=29 & grp!=31 & grp!=32
                 & grp!=34 & grp!=35 & grp!=37 & grp!=42
                 & grp!=43 & grp!=52 & grp!=53 & grp!=54
                 & grp!=57)

##########
# recode time to start from 0
PCEsub$time <- PCEsub$time -1


null <- ce(qs64 ~ 1+time, 
           ~ 1 | person, 
           ~ 1 + time | grp, 
           emergence = ~ 1, 
           method = "CEM2", 
           data = PCEsub)

# HetCEM 
CEM2 <- ce(qs64 ~ 1+time, 
           ~ 1 | person, 
           ~ 1 + time | grp, 
           emergence = ~ 1 + time, 
           method = "CEM2", 
           data = PCEsub)

# HomCEM 
homcem <- ce(qs64 ~ 1+time, 
             ~ 1 | person, 
             ~ 1 + time | grp, 
             emergence = ~ -1 + time, 
             method = "CEI2", 
             data = PCEsub,
             REML = F)

# GP 
GP <- ce(qs64 ~ 1+time, 
         ~ 1 | person, 
         ~ 1 + time | grp, 
         emergence = ~ 1, 
         method = "GP",
         time = "time",
         data = PCEsub)
########################