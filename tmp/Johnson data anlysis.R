# Johnson data JoM (2014)
# 2023-09-04
# see coding info in:
# /Dropbox/buisness/Projekt 2 Team Emergence YB JW FD/Johnson et al 2014 JOM
# /Johnson et al - 2014 - JoM - special issue Bayes - team performance with ineq constrained hypo
# /1. Data/1. Raw/PCE coding 01Jan2012.xlsx
library(ggplot2)
library(readxl)
library(dplyr)

PCEraw <- readRDS("tmp/PCEraw.Rdata")

# NOTE: group variable needs to be named group, otherwise smooth plot function
# does not work.
PCEraw$group <- PCEraw$grp

# exclude groups with only first two time points
PCEsub <- subset(PCEraw, 
                 grp!=6 & grp!=7 & grp!=11 & grp!=12
                 & grp!=15 & grp!=16 & grp!=24 & grp!=26
                 & grp!=27 & grp!=29 & grp!=31 & grp!=32
                 & grp!=34 & grp!=35 & grp!=37 & grp!=42
                 & grp!=43 & grp!=52 & grp!=53 & grp!=54
                 & grp!=57)

# recode time to start from 0
PCEsub$time <- PCEsub$time -1

# Variables to look into:
# Goal setting: qs1-4
# Interdependence: qs42-46
# Trust: qs73-75
# Cohesion: qs76-78

vars <- variable.names(PCEsub[,c(9:12,50:53,81:86)])

# fit all models for all 14 variables
output <- list()
for (i in 1:length(vars)) {
  data <- PCEsub[,c("person","group","time",vars[i])]
  names(data) <- c("person","group","time","y")
  data <- subset(data, !is.na(y))
  
  # null model
  null <- ce(y ~ 1+time, 
             ~ 1 | person, 
             ~ 1 + time | group, 
             emergence = ~ 1, 
             method = "CEM2", 
             data = data)
  
  # HetCEM 
  hetcem <- ce(y ~ 1+time, 
               ~ 1 | person, 
               ~ 1 + time | group, 
               emergence = ~ 1 + time, 
               method = "CEM2", 
               data = data)
  
  homcem <- ce(y ~ 1+time, 
               ~ 1 | person, 
               ~ 1 + time | group, 
               emergence = ~ -1 + time, 
               method = "CEI2", 
               data = data)
  
  # GP 
  GP <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 + time | group, 
           emergence = ~ 1, 
           method = "GP",
           time = "time",
           data = data)
  
  ## Nonlinear team dynamics
  # HetCEM
  hetcemNL <- ce(y ~ 1+time, 
                 ~ 1 | person, 
                 ~ 1 | group, 
                 emergence = ~ 1 + time, 
                 time = "time",
                 method = "CEM2",
                 method.team = "OU.homeostasis",
                 data = data)
  
  # HomCEM
  homcemNL <- ce(y ~ 1+time, 
                 ~ 1 | person, 
                 ~ 1 | group, 
                 emergence = ~ -1 + time, 
                 time = "time",
                 method = "CEI2", 
                 method.team = "OU.homeostasis",
                 data = data)
  
  # GP
  GPNL <- ce(y ~ 1+time, 
             ~ 1 | person, 
             ~ 1 | group, 
             emergence = ~ 1, 
             method = "GP",
             method.team = "OU.homeostasis",
             time = "time",
             data = data)
  
  # Akaike weights
  AW <- akaike.weight(list(null,hetcem,homcem,GP,hetcemNL,homcemNL,GPNL), 
                           c("null model","HetCEM","HomCEM","GP","HetCEM NL","HomCEM NL","GP NL"))
  
  # r plot 
  rp <- r.plot(list(hetcem,homcem,GP,hetcemNL,homcemNL,GPNL),
               data$y,
               data$group, 
               data$time, 
               names = c("HetCEM","HomCEM","GP","HetCEM NL","HomCEM NL","GP NL"))
  
  rp <- rp + theme(legend.position=c(0.8, 0.8),
             legend.background = element_rect(size=0.5, linetype="solid",colour ="darkgrey"))
  
  smp <- smooth.plot(models=list(hetcem,homcem,GP,hetcemNL,homcemNL,GPNL),
                     group = subset(PCEsub, !is.na(qs64))$group,
                     person = subset(PCEsub, !is.na(qs64))$person, 
                     time = subset(PCEsub, !is.na(qs64))$time, 
                     y = subset(PCEsub, !is.na(qs64))$qs64,
                     groups.to.plot = c(1,2), 
                     names = c("HetCEM", "HomCEM", "GP","HetCEM NL","HomCEM NL","GP NL")) +
    theme(legend.position = "bottom") +
    scale_y_continuous(breaks=c(0,4,8)) 
  
  output[[i]] <- list(hetcem=hetcem,homcem=homcem,GP=GP,
                      hetcemNL=hetcemNL,homcemNL=homcemNL,GPNL=GPNL,
                      AW=AW,Rplot=rp,Smoothplot=smp)
}

names(output) <- vars
##########

#saveRDS(output, file="Fitted_Models_Johnson_data.RData")
#output <- readRDS("tmp/Fitted_Models_Johnson_data.RData")

# Goal setting
ggplot(PCEsub)+
  geom_line(aes(time,qs1, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs2, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs3, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs4, col=factor(person)))+
  facet_wrap(vars(grp))

# seems to be some consensus here, especially qs3 and 4
# homcem behaves a bit strange?
output[[1]]$Rplot
output[[2]]$Rplot
output[[3]]$Rplot
output[[4]]$Rplot

# Akaike weights
output[[1]]$AW
output[[2]]$AW
output[[3]]$AW
output[[4]]$AW

output[[4]]$Smoothplot

summary.ce(output[[3]]$hetcemNL)
summary.ce(output[[3]]$homcemNL)
summary.ce(output[[3]]$GPNL)

# interdependence
# disensus!
output[[5]]$Rplot
output[[6]]$Rplot
output[[7]]$Rplot
output[[8]]$Rplot

# trust
#disensus
output[[9]]$Rplot
output[[10]]$Rplot
output[[11]]$Rplot

# cohesion
# disensus
output[[12]]$Rplot
output[[13]]$Rplot
output[[14]]$Rplot


### Other variables that may show something
# qs7,qs24,qs33,qs40,qs56,qs58,qs61,qs63,qs64
