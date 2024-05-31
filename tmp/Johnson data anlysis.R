# Johnson data JoM (2014)
# 2023-09-04
# see coding info in:
# /Dropbox/buisness/Projekt 2 Team Emergence YB JW FD/Johnson et al 2014 JOM
# /Johnson et al - 2014 - JoM - special issue Bayes - team performance with ineq constrained hypo
# /1. Data/1. Raw/PCE coding 01Jan2012.xlsx
library(MMP)
library(ggplot2)
library(readxl)
library(dplyr)
#library(misty)
data.folder <- "../../Dropbox/articles/buisness/Projekt 2 Team Emergence YB JW FD/Johnson et al 2014 JOM/data/Data/Formatted/"
dat <- read.csv(paste(data.folder,"pce_data_wide.dat",sep=""),sep=",",header=F)
varnames <- read.mplus(input =  paste(data.folder,"PCE_data_wide.inp",sep=""), return.var = TRUE)
colnames(dat) <- varnames
data.group <- dat[,c("grp","prf1","prf3","prf6","prf9","prf13","prf16","prf17","prf19","prf21")]
group.low <-unique(data.group[data.group$prf21<=12,"grp"])
process_data <- function(data, var_name) {
  # Combine selected columns from PCEsub with the scaled row means of `data`
  new_data <- cbind(PCEsub[, c("person", "group", "time")], rowMeans(scale(data, center = FALSE)))
  
  # Rename the last column to "y"
  names(new_data)[4] <- "y"
  
  # Remove rows where "y" is NA
  new_data <- subset(new_data, !is.na(y))
  
  # Calculate the empirical correlation
  r_result <- r.emperical(new_data$y, new_data$group, new_data$time)
  
  # Optionally, return the cleaned data and result as a list (if needed for further analysis)
  list(data = new_data, result = r_result)
}


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
#
vars <- variable.names(PCEsub[,c(9:12,50:53,81:86)])

# fit all models for all 14 variables
#i=3,4,5,11,14

data <- PCEsub[,c("person","group","time",vars)]
#goal
i <- 4:7 #goal
data.goal <- PCEsub[,vars[1:4]]
data.ind <- PCEsub[,vars[5:8]]
data.tru <- PCEsub[,vars[9:11]]
data.coh <- PCEsub[,vars[12:14]]

results_goal <- process_data(data.goal, "goal")
results_ind <- process_data(data.ind, "ind")
results_tru <- process_data(data.tru, "tru")
results_coh <- process_data(data.coh, "coh")


group.id <- unique(results_coh$data$group)
data.vis <- results_goal$data[results_goal$data$group%in%group.id[1:10],]
data.vis <- results_ind$data[results_ind$data$group%in%group.id[11:20],]
data.vis <- results_coh$data[results_coh$data$group%in%group.id[21:30],]
data.vis <- results_coh$data[results_coh$data$group%in%group.id[31:40],]
pl1 <- ggplot(data=data.vis, 
              aes(x=time,y=y, linetype = factor(person))) + 
  geom_line(size=0.5) + xlab("Time") + ylab("Inches") +
  guides(linetype=guide_legend(title="Subjects within each group")) 

pl1 <- pl1 +  facet_wrap(~group, ncol = 5) + 
  theme_bw() +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        #axis.text.x = element_text(size = 11),
        #axis.text.y = element_text(size = 11),  
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
print(pl1)

data <- results_tru$data

results_tru$data[,"AboveMedian"]=1.
results_tru$data[results_tru$data$group%in%group.low,"AboveMedian"]=0

#check for time series behaviour
for (i in 1:length(vars)) {

  # null model
  # fit the model with there version
  null <- ce(y ~ 1+time, 
             ~ 1 | person, 
             ~ 1 + time | group, 
             emergence = ~ 1+time, 
             method = "CEM2", 
             data = results_coh$data)
  null2 <- ce(y ~ 1+time, 
             ~ 1 + time | person, 
             ~ 1 + time | group, 
             emergence = ~ 1+time, 
             method = "CEM2", 
             data = results_coh$data)
  
  null3 <- ce(y ~ 1+time, 
              ~ 1  | person, 
              ~ 1 | group, 
              emergence = ~ 1+time, 
              method.team = "OU.homeostasis",
              method = "CEM2", 
              time = "time",
              data = results_coh$data)
  
  par0 <- null3$res$par
  par0[3:4] <- c(-6,2*0.19)
  par0 <- c(-4.25994754,-2.42425435,-5.88825259,0.34774913,-2.16421238,0.04909911,4.29736046,12.80444426)
  null3.v2 <- ce(y ~ 1+time, 
              ~ 1  | person, 
              ~ 1 | group, 
              emergence = ~ 1+time, 
              method.team = "OU.homeostasis",
              method = "CEM2", 
              time = "time",
              data = results_coh$data,
              param = par0)
  Sigma_null3 <- get.Cov(null3.v2$covariances,null3.v2$object)
  ## Nonlinear team dynamics
  # HetCEM
  hetcemNL <- ce(y ~ 1+time, 
                 ~ 1 | person, 
                 ~ 1 | group, 
                 emergence = ~ 1 + time:AboveMedian, 
                 time = "time",
                 method = "CEM2",
                 method.team = "OU.homeostasis",
                 data = results_tru$data)
  
  # HomCEM
  homcemNL <- ce(y ~ 1+time, 
                 ~ 1 | person, 
                 ~ 1 | group, 
                 emergence = ~ -1 + time, 
                 time = "time",
                 method = "CEI2", 
                 method.team = "OU.homeostasis",
                 data = data.hig)
  
  # GP
  GPNL <- ce(y ~ 1+time, 
                 ~ 1 | person, 
                 ~ 1 | group, 
                 emergence = ~ 1 + time, 
                 method = "GP",
                 method.team = "OU.homeostasis",
                 time = "time",
                 data = results_coh$data)
  data.low <- results_tru$data[results_tru$data$group%in%group.low,]
  data.hig <- results_tru$data[results_tru$data$group%in%group.low ==F,]
  GPNL.low <- ce(y ~ 1+time, 
             ~ 1 | person, 
             ~ 1 | group, 
             emergence = ~ 1, 
             method = "GP",
             method.team = "OU.homeostasis",
             time = "time",
             data = data.low)
  paramList <- GPNL.high$covariances
  Meas_error.h        <- GPNL.high$object$errorCovs[[1]]$get_Cov(paramList$error[[1]],list(E=as.matrix(rep(1,4))))
  Sigma_indv_1.h       <- GPNL.high$object$indvCovs[[1]]$get_Cov(paramList$indv[[1]], list(D=as.matrix(0:3),time=0:3))
  Sigma_indv_2_delta.h <- GPNL.high$object$indvCovs[[2]]$get_Cov(paramList$indv[[2]], list(D=as.matrix(0:3),time=0:3))
  
  Sigma_team_1.h       <- GPNL.high$object$teamCovs[[1]]$get_Cov(paramList$team[[1]], list(D=as.matrix(0:3),time=0:3))
  Sigma_team_2_delta.h <- GPNL.high$object$teamCovs[[2]]$get_Cov(paramList$team[[2]], list(D=as.matrix(0:3),time=0:3))
  
  high_tot  = Meas_error.h + c(Sigma_indv_1.h+Sigma_team_1.h)*outer(rep(1,4),rep(1,4))+Sigma_indv_2_delta.h + Sigma_team_2_delta.h
  paramList <- GPNL.low$covariances
  Meas_error.l        <- GPNL.low$object$errorCovs[[1]]$get_Cov(paramList$error[[1]],list(E=as.matrix(rep(1,4))))
  Sigma_indv_1.l       <- GPNL.low$object$indvCovs[[1]]$get_Cov(paramList$indv[[1]], list(D=as.matrix(0:3),time=0:3))
  Sigma_indv_2_delta.l <- GPNL.low$object$indvCovs[[2]]$get_Cov(paramList$indv[[2]], list(D=as.matrix(0:3),time=0:3))
  
  Sigma_team_1.l       <- GPNL.low$object$teamCovs[[1]]$get_Cov(paramList$team[[1]], list(D=as.matrix(0:3),time=0:3))
  Sigma_team_2_delta.l <- GPNL.low$object$teamCovs[[2]]$get_Cov(paramList$team[[2]], list(D=as.matrix(0:3),time=0:3))
  low_tot  = Meas_error.l + c(Sigma_indv_1.l+Sigma_team_1.l)*outer(rep(1,4),rep(1,4)) +Sigma_indv_2_delta.l+ Sigma_team_2_delta.l
  
  
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
