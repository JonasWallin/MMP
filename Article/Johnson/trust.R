


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


data.tru <- PCEsub[,vars[9:11]]


results_tru <- process_data(data.tru, "tru")

data.vis <- results_tru$data[results_tru$data$group%in%group.id[1:10],]

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

par.CEM2  <- c(-3.84852279,-2.404715839,-30.693516114,6.77164903,-2.462483834,-3.188759865,0.004174169)
CEM2 <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 + time | group, 
           emergence = ~ 1+time, 
           method = "CEM2", 
           data = results_tru$data,
           param = par.CEM2)
Sigma_CEM2 <- get.Cov(CEM2$covariances, CEM2$object)

par.CEM2.v2 <- c(-3.848522793,-2.404715839,-30.693516114,
                 6.771649032,-2.462483834 ,-3.188759865, 0.004174169)
CEM2.v2 <- ce(y ~ 1+time, 
            ~ 1 + time | person, 
            ~ 1  | group, 
            emergence = ~ 1+time, 
            method = "CEM2", 
            data = results_tru$data,
            param = par.CEM2.v2)
Sigma_CEM2.v2 <- get.Cov(CEM2.v2$covariances, CEM2.v2$object)
par.GP <- c(-4.71486561,-2.40836018,0.09448087,0.58832019,5.31555305,-2.29690975,0.08903990,-0.72426726 ,3.47921197)

GPNL <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 | group, 
           emergence = ~ 1 , 
           method = "GP",
           method.team = "OU.homeostasis",
           time = "time",
           param = par.GP,
           data =results_tru$data)