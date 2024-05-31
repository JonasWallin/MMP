

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



data.goal <- PCEsub[,vars[1:4]]


#compute the eigenvalues and eigenvectors of the correlation matrix
SigmaE <- eigen(cov2cor(cov(data.goal,use="complete")))
cat('largest eigenvalue explains =',100*round(SigmaE$values[1]/sum(SigmaE$values),2),'% of the variability\n')
cat('the corresponding eigenvector is:', round(SigmaE$vectors[,1],2),'\n')
results_goal <- process_data(data.goal, "goal")

data.vis <- results_goal$data[results_goal$data$group%in%group.id[1:10],]

pl1 <- ggplot(data=data.vis, 
              aes(x=time,y=y, linetype = factor(person))) + 
  geom_line(size=0.5) + xlab("Time") + ylab("Goal") +
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

par.null  <- c(-4.2122558, -2.4602522, -4.8964796, -2.0860085, -2.9267910, -0.0183323)
null.model <- ce(y ~ 1+time, 
                 ~ 1 | person, 
                 ~ 1 + time | group, 
                 emergence = ~ 1, 
                 method = "CEM2", 
                 data = results_goal$data,
                 param = par.null)


GPNL <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 | group, 
           emergence = ~ 1 , 
           method = "GP",
           method.team = "OU.homeostasis",
           time = "time",
           data =results_goal$data)
Sigma_GPNL <- get.Cov(GPNL$covariances, GPNL$object)