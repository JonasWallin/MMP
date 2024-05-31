


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
library(MASS)
library(misty)
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




data.coh <- PCEsub[,vars[12:14]]

#compute the eigenvalues and eigenvectors of the correlation matrix
SigmaE <- eigen(cov2cor(cov(data.coh,use="complete")))
cat('largest eigenvalue explains =',100*round(SigmaE$values[1]/sum(SigmaE$values),2),'% of the variability\n')
cat('the corresponding eigenvector is:', round(SigmaE$vectors[,1],2),'\n')
results_coh <- process_data(data.coh, "tru")

group.id <- unique(results_coh$data$group)
data.vis <- results_coh$data[results_coh$data$group%in%group.id[21:30],]

pl1 <- ggplot(data=data.vis, 
              aes(x=time,y=y, linetype = factor(person))) + 
  geom_line(size=0.5) + xlab("Time") + ylab("Cohesion") +
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

if(save.fig){
  ggsave("Cohesion.traj.pdf",pl1)
}
print(pl1)


par.null  <- c(-4.2122558, -2.4602522, -4.8964796, -2.0860085, -2.9267910, -0.0183323)
null.model <- ce(y ~ 1+time, 
              ~ 1 | person, 
              ~ 1 + time | group, 
              emergence = ~ 1, 
              method = "CEM2", 
              data = results_coh$data,
              param = par.null)
Sigma.null <- get.Cov(null.model$covariances, null.model$object)
par.CEM  <- c(-3.81459831, -2.46355658, -27.31109528,   6.81364972,
              -2.08420523,  -2.92549732,  -0.01843801)
CEM.hem <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 + time | group, 
           emergence = ~ 1+time, 
           method = "CEM2", 
           data = results_coh$data,
           param = par.CEM)
Sigma_CEM.hem <- get.Cov(CEM.hem$covariances, CEM.hem$object)
par_CEM.hom <- c(-3.84398510, -4.67598509 , 0.24852295,
                 -2.85337209, -2.06347094, -2.95799850, -0.01947531)
CEM.hom <- ce(y ~ 1+time, 
            ~ 1  | person, 
            ~ 1  + time| group, 
            emergence = ~ -1+time, 
            method = "CEI2", 
            time = "time",
            data = results_coh$data,
            par= par_CEM.hom)
Sigma_CEM.hom <- get.Cov(CEM.hom$covariances, CEM.hom$object)
par.GP <- c(-3.90021116, -4.34476048,  0.21589538,
            -3.49020635, -2.23371314, -2.06447437,
            -2.95203464, -0.01863654, -7.90253494)
GP <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 +time | group, 
           emergence = ~ 1 , 
           method = "GP",
           time = "time",
           data =results_coh$data,
           param=par.GP)
Sigma_GP <- get.Cov(GP$covariances, GP$object)
par.GP <- c(-4.2035533, -5.2085661 , 0.1726854 ,-2.9251160,
            -1.5413470 ,-4.3169783 , 0.0823050 ,-2.2834835 ,-1.1748932)
GPNL <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 | group, 
           emergence = ~ 1 , 
           method = "GP",
           method.team = "OU.homeostasis",
           time = "time",
           data =results_coh$data,
           param=par.GP)
Sigma_GPNL <- get.Cov(GPNL$covariances, GPNL$object)
par.CEM.NL <- c(-4.44916866, -2.42968033, -5.27093651,  0.22947725,
                -6.17190719,  0.08602995, -2.30145879, -1.20964578)
CEM.hem.NL <- ce(y ~ 1+time, 
               ~ 1  | person, 
               ~ 1 | group, 
               emergence = ~ 1+time, 
               method.team = "OU.homeostasis",
               method = "CEM2", 
               time = "time",
               data = results_coh$data,
               param = par.CEM.NL)
Sigma_CEM_NL <- get.Cov(CEM.hem.NL$covariances, CEM.hem.NL$object)

par.Random.NL <- c(-4.280932121, -2.582142403, -3.232578980 , 0.002251664 ,-5.924047557,
                   -2.147554920, -0.110134175 ,-0.477223419,2.843943268)
Randomslope.NL <- ce(y ~ 1+time, 
                 ~ 1 + time  | person, 
                 ~ 1 | group, 
                 emergence = ~ 1, 
                 method.team = "OU.homeostasis",
                 method = "CEM2", 
                 time = "time",
                 data = results_coh$data,
                 param = par.Random.NL)
Sigma_RandomSlope_NL <- get.Cov(Randomslope.NL$covariances, Randomslope.NL$object)


par_CEM.hom.NL <- c(-4.01322217, -3.88106536, 0.25777991, -2.83494409,
                    -4.41923045 , 0.07912546, -2.29050619, -1.19022169)
CEM.hom.NL <- ce(y ~ 1+time, 
                 ~ 1  | person, 
                 ~ 1 | group, 
                 emergence = ~ -1+time, 
                 method.team = "OU.homeostasis",
                 method = "CEI2", 
                 time = "time",
                 data = results_coh$data,
                 param = par_CEM.hom.NL)
Sigma_hom_NL <- get.Cov(CEM.hom.NL$covariances, CEM.hom.NL$object)

par.nul.NL <-c(-4.06906401, -2.40471177, -6.28582685,
               -6.69748389,  0.09287034, -2.31271413, -1.21208195)
nul.NL <- ce(y ~ 1+time, 
                 ~ 1  | person, 
                 ~ 1 | group, 
                 emergence = ~ 1, 
                 method.team = "OU.homeostasis",
                 method = "CEM2", 
                 time = "time",
                 data = results_coh$data,
             param = par.nul.NL)
Sigma_NUL_NL <- get.Cov(nul.NL$covariances, nul.NL$object)


#table Loglik AIC 


# bootstrap
#' simulate data using get.Cov covariance
simulate.data <- function(data, Sigmas){
  
  groups <- unique(data$group)
  data.sim <- data
  for(group in groups){
    
    group.data <- data[data$group==group,]
    y.group <-  mvrnorm(n = 1, mu = rep(0,dim(Sigmas$SigmaT)[1]), Sigma = Sigmas$SigmaT)
    indvs <- unique(group.data$person)
    for(indv in indvs){
      y.indv <-  mvrnorm(n = 1, mu = rep(0,dim(Sigmas$SigmaT)[1]), Sigma = Sigmas$SigmaI+Sigmas$SigmaE)
      Y <-  y.indv + y.group
      group.data$y[group.data$person==indv] = Y[match(group.data$time[group.data$person==indv],c(0,1,2,3))] 
      
    }
    data.sim[data$group==group,] <- group.data
  }
  
  
  return(data.sim)
  
}
model.select.lin <- data.frame(name=c("NULL","HEM","HOM","GP"),
                               loglik = c(null.model$loglik,CEM.hem$loglik, CEM.hom$loglik,GP$loglik),
                               AIC    = c(null.model$AIC,CEM.hem$AIC, CEM.hom$AIC,GP$AIC),
                               group=c("multi","multi","multi","multi"))
model.select.non <- data.frame(name=c("HEM","HOM","slope","GP"),
                           loglik = c(nul.NL$loglik, CEM.hom.NL$loglik,Randomslope.NL$loglik, GPNL$loglik),
                           AIC    = c(nul.NL$AIC, CEM.hom.NL$AIC,Randomslope.NL$AIC, GPNL$AIC),
                           group  = c("GP non stationary","GP non stationary", "GP non stationary","GP non stationary") )
model.select <- rbind(model.select.lin,model.select.non)
model.select$wLik= round(exp(model.select$loglik)/sum(exp(model.select$loglik)),2)
print(model.select)
var.Indv.GP <- round(diag(Sigma_GPNL$SigmaI)/Sigma_GPNL$SigmaI[1,1],2)
var.Indv.GP.rel <- round(diag(Sigma_GPNL$SigmaI)/diag(Sigma_GPNL$Sigma),2)
var.Team.GP.rel <- round(diag(Sigma_GPNL$SigmaT)/diag(Sigma_GPNL$Sigma),2)
cat('GP:  V[I_t]/V[I_0] = ',var.Indv.GP,'\n')
cat('GP:  V[I_t]/V[Y_t] = ',var.Indv.GP.rel,'\n')
cat('GP:  V[G_t]/V[Y_t] = ',var.Team.GP.rel,'\n')

var.Indv.HEM <- round(diag(Sigma_hom_NL$SigmaI)/Sigma_hom_NL$SigmaI[1,1],2)
var.Indv.HEM.rel <- round(diag(Sigma_hom_NL$SigmaI)/diag(Sigma_hom_NL$Sigma),2)
var.Team.HEM.rel <- round(diag(Sigma_hom_NL$SigmaT)/diag(Sigma_hom_NL$Sigma),2)
var.Meas.HEM.rel <- round(diag(Sigma_hom_NL$SigmaE)/diag(Sigma_hom_NL$Sigma),2)
cat('HOM: V[I_t]/V[I_0] = ',var.Indv.HEM,'\n')
cat('HOM: V[I_t]/V[Y_t] = ',var.Indv.HEM.rel,'\n')
cat('HOM: V[E_t]/V[Y_t] = ',var.Meas.HEM.rel,'\n')
cat('HOM:  V[G_t]/V[Y_t] = ',var.Team.HEM.rel,'\n')

