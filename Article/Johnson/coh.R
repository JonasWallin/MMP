
##
# move to data processing
## 

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
library(ggpubr)

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

###########################
#File starts



#######
# Figure 5
######
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
print(pl1)

####
# Fitting all the models parameters 
# (using known starting point to imporve speed of convergence)
###
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





###
# TABLE 2
##
names <- c("L:null","L:HetCEM","L:HomCEM","L:GP","GP:HetCEM","GP:HomCEM","GP:GP")
models <- list(null.model,CEM.hem,CEM.hom,GP,CEM.hem.NL,CEM.hom.NL,GPNL)
Table <- akaike.weight(models, 
                       names)
Table <- data.frame(Model = Table$names,
                    AIC =  sapply(models, function(i) i$AIC),
                    loglik = sapply(models, function(i) i$loglik),
                    "Akaike weight" =  Table$weight)
Table$Akaike.weight <- round(Table$Akaike.weight,2)
print(xtable::xtable(Table))


###
# TABLE 4
##

Table4 <- c(GPNL$betas,  #betas
            exp(2 * GPNL$covariances$indv[[1]]), #sigma^2_v0
            exp(2 * GPNL$covariances$indv[[2]][2] - GPNL$covariances$indv[[2]][3])/2, #sigma^2_v1
            GPNL$covariances$indv[[2]][1], #delta_v
            exp( -GPNL$covariances$indv[[2]][3]), #kappa_v
            exp(2 * GPNL$covariances$team[[1]]), # sigma^2_tau0
            exp(2 * GPNL$covariances$team[[2]][2]- GPNL$covariances$team[[2]][3])/2, # sigma^2_tau1
            GPNL$covariances$team[[2]][1], #delta_tau
            exp( -GPNL$covariances$team[[2]][3]), #kappa_tau
            exp(2 * GPNL$covariances$error[[1]]) #sigma2_eps
) 
Table4 <- round(Table4,3)
cat('Table4:\n')
print(Table4)

###
# Figure 6
###
t <- 0:3
dat.fig <- data.frame(t = 0:3, 
                  VeY = diag(Sigma_GPNL$SigmaE),
                  VP = diag(Sigma_GPNL$SigmaI),
                  VG = diag(Sigma_GPNL$SigmaT))
dat.fig$tot <- rowSums(dat.fig[,2:4])
dat.fig$PVeY <- dat.fig$VeY/dat.fig$tot
dat.fig$PVP <- dat.fig$VP/dat.fig$tot
dat.fig$PVG <- dat.fig$VG/dat.fig$tot
dat.fig$PVY <- dat.fig$tot/dat.fig$tot
dat.fig$VY <- dat.fig$tot


pvar <- dat.fig %>% 
  pivot_longer(c(VeY,VP,VG,VY),names_to = "var",values_to = "value") %>% 
  ggplot(aes(x=t,y=value,col=var,linetype=var)) +
  geom_line() +
  geom_point() +
  theme_bw()+
  ylab("Variance") +
  xlab("Time") +
  #  scale_color_discrete(name="",labels = c("VeY" = TeX('$\\V[e^Y_{tij}]$'),"VG" = TeX('$\\V[G_{tij}]$'),"VP" = TeX('$\\V[P_{tij}]$'),"VY" = TeX('$\\V[Y_{tij}]$')))+
  scale_color_discrete(name="Level",labels = c("VeY" = "Measurement","VP" = "Individual","VG" = "Group","VY" = "Total"),breaks=c("VeY","VP","VG","VY"))+
  scale_linetype_discrete(name="Level",labels = c("VeY" = "Measurement","VP" = "Individual","VG" = "Group","VY" = "Total"),breaks=c("VeY","VP","VG","VY"))+
  theme(legend.position = "none")
pvar

pperc <- dat.fig %>% 
  pivot_longer(c(PVeY,PVP,PVG,PVY),names_to = "var",values_to = "value") %>% 
  ggplot(aes(x=t,y=value,col=var,linetype=var)) +
  geom_line() +
  geom_point() +
  theme_bw()+
  ylab("Proportion") +
  xlab("Time") +
  # scale_color_discrete(name="",labels = c("PVeY" = TeX('$\\V[e^Y_{tij}]$'),"PVG" = TeX('$\\V[G_{tij}]$'),"PVP" = TeX('$\\V[P_{tij}]$'),"PVY" = TeX('$\\V[Y_{tij}]$')))+
  scale_color_discrete(name="Level",labels = c("PVeY" = "Measurement","PVP" = "Individual","PVG" = "Group","PVY" = "Total"),breaks=c("PVeY","PVP","PVG","PVY"))+
  scale_linetype_discrete(name="Level",labels = c("PVeY" = "Measurement","PVP" = "Individual","PVG" = "Group","PVY" = "Total"),breaks=c("PVeY","PVP","PVG","PVY"))+
  theme(legend.position = "right")
pperc

varplot2 <- ggpubr::ggarrange(pvar, pperc, 
                      #labels = c("A", "B"),
                      ncol = 2, nrow = 1,
                      common.legend = T,
                      legend="bottom")

