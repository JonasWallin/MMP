### reproduction of Sherif data analysis 



library(MMP)
library(tidyverse)
library(ggplot2)
library(xtable)
library(latex2exp)
library(ggpubr)

## Data preparation
data("sherifdat")
sherifdat$time <- sherifdat$time + 1

#######
#Figure 3
#######
fig1 <- ggplot(data=sherifdat, 
               aes(x=time,y=y, linetype = factor(person, labels = c("Subject 1", "Subject 2", "Subject 3")))) + 
  geom_line(size=0.5) + xlab("Time") + ylab("Inches") +
  guides(linetype=guide_legend(title="Subjects within each group")) +
  scale_x_continuous(sec.axis = sec_axis(~.*1,labels = c("indv", "group", "group", "group", "indv" )))

fig1 <- fig1 +  facet_wrap(~group, ncol = 4) + 
  # labs(title = "Sherif (1935) autokinetic data") + 
  theme_bw() +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        #axis.text.x = element_text(size = 11),
        #axis.text.y = element_text(size = 11),  
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))


print(fig1)



#### 
# Fitting all the models
####
# excluding last time point (individual measurement)
sherifdat <- subset(sherifdat, time <= 3)

## Linear group dynamics
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

## Gaussian process team dynamics
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

#########
#Table 1
#########
names <- c("L:null","L:HetCEM","L:HomCEM","L:GP","GP:HetCEM","GP:HomCEM","GP:GP")
models <- list(null,CEM2,CEI2,GP,CEM2.h,CEI2.h,GP.h)
Table <- akaike.weight(models, 
                         names)
Table <- data.frame(Model = Table$names,
                    AIC =  sapply(models, function(i) i$AIC),
                    loglik = sapply(models, function(i) i$loglik),
                    "Akaike weight" =  Table$weight)
Table$Akaike.weight <- round(Table$Akaike.weight,2)
print(xtable::xtable(Table))

###########
# Table 3 (appendix)
##########

Table3 <- c(CEI2.h$betas,  #betas
            exp(2 * CEI2.h$covariances$indv[[1]]), #sigma^2_v0
            exp(2 * CEI2.h$covariances$indv[[2]][2]), #sigma^2_v1
            CEI2.h$covariances$indv[[2]][1], #delta_v
            exp(2 * CEI2.h$covariances$team[[1]]), # sigma^2_tau0
            exp(2 * CEI2.h$covariances$team[[2]][2] - CEI2.h$covariances$team[[2]][3])/2, # sigma^2_tau1
            CEI2.h$covariances$team[[2]][1], #delta_tau
            exp(- CEI2.h$covariances$team[[2]][3]), #kappa_tau
            exp(2 * CEI2.h$covariances$error[[1]]) #sigma2_eps
            )  
Table3 <- round(Table3,3)
cat('Table3:\n')
print(Table3)


##########
#Figure 4
##########
#generating the variances of each component
Covs <- get.Cov(CEI2.h$covariances,CEI2.h$object)

dat <- data.frame(t = 0:3, 
                  VeY = diag(Covs$SigmaE),
                  VP = diag(Covs$SigmaI),
                  VG = diag(Covs$SigmaT))
dat$tot <- rowSums(dat[,2:4])
dat$PVeY <- dat$VeY/dat$tot
dat$PVP <- dat$VP/dat$tot
dat$PVG <- dat$VG/dat$tot
dat$PVY <- dat$tot/dat$tot
dat$VY <- dat$tot


pvar <- dat %>% 
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

pperc <- dat %>% 
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

Figure4 <- ggpubr::ggarrange(pvar, pperc, 
                      #labels = c("A", "B"),
                      ncol = 2, nrow = 1,
                      common.legend = T,
                      legend="bottom")
print(Figure4)