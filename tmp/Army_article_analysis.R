## Army data
## 2021-03-25
library(ggplot2)
library(tidyverse)
library(MMP)
# the following is code from Lang et al. 2018 article, with my comments on what
# it does

## get and prepare data
library(multilevel)
data(univbct)

## prepare variables
univbct2<-univbct

# Y: create UNIT variable (batalion + company)
univbct2$UNIT <- paste(univbct2$BTN,univbct2$COMPANY,sep="") 

# Y: Calculate the mean of READY1 (readiness at time point 1) for each unit at 
# the first time point (=0)
RMEANS <- aggregate(READY1 ~ UNIT,univbct2[univbct2$TIME==0,],mean)
names(RMEANS) <- c("UNIT","CREAD")
univbct2 <- merge(univbct2,RMEANS,by="UNIT")

# delete persons with single observations
#univbct2 <- univbct2[rowSums(sapply(subset(univbct2,select=c(JOBSAT1,JOBSAT2,JOBSAT3)),is.na))<2,]

# not in original code: delete persons with missing obs on JSAT 
univbct2$nas.per.row <- rowSums(sapply(subset(univbct2,select=c(JOBSAT1,JOBSAT2,JOBSAT3)),is.na))

totnas <- aggregate(nas.per.row ~ SUBNUM,univbct2,sum)
names(totnas) <- c("SUBNUM","TOTNAS")

univbct2 <- merge(univbct2,totnas,by="SUBNUM")

univbct2 <- univbct2[univbct2$TOTNAS ==0,]

# Y: only use units with at least three members 
members <- table(univbct2$UNIT)/3
univbct2 <- univbct2[univbct2$UNIT %in% names(members[members>2]), ]

# end of code from Lang et al

## Centering
# The within-group means for each time point are
# subtracted from the observed y values.

army <- univbct2[,c(1:10,20:24)]

means <- army %>%
  group_by(UNIT, TIME) %>%
  summarise_at(vars(JSAT), list(~mean(., na.rm=TRUE)))

army <- left_join(army, means, by =c("UNIT", "TIME"))

army$JSAT.centered <- army$JSAT.x-army$JSAT.y

army <- ungroup(army)

army <- rename(army, JSAT = JSAT.x, JSAT.mean = JSAT.y)

army2 <- army[,c(1:4,11:17)]


## recode id variable SUBNUM
units <- unique(army2$UNIT)
n <- army2 %>% 
  count(UNIT)
n$m <- n$n/3

# sort obs
army2 <- army2[order(army2$UNIT,army2$SUBNUM),]
id <- NULL
for (i in 1:nrow(n)) {
  # number from 1-n
  tempid <- rep(1:n$m[i], each=3)
  id <- c(id,tempid)
}

army2 <- cbind(army2,id)

# check, seems to have worked
HHC <- subset(army2,army2$UNIT=="4042SVC")
table(HHC$SUBNUM,HHC$id)

## ANALYSIS

table(army2$UNIT)/3 

# plot of some of the units with not so many members
ggplot(subset(army2,army2$UNIT=="1010F" |army2$UNIT=="1022D" |army2$UNIT=="2004D" 
              |army2$UNIT=="2004HHC" |army2$UNIT=="4000REC" |army2$UNIT=="4042C" 
              |army2$UNIT=="404B" |army2$UNIT=="404HHC" |army2$UNIT=="4B" 
              |army2$UNIT=="4C" |army2$UNIT=="4HHC"),
       aes(x=TIME, y=JSAT.centered)) +
  geom_point(aes(shape=as.factor(id)), show.legend = T) +
  facet_wrap(vars(UNIT), ncol =3) +
  # guides(shape=guide_legend(title="Subjects within each group")) +
  geom_line(aes(col=as.factor(id)), show.legend = F) +
  theme_bw()

# plot of some of the units with  many members
ggplot(subset(army2,army2$UNIT=="1000HHC" |army2$UNIT=="1022A" |army2$UNIT=="1022B" 
              |army2$UNIT=="1022HHC" |army2$UNIT=="124A" |army2$UNIT=="144A" 
              |army2$UNIT=="2008D" |army2$BTN=="299" | army2$BTN=="3066" 
              |army2$UNIT=="4042A" |army2$UNIT=="4042B" |army2$UNIT=="4042SVC" 
              |army2$UNIT=="404A" |army2$UNIT=="4A" |army2$UNIT=="4D"),
       aes(x=TIME, y=JSAT.centered)) +
  geom_point(aes(col=as.factor(id),group=SUBNUM), show.legend = F) +
  facet_wrap(vars(UNIT), ncol = 4) +
  #scale_fill_hue(l=20) +
  geom_line(aes(col=as.factor(id)), show.legend = F) +
  theme_bw()

#########

# check of variance at different time points
var(army2$JSAT.centered[army2$TIME==0])
var(army2$JSAT.centered[army2$TIME==1])
var(army2$JSAT.centered[army2$TIME==2])

# CEM from article
CEM.army <- ce(JSAT ~ 1+TIME, 
          ~ 1 | SUBNUM, 
          ~ 1 + TIME | UNIT, 
          emergence = ~ 1 + TIME, 
          method = "CEM", 
          # method.team = "OU",
          data = army2)

summary.ce(CEM.army)



# null model
null.army <- ce(JSAT ~ 1+TIME, 
               ~ 1 | SUBNUM, 
               ~ 1 + TIME | UNIT, 
               emergence = ~ 1, 
               method = "CEM2", 
               # method.team = "OU",
               data = army2)
summary.ce(null.army)

# HeCEM
HeCEM <- ce(JSAT ~ 1+TIME, 
               ~ 1 | SUBNUM, 
               ~ 1 + TIME | UNIT, 
               emergence = ~ 1 + TIME+I(TIME^2), 
               method = "CEM2", 
               # method.team = "OU",
               data = army2)

r <- r.emperical(army2$JSAT, army2$UNIT, army2$TIME)

time <- seq(min(r$time), max(r$time), length.out = 100)


re <- exp(2* HeCEM$covariances$indv[[1]]) + exp(HeCEM$covariances$indv[[2]][1] + HeCEM$covariances$indv[[2]][2]*time+HeCEM$covariances$indv[[2]][3]*time^2) 
re <- sqrt(re/re[1])

plot(r$time, r$r, ylim=c(0,1.2), cex=1, pch=19)
lines(time, re,col='red')
legend("right", 
       legend = c("HeCEM", "HoCEM"), 
       lty=1, 
       col = c("red","blue"), 
       title = "Model",
       cex=0.5)

# CEI2

  ro <- exp(2* HoCEM$covariances$indv[[1]]) + exp(2*HoCEM$covariances$indv[[2]][3] +(2*HoCEM$covariances$indv[[2]][1]*time)+(2*HoCEM$covariances$indv[[2]][2]*time^2))
  ro <- sqrt(ro/ro[1])
  lines(time, ro,col='blue')

# HoCEM
HoCEM <- ce(JSAT ~ 1+TIME, 
               ~ 1 | SUBNUM, 
               ~ 1 + TIME | UNIT, 
               emergence = ~ -1 + TIME+I(TIME^2), 
               method = "CEI2", 
               # method.team = "OU",
               data = army2)
 
summary.ce(CEI.army)

# GP
GP.army <- ce(JSAT ~ 1+TIME, 
              ~ 1 | SUBNUM, 
              ~ 1 + TIME | UNIT,
              emergence = ~ TIME + I(TIME^2), ## NOT POSSIBLE TO DO THIS
              method = "GP",
              time = "TIME",
              data = army2)

GP <- ce(JSAT ~ 1+TIME, 
              ~ 1 | SUBNUM, 
              ~ 1 + TIME | UNIT,
              emergence = ~ 1, 
              method = "GP",
              time = "TIME",
              data = army2)

## with team GP
# HeCEM
HeCEMGP <- ce(JSAT.centered ~ 1+TIME, 
            ~ 1 | SUBNUM, 
            ~ 1 | UNIT, 
            emergence = ~ 1 + TIME, 
            method = "CEM2", 
            method.team = "OU.homeostasis",
            data = army2)
#Error in optim(param, function(x) { : 
#function cannot be evaluated at initial parameters

# HoCEM
HoCEMGP <- ce(JSAT ~ 1+TIME, 
            ~ 1 | SUBNUM, 
            ~ 1 | UNIT, 
            emergence = ~ -1 + TIME, 
            method = "CEI2", 
            method.team = "OU.homeostasis",
            data = army2)

## Error in optim(param, function(x) { : 
# function cannot be evaluated at initial parameters
# GP
GP.armyGP <- ce(JSAT ~ 1+TIME, 
              ~ 1 | SUBNUM, 
              ~ 1 | UNIT,
              emergence = ~ 1, 
              method = "GP",
              method.team = "OU.homeostasis",
              time = "TIME",
              data = army2)



r.emperical(army2$JSAT.centered,army2$UNIT,army2$TIME)

# r plot
r.plot(CEM.army, CEI.army, GP.army, army2$JSAT,army2$UNIT, army2$TIME)
r.plot(list(HoCEM,HeCEM), army2$JSAT,army2$UNIT, army2$TIME)

# smoothing plot
smooth.plot(CEM.army, CEI.army, GP.army, army2$JSAT, army2$TIME, 
            groups.to.plot = c(1,8))

