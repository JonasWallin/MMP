## Army data
## 2021-03-25
library(ggplot2)

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
univbct2 <- univbct2[rowSums(sapply(subset(univbct2,select=c(JOBSAT1,JOBSAT2,JOBSAT3)),is.na))<2,]

# Y: only use units with at least three members 
members <- table(univbct2$UNIT)/3
univbct2 <- univbct2[univbct2$UNIT %in% names(members[members>2]), ]

# end of code from Lang et al

## Centering
# The within-group means for each time point are
# subtracted from the observed y values.

army <- univbct2[,c(1:8,19:24)]

means <- army %>%
  group_by(UNIT, TIME) %>%
  summarise_at(vars(JSAT), list(~mean(., na.rm=TRUE)))

army <- left_join(army, means, by =c("UNIT", "TIME"))

army$JSAT.centered <- army$JSAT.x-army$JSAT.y

army <- ungroup(army)

army <- rename(army, JSAT = JSAT.x, JSAT.mean = JSAT.y)


army2 <- army[,c(1:3,9:16)]
army2 <- na.omit(army2)
army <- na.omit(army) # obs maybe remove some variables before doing this



## ANALYSIS

table(army2$UNIT)/3 # something might be off here, should give integers?

# plot of some of the units with not so many memebers
ggplot(subset(army2,army2$BTN=="299"),aes(x=TIME, y=JSAT.centered)) +
  geom_point(aes(col=as.factor(SUBNUM),group=SUBNUM), show.legend = F) +
  facet_wrap(vars(UNIT)) +
  #scale_fill_hue(l=20) +
  #geom_smooth(method = lm, se=F,aes(col=as.factor(SUBNUM)), show.legend = F) +
  theme_bw()

# CEM from article
CEM.army <- ce(JSAT.centered ~ 1+TIME, 
          ~ 1 | SUBNUM, 
          ~ 1 + TIME | UNIT, 
          emergence = ~ 1 + TIME, 
          method = "CEM", 
          # method.team = "OU",
          data = army2)

summary.ce(CEM.army)

CEM.army.null <- ce(JSAT ~ 1+TIME, 
               ~ 1 | SUBNUM, 
               ~ 1 + TIME | UNIT, 
               emergence = ~ 1, 
               method = "CEM", 
               data = army2)

summary.ce(CEM.army.null)


# CEI ## something is wrong here
CEI.bridge.army <- ce(JSAT.centered ~ 1+TIME, 
                 ~ 1 | SUBNUM, 
                 ~ 1 | UNIT, 
                 emergence = ~ -1 + TIME, 
                 time = "TIME",
                 method = "CEI", 
                 method.team = "OU",
                 data = army2)

 summary.ce(CEI.bridge)
 
 # GP
 GP.army <- ce(JSAT.centered ~ 1+TIME, 
          ~ 1 | SUBNUM, 
          ~ 1 | UNIT, 
          emergence = ~ 1, 
          method = "GP",
          method.team = "OU",
          time = "TIME",
          data = army2)
 
 summary.ce(GP.army)

 r.emperical(army2$JSAT,army2$UNIT,army2$TIME)


