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

grps <- PCEraw %>% 
  group_by(grp) %>% 
  summarise(length(unique(ids)))

idnumbers <- PCEraw %>% 
  group_by(grp) %>% 
  summarise(min(ids),max(ids))

PCEraw$person <- PCEraw$ids-6*(PCEraw$grp-1)

idnumbers <- PCEraw %>% 
  group_by(grp) %>% 
  summarise(min(ids),max(ids),min(person),max(person))

# exclude groups with only first two time points
PCEsub <- subset(PCEraw, 
                 grp!=6 & grp!=7 & grp!=11 & grp!=12
                 & grp!=15 & grp!=16 & grp!=24 & grp!=26
                 & grp!=27 & grp!=29 & grp!=31 & grp!=32
                 & grp!=34 & grp!=35 & grp!=37 & grp!=42
                 & grp!=43 & grp!=52 & grp!=53 & grp!=54
                 & grp!=57)

idn <- PCEsub %>% 
  group_by(grp) %>% 
  summarise(min(ids),max(ids),min(person),max(person))

# goal setting
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

# interdependence
ggplot(PCEsub)+
  geom_line(aes(time,qs42, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs43, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs44, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs45, col=factor(person)))+
  facet_wrap(vars(grp))

# trust
ggplot(PCEsub)+
  geom_line(aes(time,qs73, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs74, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs75, col=factor(person)))+
  facet_wrap(vars(grp))

# cohesion
ggplot(PCEsub)+
  geom_line(aes(time,qs76, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs77, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs78, col=factor(person)))+
  facet_wrap(vars(grp))

# The rest of the variables
ggplot(PCEsub)+
  geom_line(aes(time,qs5, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs6, col=factor(person)))+
  facet_wrap(vars(grp))

# maybe this one
ggplot(PCEsub)+
  geom_line(aes(time,qs7, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs8, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs9, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs10, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs11, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs12, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs13, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs14, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs15, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs16, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs17, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs18, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs19, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs20, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs21, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs22, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs23, col=factor(person)))+
  facet_wrap(vars(grp))

# maybe this one: (probably not)
ggplot(PCEsub)+
  geom_line(aes(time,qs24, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs25, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs26, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs27, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs28, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs29, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs30, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs31, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs32, col=factor(person)))+
  facet_wrap(vars(grp))

# maybe this one:
ggplot(PCEsub)+
  geom_line(aes(time,qs33, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs34, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs35, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs36, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs37, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs38, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs39, col=factor(person)))+
  facet_wrap(vars(grp))

# maybe this one:
ggplot(PCEsub)+
  geom_line(aes(time,qs40, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs41, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs47, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs48, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs49, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs50, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs51, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs52, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs53, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs54, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs55, col=factor(person)))+
  facet_wrap(vars(grp))

# maybe this one:
ggplot(PCEsub)+
  geom_line(aes(time,qs56, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs57, col=factor(person)))+
  facet_wrap(vars(grp))

# maybe this one:
ggplot(PCEsub)+
  geom_line(aes(time,qs58, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs59, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs60, col=factor(person)))+
  facet_wrap(vars(grp))

# maybe this one
ggplot(PCEsub)+
  geom_line(aes(time,qs61, col=factor(person)))+
  facet_wrap(vars(grp))


ggplot(PCEsub)+
  geom_line(aes(time,qs62, col=factor(person)))+
  facet_wrap(vars(grp))

# maybe:
ggplot(PCEsub)+
  geom_line(aes(time,qs63, col=factor(person)))+
  facet_wrap(vars(grp))

# maybe this one!
ggplot(PCEsub)+
  geom_line(aes(time,qs64, col=factor(person)))+
  facet_wrap(vars(grp))


ggplot(PCEsub)+
  geom_line(aes(time,qs65, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs66, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs67, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs68, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs69, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs70, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs71, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs72, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs79, col=factor(person)))+
  facet_wrap(vars(grp))


ggplot(PCEsub)+
  geom_line(aes(time,qs80, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs81, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs82, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs83, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs84, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs85, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs86, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs87, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs88, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs89, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs90, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs91, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs92, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs93, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs94, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs95, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs96, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs97, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs98, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs99, col=factor(person)))+
  facet_wrap(vars(grp))

ggplot(PCEsub)+
  geom_line(aes(time,qs100, col=factor(person)))+
  facet_wrap(vars(grp))
