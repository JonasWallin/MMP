## Generate data for illustration of delta_Y and delta_G
library(ggplot2)
library(gridExtra)
library(multilevel)
library(psych)
library(dplyr)




gendat <- function(timepoints, persons, groups, g0, sg, sy, dg, dy) {
  dat=expand.grid(time=0:(timepoints-1),
                  person=1:persons,
                  group=1:groups)
  eg <- rep(rnorm(groups, 0, sg), each = timepoints*persons) 
  ey <- rnorm(timepoints*persons*groups, 0, sd=sy*exp(dy*dat$time))
  dat$y <- g0+exp(dg*dat$time)*eg+ey 
  dat$icc1t <- sg^2 / (sg^2 + (sy^2*exp(2*(dy-dg)*dat$time)))
  return(dat)
}



set.seed(123456)
set.seed(121314)
## case 1
data1 <- gendat(timepoints = 4, persons = 4, groups = 4, g0=0, sg=1, sy=1, dg=0, dy=-.6)



##############

data1$person2 <- "0"
data1$person2 <- ifelse(data1$person == 1, "1", data1$person2)
data1$person2 <- ifelse(data1$person == 2, "2", data1$person2)
data1$person2 <- ifelse(data1$person == 3, "3", data1$person2)



pl1 <- ggplot(data=data1, aes(x=time,y=y)) + 
  geom_point(size=1.5, aes(x=time,y=y, shape = person2)) +
  # geom_line(aes(col=person2))+
  xlab("time") + 
  ylab("y") + 
  ggtitle("Case 1: "~delta[P]==-.6~", "~delta[G]==0) + 
  labs(shape = "Person within group") #+
#scale_shape_manual(values=c(21:24))
pl1 <- pl1 +  
  facet_grid(cols=vars(group))  + 
  theme_bw() + 
  guides(shape="none" ) + 
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8))

print(pl1)


## case 2
data2 <- gendat(timepoints = 4, persons = 4, groups = 4, g0=0, sg=1, sy=1, dg=.6, dy=0)

data2$person2 <- "0"
data2$person2 <- ifelse(data2$person == 1, "1", data2$person2)
data2$person2 <- ifelse(data2$person == 2, "2", data2$person2)
data2$person2 <- ifelse(data2$person == 3, "3", data2$person2)

pl2 <- ggplot(data=data2, aes(x=time,y=y)) + 
  geom_point(size=1.5, aes(x=time,y=y, shape = person2)) +
  xlab("time") + 
  ylab("y") + 
  ggtitle("Case 2: "~delta[P]==0~", "~delta[G]==.6) + 
  labs(shape = "Person within group")
pl2 <- pl2 +  
  facet_wrap(~group)  + 
  theme_bw() + guides(shape="none" ) + 
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8))

print(pl2)

grid.arrange(pl1, pl2, ncol=2)
