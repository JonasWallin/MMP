## Generate data for illustration of delta_Y and delta_G

library(ggplot2)
library(gridExtra)
library(ggpubr)

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



set.seed(3)

## case 1
data1 <- gendat(timepoints = 4, persons = 4, groups = 4, g0=0, sg=1, sy=1, dg=0, dy=-.6)

round(unique(data1$icc1t),2)

##############

data1$person2 <- "0"
data1$person2 <- ifelse(data1$person == 1, "1", data1$person2)
data1$person2 <- ifelse(data1$person == 2, "2", data1$person2)
data1$person2 <- ifelse(data1$person == 3, "3", data1$person2)
# data1$person2 <- as.character(data1$person)


pl1 <- ggplot(data=data1, aes(x=time,y=y)) + 
  geom_point(size=1.5, aes(x=time,y=y, shape = person2)) +
  #geom_line(size=0.5, aes(linetype=person2))+
  xlab("Time") + 
  ylab("y") + 
 ggtitle(expression("Case 1: "~delta[P]==-.6~", "~delta[G]==0), 
       subtitle = expression("ICC: "~t[0]==.50~", "~t[1]==.77~", "~t[2]==.92~", "~t[3]==.97)) +
  #        subtitle = "ICC: "~t[0]==.50~", "~t[1]==.77~", "~t[2]==.92~", "~t[3]==.97~) + 
  labs(shape = "Person within group") #+
#scale_shape_manual(values=c(21:24))
pl1 <- pl1 +  
  facet_wrap(vars(group))  + 
  theme_bw() + 
  guides(shape="none" ) + 
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8)) +
  theme(legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        #axis.text.x = element_text(size = 11),
       # axis.text.y = element_text(size = 11),  
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size=12),
        plot.subtitle = element_text(size=10))
        

print(pl1)
pl1 <- ggplot(data=data1, aes(x=time,y=y)) + 
  geom_line(aes(linetype=person2))+
  geom_point(size=1.5, aes(x=time,y=y)) +
  xlab("time") + 
  ylab("y")  
#  ggtitle("Case 2: "~delta[P]==0~", "~delta[G]==.6, 
#          subtitle = expression("ICC: "~t[0]==.50~", "~t[1]==.77~", "~t[2]==.92~", "~t[3]==.97)) + 
#  labs(shape = "Person within group")
pl1 <- pl1 +  
  facet_wrap(~group)  + 
  theme_classic(base_size = 12) + #guides(shape="none" ) + 
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8))+
  scale_color_discrete(guide="none")+
  scale_linetype_discrete(guide="none")
#theme(legend.text=element_text(size=12),
#      legend.title=element_text(size=12),
#axis.text.x = element_text(size = 11),
#axis.text.y = element_text(size = 11),  
#     axis.title.x = element_text(size = 12),
#      axis.title.y = element_text(size = 12),
#      plot.title = element_text(size=12),
#      plot.subtitle = element_text(size=10))
pl1

## case 2
data2 <- gendat(timepoints = 4, persons = 4, groups = 4, g0=0, sg=1, sy=1, dg=.6, dy=0)

round(unique(data2$icc1t),2)

data2$person2 <- "0"
data2$person2 <- ifelse(data2$person == 1, "1", data2$person2)
data2$person2 <- ifelse(data2$person == 2, "2", data2$person2)
data2$person2 <- ifelse(data2$person == 3, "3", data2$person2)

pl2 <- ggplot(data=data2, aes(x=time,y=y)) + 
  geom_point(size=1.5, aes(x=time,y=y, shape = person2)) +
  xlab("Time") + 
  ylab("y") + 
  ggtitle("Case 2: "~delta[P]==0~", "~delta[G]==.6, 
          subtitle = expression("ICC: "~t[0]==.50~", "~t[1]==.77~", "~t[2]==.92~", "~t[3]==.97)) + 
  labs(shape = "Person within group")
pl2 <- pl2 +  
  facet_wrap(~group)  + 
  theme_bw() + guides(shape="none" ) + 
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8))+
  theme(legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        #axis.text.x = element_text(size = 11),
        #axis.text.y = element_text(size = 11),  
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size=12),
        plot.subtitle = element_text(size=10))

print(pl2)

grid.arrange(pl1, pl2, ncol=2)

pl2 <- ggplot(data=data2, aes(x=time,y=y)) + 
  #geom_line(aes(linetype=person2,col=person2))+
  geom_line(aes(linetype=person2))+
#  geom_point(size=1.5, aes(x=time,y=y, col=person2)) +
  geom_point(size=1.5, aes(x=time,y=y)) +
  xlab("time") + 
  ylab("y")  
#  ggtitle("Case 2: "~delta[P]==0~", "~delta[G]==.6, 
#          subtitle = expression("ICC: "~t[0]==.50~", "~t[1]==.77~", "~t[2]==.92~", "~t[3]==.97)) + 
#  labs(shape = "Person within group")
pl2 <- pl2 +  
  facet_wrap(~group)  + 
  theme_classic(base_size = 12) + #guides(shape="none" ) + 
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8))+
  scale_color_discrete(guide="none")+
  scale_linetype_discrete(guide="none")
  #theme(legend.text=element_text(size=12),
  #      legend.title=element_text(size=12),
        #axis.text.x = element_text(size = 11),
        #axis.text.y = element_text(size = 11),  
  #     axis.title.x = element_text(size = 12),
  #      axis.title.y = element_text(size = 12),
  #      plot.title = element_text(size=12),
  #      plot.subtitle = element_text(size=10))
pl2


fig <- ggarrange(pl1, pl2, labels = c("A","B"),
          common.legend = TRUE, legend = "bottom")
ggsave("ICC_conflation_v2.pdf",fig)
