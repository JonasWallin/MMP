## illustration of the different types of consensus assumed by the models
# date: 2021-06-16
# last updated: 2022-03-16

library(dplyr)
library(ggplot2)

# cei = Homogeneous Consensus Emergence
# cem = Heterogeneous Consensus Emergence

x4 <- rep(8,4)
x3 <- c(10,9,7,6)
x2 <- c(12,10,6,4)
x1 <- c(14,11,5,2)


cei <- data.frame(cbind(y=c(x1,x2,x3,x4),p=rep(1:4, 4),time=rep(1:4,each=4),type="Homogeneous CE"))

set.seed(123)
y4 <- x4
y3 <- sample(x3,4)
y2 <- sample(x2,4)
y1 <- sample(x1,4)

cem <- data.frame(cbind(y=c(y1,y2,y3,y4),p=rep(1:4, 4),time=rep(1:4,each=4),type="Heterogeneous CE"))

dat <- bind_rows(cei,cem)
dat$y <- as.integer(dat$y)

# Names for facets
time.labs <- c("Time point 1","Time point 2","Time point 3","Time point 4")
names(time.labs) <- c(1,2,3,4)

ggplot(dat) + 
  geom_col(aes(p,y)) + 
  facet_grid(cols=vars(time),rows = vars(type),labeller = labeller(time=time.labs)) + 
  theme_bw() +
  labs(title ="Consensus emergence patterns") +
  xlab("Individual") +
  ylab("Y")


# for presentation
ggplot(dat) + 
  geom_col(aes(p,y)) + 
  facet_grid(cols=vars(time),rows = vars(type),labeller = labeller(time=time.labs)) + 
  theme_bw() +
  #labs(title ="Consensus emergence patterns") +
  xlab("Individual") +
  ylab("Y")

ggplot(dat) + 
  geom_col(aes(p,y)) + 
  geom_abline(aes(p,y), intercept = 8, slope = 0) +
  facet_grid(cols=vars(time),rows = vars(type),labeller = labeller(time=time.labs)) + 
  theme_bw() +
  #labs(title ="Consensus emergence patterns") +
  xlab("Individual") +
  ylab("Y")

