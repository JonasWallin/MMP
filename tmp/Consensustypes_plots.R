## illustration of the different types of consensus assumed by the models
# date: 2021-06-16

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



# Date: 2021-04-23 (old versions)

y <- c(8,7,3,1,
       7,5,4,2,
       6,4.75,4,3,
       4.5,4.5,4.5,4.5)
i <- rep(c("1","2", "3","4"),4)
t <- rep(1:4, each = 4)

cei <- as.data.frame(cbind(i,y,t))


ggplot(cei,aes(t,y, group=i)) + 
  geom_point(shape = as.factor(i),show.legend = T) + 
  geom_line(col = t) + 
  theme_bw() +
  labs(title ="CEI type of consensus") +
  xlab("Time") +
  ylab("Y")

ggplot(cei) + 
  geom_col(aes(i,y)) + 
  facet_grid(cols=vars(t)) + 
  theme_bw() +
  labs(title ="CEI type of consensus") +
  xlab("Individual") +
  ylab("Y")


z <- c(8,7,3,1,
       3,2,5,6,
       5,6,4,3,
       4.5,4.5,4.5,4.5)

cem <- as.data.frame(cbind(i,y,t))

ggplot(cem,aes(t,z, group=i)) + 
  geom_point(shape = as.factor(i))+ 
  geom_line(col = t) + 
  theme_bw() +
  labs(title ="CEM type of consensus") +
  xlab("Time") +
  ylab("Y")

ggplot(cem) + 
  geom_col(aes(i,z)) + 
  facet_grid(cols=vars(t)) + 
  theme_bw() +
  labs(title ="CEM type of consensus") +
  xlab("Individual") +
  ylab("Y")


