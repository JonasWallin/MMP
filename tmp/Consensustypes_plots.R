## illustration of the different types of consensus assumed by the models
# Date: 2021-04-23

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


