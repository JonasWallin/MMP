## Plot of Sherifdata
## 2021-04-06

data("sherifdat")

sherifdat$time <- sherifdat$time + 1


pl1 <- ggplot(data=sherifdat, 
              aes(x=time,y=y, 
                  colour = as.factor(person), 
                  shape = factor(person, labels = c("Subject 1", "Subject 2", "Subject 3")))) + 
  geom_line(size=0.5) +geom_point(size=2)+ xlab("time") + ylab("inches") 


pl1 <- ggplot(data=sherifdat, 
              aes(x=time,y=y, linetype = factor(person, labels = c("Subject 1", "Subject 2", "Subject 3")))) + 
  geom_line(size=0.5) + xlab("time") + ylab("inches") +
  guides(linetype=guide_legend(title="Subjects within each group")) +
  scale_x_continuous(sec.axis = sec_axis(~.*1,labels = c("indv", "group", "group", "group", "indv" )))



pl1 <- pl1 +  facet_wrap(~group, ncol = 2) + 
  labs(title = "Sherif (1935) autokinetic data") + 
  theme_bw() 


print(pl1)
