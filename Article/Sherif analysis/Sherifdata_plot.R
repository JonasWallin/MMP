## Plot of Sherifdata
## 2021-04-06

data("sherifdat")

sherifdat$time <- sherifdat$time + 1

pl1 <- ggplot(data=sherifdat, 
              aes(x=time,y=y, linetype = factor(person, labels = c("Subject 1", "Subject 2", "Subject 3")))) + 
  geom_line(size=0.5) + xlab("Time") + ylab("Inches") +
  guides(linetype=guide_legend(title="Subjects within each group")) +
  scale_x_continuous(sec.axis = sec_axis(~.*1,labels = c("indv", "group", "group", "group", "indv" )))

pl1 <- pl1 +  facet_wrap(~group, ncol = 4) + 
 # labs(title = "Sherif (1935) autokinetic data") + 
  theme_bw() +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        #axis.text.x = element_text(size = 11),
        #axis.text.y = element_text(size = 11),  
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))


print(pl1)

## for presentation
pl2 <- ggplot(data=subset(sherifdat, time <= 3), 
              aes(x=time,y=y, linetype = factor(person, labels = c("Subject 1", "Subject 2", "Subject 3")))) + 
  geom_line(size=0.5) + xlab("Time") + ylab("Inches") +
  guides(linetype=guide_legend(title="Subjects within each group")) 
  #scale_x_continuous(sec.axis = sec_axis(~.*1,labels = c("indv", "group", "group", "group", "indv" )))

pl2 <- pl2 +  facet_wrap(~group, ncol = 4) + 
  theme_bw() +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        #axis.text.x = element_text(size = 11),
        #axis.text.y = element_text(size = 11),  
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))


print(pl2)
