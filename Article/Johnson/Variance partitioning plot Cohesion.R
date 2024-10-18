## Cohesion variance partitioning plot


library(ggplot2)
library(tidyverse)
library(latex2exp)
library(ggpubr)

t <- 0:3
VeY <- rep(0.015,4)
VP <- c(0.007,0.01,0.014,0.02)
VG <- c(0.017,0.02,0.023,0.028)
dat <- data.frame(t, VeY,VP,VG)
dat$tot <- rowSums(dat[,2:4])
dat$PVeY <- dat$VeY/dat$tot
dat$PVP <- dat$VP/dat$tot
dat$PVG <- dat$VG/dat$tot
dat$PVY <- dat$tot/dat$tot
dat$VY <- dat$tot

pvar <- dat %>% 
  pivot_longer(c(VeY,VP,VG,VY),names_to = "var",values_to = "value") %>% 
  ggplot(aes(x=t,y=value,col=var,linetype=var)) +
  geom_line() +
  geom_point() +
  theme_bw()+
  ylab("Variance") +
  xlab("Time") +
  #  scale_color_discrete(name="",labels = c("VeY" = TeX('$\\V[e^Y_{tij}]$'),"VG" = TeX('$\\V[G_{tij}]$'),"VP" = TeX('$\\V[P_{tij}]$'),"VY" = TeX('$\\V[Y_{tij}]$')))+
  scale_color_discrete(name="Level",labels = c("VeY" = "Measurement","VP" = "Individual","VG" = "Group","VY" = "Total"),breaks=c("VeY","VP","VG","VY"))+
  scale_linetype_discrete(name="Level",labels = c("VeY" = "Measurement","VP" = "Individual","VG" = "Group","VY" = "Total"),breaks=c("VeY","VP","VG","VY"))+
  theme(legend.position = "none")
pvar

pperc <- dat %>% 
  pivot_longer(c(PVeY,PVP,PVG,PVY),names_to = "var",values_to = "value") %>% 
  ggplot(aes(x=t,y=value,col=var,linetype=var)) +
  geom_line() +
  geom_point() +
  theme_bw()+
  ylab("Proportion") +
  xlab("Time") +
  # scale_color_discrete(name="",labels = c("PVeY" = TeX('$\\V[e^Y_{tij}]$'),"PVG" = TeX('$\\V[G_{tij}]$'),"PVP" = TeX('$\\V[P_{tij}]$'),"PVY" = TeX('$\\V[Y_{tij}]$')))+
  scale_color_discrete(name="Level",labels = c("PVeY" = "Measurement","PVP" = "Individual","PVG" = "Group","PVY" = "Total"),breaks=c("PVeY","PVP","PVG","PVY"))+
  scale_linetype_discrete(name="Level",labels = c("PVeY" = "Measurement","PVP" = "Individual","PVG" = "Group","PVY" = "Total"),breaks=c("PVeY","PVP","PVG","PVY"))+
  theme(legend.position = "right")
pperc

varplot2 <- ggarrange(pvar, pperc, 
                      #labels = c("A", "B"),
                      ncol = 2, nrow = 1,
                      common.legend = T,
                      legend="bottom")

#ggsave("cohesion_varianceplot2.pdf",varplot2)  