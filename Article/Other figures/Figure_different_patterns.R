##
# Reproducing Figure 1
##
graphics.off()
library(ggplot2)
library(MMP)
library(ggpubr)
library(tidyr)
set.seed(22)
save.fig = F
T_     <- 4
delta = -0.4
n.grid <- 500
n.obs  <- 5
Maternparam <- log(c(1,3,2.4)) #sigma, kappa, nu

# creating grid for continous line, and observations location
grid       <- seq(0,T_, length.out = n.grid)
grid.obs   <- c(seq(1,n.grid, by=ceiling(n.grid/n.obs)), n.grid) #location of the grid points!
Sigma.grid <- Materncov(d = as.matrix(dist(grid)), Maternparam)
L.grid     <- t(chol(Sigma.grid))
# create three different
Matern_cov     <- MaternConsenus$new(1)
Matern_input   <-  list(D =    grid,
                        time = grid) 
Sigma_con  <- Matern_cov$get_Cov(c(delta,Maternparam), Matern_input)

Z1 <- rnorm(n.grid)
Z2 <- rnorm(n.grid)


# Plot GP
L.con      <- t(chol(Sigma_con))
x1 <-  L.con%*%Z1
x2 <- L.con%*%Z2
fig.data = data.frame(time = grid, lower.CI = -1.96*exp(delta*grid), upper.CI =1.96*exp(delta*grid),
                      x1 = x1, x2 = x2)
fig.obs <- data.frame(time = grid[grid.obs], x1 = x1[grid.obs], x2 = x2[grid.obs])
fig.obs <- gather(fig.obs, "individual", "obs", -time)

pl1 <- ggplot(data=fig.data) + 
  geom_line(aes(x=time,y=x1 ),size=1,alpha=0.3,col='red')+ geom_line(aes(x=time,y=x2 ),size=1,alpha=0.3,col='red')  +
#  geom_point(data=fig.obs,aes(x=time,y=x1 ),size=2.5,col='red',shape=15) +geom_point(data=fig.obs,aes(x=time,y=x2 ),size=2.5,col='red')+
  geom_point(data=fig.obs,aes(x=time,y=obs, shape=individual ),size=2.5,col='red') +
  geom_line(aes(x = time,  y = lower.CI),color='blue',size=1.5,linetype = "dashed")+
  geom_line(aes(x = time, y = upper.CI),color='blue',size=1.5,linetype = "dashed")+
  xlab("Time") + 
  theme_bw()+
  theme(legend.text=element_text(size=12))+
  scale_shape_discrete(name="",labels = c("Individual 1", "Individual 2"))+
  ylab("")
print(pl1)

# plot Het
x1.het <-  exp(delta*grid)*Z1
x2.het <- exp(delta*grid)*Z2
fig.data = data.frame(time = grid, lower.CI = -1.96*exp(delta*grid), upper.CI =1.96*exp(delta*grid),
                      x1 = x1.het, x2 = x2.het)
fig.obs <- data.frame(time = grid[grid.obs], x1 = x1.het[grid.obs], x2 = x2.het[grid.obs])
fig.obs <- gather(fig.obs, "individual", "obs", -time)

pl2 <- ggplot(data=fig.data) + 
  geom_point(data=fig.obs,aes(x=time,y=obs, shape=individual ),size=2.5,col='red') +
#  geom_point(data=fig.obs,aes(x=time,y=x1 ),size=2.5,col='red',shape=15) +geom_point(data=fig.obs,aes(x=time,y=x2 ),size=2.5,col='red')+
  geom_line(aes(x = time,  y = lower.CI),color='blue',size=1.5,linetype = "dashed")+
  geom_line(aes(x = time, y = upper.CI),color='blue',size=1.5,linetype = "dashed")+
  xlab("Time") + 
  theme_bw()+
  theme(legend.text=element_text(size=12))+
  scale_shape_discrete(name="",labels = c("Individual 1", "Individual 2"))+
  ylab("")
print(pl2)

# plot Hom
x1.hom <-  exp(delta*grid)*Z1[1]
x2.hom <- exp(delta*grid)*Z2[1]
fig.data = data.frame(time = grid, lower.CI = -1.96*exp(delta*grid), upper.CI =1.96*exp(delta*grid),
                      x1 = x1.hom, x2 = x2.hom)
fig.obs <- data.frame(time = grid[grid.obs], x1 = x1.hom[grid.obs], x2 = x2.hom[grid.obs])

fig.obs <- gather(fig.obs, "individual", "obs", -time)
pl3 <- ggplot(data=fig.data) + 
  geom_line(aes(x=time,y=x1 ),size=1,alpha=0.3,col='red')+ 
  geom_line(aes(x=time,y=x2 ),size=1,alpha=0.3,col='red')  +
  geom_point(data=fig.obs,aes(x=time,y=obs, shape=individual ),size=2.5,col='red') +
#  geom_point(data=fig.obs,aes(x=time,y=x1 ),size=2.5,col='red',shape=15) +
#  geom_point(data=fig.obs,aes(x=time,y=x2 ),size=2.5,col='red')+
  geom_line(aes(x = time,  y = lower.CI),color='blue',size=1.5,linetype = "dashed")+
  geom_line(aes(x = time, y = upper.CI),color='blue',size=1.5,linetype = "dashed")+
  xlab("Time") + 
  theme_bw()+
  theme(legend.text=element_text(size=12))+
  scale_shape_discrete(name="",labels = c("Individual 1", "Individual 2"))+
  ylab("")


plot_all <- ggarrange(pl3,pl1,pl2, 
                      labels = c("Homogeneous", "Gaussian Process", "Heterogeneous"),
                      ncol = 3, nrow = 1,
                      common.legend = T,
                      vjust = 2,
                      hjust = c(-0.7,-0.6,-0.7),
                      font.label = list(size = 14, color = "black", face = "bold", family = NULL),
                      legend="bottom")
print(plot_all)

