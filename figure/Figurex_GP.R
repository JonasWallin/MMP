##
# simulate n indivual GP then fit 
# then same processes with conesus modeling
# 
##
graphics.off()
library(ggplot2)
library(MMP)
set.seed(6)
save.fig = T
n.indv <- 4
T_     <- 4
n.grid <- 500
n.obs  <- 5
Maternparam <- log(c(1.,1.,2.5))
sigma_Y <- 0.1
grid       <- seq(0,T_, length.out = n.grid)
grid.obs   <- seq(1,n.grid, by=ceiling(n.grid/n.obs))
Sigma.grid <- Materncov(d = as.matrix(dist(grid)), Maternparam)
L.grid     <- t(chol(Sigma.grid))
Matern_cov     <- MaternBridge$new(0, T_,1)
Matern_input   <-  list(D =    as.matrix(1),
                    time = grid) 
Sigma_con  <- Matern_cov$get_Cov(c(log(T_+1),Maternparam), Matern_input)
L.con      <- t(chol(Sigma_con))
Indivauls <- list()
Z <- c(0,0.2)
Latent      <- c()
Latent.con  <- c()
Obs         <- c()
Obs.con         <- c()
for(i in 1:n.indv){
  e <- sigma_Y*rnorm(n.obs)
  eps <- rnorm(n= n.grid)
  Indivauls$W[[i]]      <- L.grid%*%eps
  Indivauls$Wcon[[i]]   <- L.con%*%eps
  Indivauls$Z[[i]]   <- Z[1] + grid*Z[2]
  Indivauls$obs[[i]] <- sample(1:n.grid,n.obs)
  Indivauls$Y[[i]]   <- Indivauls$W[[i]][Indivauls$obs[[i]]] + 
                        Indivauls$Z[[i]][Indivauls$obs[[i]]] + e
  Indivauls$Ycon[[i]]   <- Indivauls$Wcon[[i]][Indivauls$obs[[i]]] + 
    Indivauls$Z[[i]][Indivauls$obs[[i]]] + e
                        
  Latent <- rbind(Latent, cbind(grid,
                                Indivauls$W[[i]] + Indivauls$Z[[i]],
                                i,0))
  Obs     <- rbind(Obs, cbind(grid[Indivauls$obs[[i]]], 
                              Indivauls$Y[[i]],
                              i,0))
  Latent.con <- rbind(Latent.con, cbind(grid,
                                Indivauls$Wcon[[i]] + Indivauls$Z[[i]],
                                i,1))
  Obs.con     <- rbind(Obs.con, cbind(grid[Indivauls$obs[[i]]], 
                              Indivauls$Ycon[[i]],
                              i,1))
}
Obs.data   <- data.frame(time  = Obs[,1],
                         y  = Obs[,2],
                         I  = Obs[,3],
                         c  = Obs[,4])
Obs.data.con   <- data.frame(time  = Obs.con[,1],
                         y  = Obs.con[,2],
                         I  = Obs.con[,3],
                         c  = Obs.con[,4])
Latent.data <- data.frame(time  = Latent[,1],
                          U  = Latent[,2],
                          I  = Latent[,3],
                          c  = Latent[,4])
Latent.data.con <- data.frame(time  = Latent.con[,1],
                          U  = Latent.con[,2],
                          I  = Latent.con[,3],
                          c  = Latent.con[,4])



pl1 <- ggplot(data=Latent.data, 
              aes(x=time,y=U )) + 
  geom_line(size=0.5,alpha=0.2) + 
  geom_point(data=Obs.data, aes(x=time,y=y) ,size=3)+ xlab("time") + ylab("y,u") 
pl1 <- pl1 +  facet_wrap(~I, ncol = 2) 
if(save.fig)
  ggsave('Figure_GP_obs.pdf', pl1)
print(pl1)
con.data <- rbind(Latent.data,Latent.data.con)
con.data$c <- as.factor(con.data$c)
pl2 <- ggplot(data=con.data, 
              aes(x=time,y=U, color=c,
              linetype=c)) + 
       geom_line(size=1,alpha=1) +
  scale_color_manual(name=expression(delta),values = c("black", "red"),
                     labels = c(expression(infinity), T_+1))+
       scale_linetype_manual(name=expression(delta),
                             values=c("dashed", "solid"), 
                             labels = c(expression(infinity), T_+1)) +
       xlab("time") + ylab("y,u") 
pl2 <- pl2 +  facet_wrap(~I, ncol = 2) 
print(pl2)
if(save.fig)
  ggsave('Figure_GP_cons.pdf', pl2)
