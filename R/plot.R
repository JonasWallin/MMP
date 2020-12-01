#### Plots

r.plot <- function(CEM, CEI, GP, y, group, time, type = "total") {
  
  paramList1 <- CEM$covariances # CEM
  paramList2 <- CEI$covariances # CEI
  paramList3 <- GP$covariances # GP
  
  y <- y
  group <- group
  time <- time
  
  type <- type # total or level2
  
  ###
  #'
  #' computing r for emperical and the three models
  #'
  ###
  
  r <- r.emperical(y, group, time)
  
  plot(r$time, r$r, ylim=c(0,1.2), cex=2, pch=19) # possibly change the range of y
  #
  time <- seq(min(r$time), max(r$time), length.out = 100)
  
  # CEM
  
  r_1 <- exp(paramList1$error[[1]][1] + paramList1$error[[1]][2]*time) + exp(2* paramList1$indv[[1]]) #indv term, same as pdlogchol
  r_1 <- sqrt(r_1/r_1[1])
  lines(time, r_1, col='red')

  
  # CEI

  r_2 <- exp(2*paramList2$indv[[1]][1] + paramList2$indv[[1]][2]*time) + exp(paramList2$error[[1]]) # this seems more reasonable, kolla med JOnas
  r_2 <- sqrt(r_2/r_2[1])
  lines(time, r_2, col='blue')
  
  
  # GP

  Cov <- rep(0, length(time))
  for(i in 1:length(time)){
    indv <- list(D=as.matrix(1),time=time[i])
    Cov[i] <- GP$object$indvCovs[[2]]$get_Cov(paramList3$indv[[2]], indv)  + exp(paramList3$error[[1]])
  }
  r_3 <- sqrt(Cov/Cov[1])
  lines(time, r_3, col='green',cex=2)
}

###### try to find what's wrong
Cov <- rep(0, length(time))
for(i in 1:length(time)){
  indv <- list(D=as.matrix(1),time=time[4])
  Cov[4] <- GP$object$indvCovs[[2]]$get_Cov(paramList3$indv[[2]], indv)  + exp(paramList3$error[[1]])
}



r_3 <- sqrt(Cov/Cov[1])
lines(time, r_3, col='green',cex=2)


#
# r.plot(sherif1, sherif4, sherif8, sherifdat$y, sherifdat$group, sherifdat$time)
##
y <- sherifdat$y
group <- sherifdat$group
time <- sherifdat$time
CEM <- sherif1
CEI <- sherif4
GP <- sherif8
paramList3 <- GP$covariances # GP

#########
## Smooth plot
  betas <-getbeta(res$par, TeamObj)
  Smooth1 <- smoothIndivual(res1$par, object1)
  Smooth2 <- smoothIndivual(res2$par, object2)
  Smooth3 <- smoothIndivual(res3$par, object3)
  
  
  
  par(mfrow=c(2,3))
  time <- c(0,1,2)
  for(j in 1:2){
    for(i in 1:3){
      plot(time, as.vector(Matrix::t(object2$teams[[j]]$indv[[i]]$A)%*%object2$teams[[j]]$data$Y), ylab='obs',ylim=c(-1.5,1.5) ,cex=2, pch=19)
      
      # CEI
      m <- Smooth2$teams[[j]]$indv[[i]]$mean
      s <-  sqrt(Smooth2$teams[[j]]$indv[[i]]$var  + exp(paramList2$error[[1]]))
      lines(time,m + 2*s, col='blue' )
      lines(time,m - 2*s, col='blue' )
      
      # CEM
      m <- Smooth1$teams[[j]]$indv[[i]]$mean
      s <-  sqrt(Smooth1$teams[[j]]$indv[[i]]$var  +  exp(paramList1$error[[1]][1] + paramList1$error[[1]][2]*time))
      lines(time,m + 2*s, col='red' )
      lines(time,m - 2*s, col='red' )
      
      # GP
      m <- Smooth3$teams[[j]]$indv[[i]]$mean
      s <-  sqrt(Smooth3$teams[[j]]$indv[[i]]$var  +  exp(paramList3$error[[1]][1] ))
      lines(time,m + 2*s, col='green' )
      lines(time,m - 2*s, col='green' )
    }
  }
  


