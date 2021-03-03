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
  
  # Add extra space to right of plot area
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  
  plot(r$time, r$r, ylim=c(0,1.2), cex=1, pch=19) # possibly change the range of y
  
  legend("right", 
         inset=c(-0.2,0), 
         legend = c("CEM", "CEI", "GP"), 
         lty=1, 
         col = c("red","blue", "green"), 
         title = "Model",
         cex=0.5)
  #
  time <- seq(min(r$time), max(r$time), length.out = 100)
  
  # CEM
  
  r_1 <- exp(paramList1$error[[1]][1] + paramList1$error[[1]][2]*time) + exp(2* paramList1$indv[[1]]) #indv term, same as pdlogchol
  r_1 <- sqrt(r_1/r_1[1])
  lines(time, r_1, col='red')

  
  # CEI

  r_2 <- exp(2*paramList2$indv[[1]][2] + paramList2$indv[[1]][1]*time)  + exp(paramList2$error[[1]])
  r_2 <- sqrt(r_2/r_2[1])
  lines(time, r_2, col='blue')
  
  
  # GP

  Cov <- rep(0, length(time))
  for(i in 1:length(time)){
    indv <- list(D=as.matrix(1),time=time[i])
    Cov[i] <- GP$object$indvCovs[[1]]$get_Cov(paramList3$indv[[1]], indv) + GP$object$indvCovs[[2]]$get_Cov(paramList3$indv[[2]], indv)  
  }
  r_3 <- sqrt(Cov/Cov[1])
  lines(time, r_3, col='green',cex=2)
}

### Smoothing distribution where new observations could land


smooth.plot <- function(CEM, CEI, GP, y, time, type = "total") {
  
  paramList1 <- CEM$covariances # CEM
  paramList2 <- CEI$covariances # CEI
  paramList3 <- GP$covariances # GP
  
  y <- y
  time <- time
  
  range_y <- range(y)
  
  type <- type # total or level2
  
 # betas <-getbeta(CEM$covariances, CEM$object) # not needed?
  Smooth1 <- smoothIndivual(CEM$unlisted_covariances, CEM$object)
  Smooth2 <- smoothIndivual(CEI$unlisted_covariances, CEI$object)
  Smooth3 <- smoothIndivual(GP$unlisted_covariances, GP$object)


  par(mfrow=c(2,3), mai = c(0.5, 0.2, 0.1, 0.2))
  time <- unique(time)
  for(j in 1:2){
    for(i in 1:3){
    
      plot(time, as.vector(Matrix::t(CEM$object$teams[[j]]$indv[[i]]$A)%*%CEM$object$teams[[j]]$data$Y), 
           ylab='obs',ylim=c(range_y[1]-1,range_y[2]+1),cex=1, pch=19)
    
      # CEM
      m <- Smooth1$teams[[j]]$indv[[i]]$mean
      s <-  sqrt(Smooth1$teams[[j]]$indv[[i]]$var  +  exp(paramList1$error[[1]][1] + paramList1$error[[1]][2]*time))
      lines(time,m + 2*s, col='red' )
      lines(time,m - 2*s, col='red' )
    
      # CEI
      m <- Smooth2$teams[[j]]$indv[[i]]$mean
      s <-  sqrt(Smooth2$teams[[j]]$indv[[i]]$var  + exp(paramList2$error[[1]]))
      lines(time,m + 2*s, col='blue' )
      lines(time,m - 2*s, col='blue' )
    
      # GP
      m <- Smooth3$teams[[j]]$indv[[i]]$mean
      s <-  sqrt(Smooth3$teams[[j]]$indv[[i]]$var  +  exp(paramList3$error[[1]][1]))
      lines(time,m + 2*s, col='green' )
      lines(time,m - 2*s, col='green' )
      }
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend('bottom',legend = c("CEM", "CEI", "GP"), col = c("red","blue", "green"), lwd = 5, xpd = TRUE, horiz = TRUE, cex = 1, seg.len=1, bty = 'n')
  }



