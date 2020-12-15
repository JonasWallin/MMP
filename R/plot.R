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



