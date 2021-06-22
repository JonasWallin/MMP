#### Plots


library(dplyr)
library(ggplot2)

rt <- function(x, time) {
  model <- x$model$Method
  # CEM2
  if (model == "CEM2") {
    paramList <- x$covariances
    rx <- exp(2* paramList$indv[[1]]) + exp(paramList$indv[[2]][1] + paramList$indv[[2]][2]*time) 
    rx <- sqrt(rx/rx[1])
  }
  # CEI2
  if (model == "CEI2") {
    paramList <- x$covariances
    rx <- exp(2* paramList$indv[[1]]) + exp(2*paramList$indv[[2]][2] +(2*paramList$indv[[2]][1]*time))
    rx <- sqrt(rx/rx[1])
  }
  # GP
  if (model == "GP") {
    paramList <- x$covariances
    Cov <- rep(0, length(time))
    for(i in 1:length(time)){
      indv <- list(D=time[i],time=time[i])
      Cov[i] <- x$object$indvCovs[[1]]$get_Cov(paramList$indv[[1]], indv) + x$object$indvCovs[[2]]$get_Cov(paramList$indv[[2]], indv)  
    }
    rx <- sqrt(Cov/Cov[1])
  }
  r_temp <- as.data.frame(cbind(time, r = rx))
  return(r_temp)
}


r.plot <- function(models, y, group, time, names = NULL) {
  
  
  r <- r.emperical(y, group, time)

  time <- seq(min(r$time), max(r$time), length.out = 100)
  
  dat <- sapply(models, simplify= FALSE, function(i,...) {
    
    temp <- rt(i,time)
    
  })
  

  dat <- bind_rows(dat, .id = "column_label")

  if (!is.null(names)) {
    dat$Model <- rep(names,each=100)
  } else {
    dat$Model <- paste("Model",dat$column_label)
  }
  
  r <- rename(r, r_empirical = r)
  dat <- bind_rows(dat, r)

  # plot
  ggplot(dat) +
    geom_line(mapping = aes(time, r, linetype = Model, color = Model), na.rm = T) +
    geom_point(mapping = aes(time, r_empirical),na.rm = T) +
    scale_linetype_discrete(na.translate=FALSE) +
    scale_color_discrete(na.translate=FALSE) +
    theme_bw() +
    labs(title = "r(t) plot") +
    xlab("Time") 
}





### Smoothing distribution where new observations could land
smooth.plot <- function(CEM, CEI, GP, y, time, type = "total", groups.to.plot = c(1,2)) {
  
  paramList1 <- CEM$covariances # CEM
  paramList2 <- CEI$covariances # CEI
  paramList3 <- GP$covariances # GP
  
  y <- y
  time <- time
  
  Smooth1 <- smoothIndivual(CEM$unlisted_covariances, CEM$object)
  Smooth2 <- smoothIndivual(CEI$unlisted_covariances, CEI$object)
  Smooth3 <- smoothIndivual(GP$unlisted_covariances, GP$object)
  
  time <- unique(time)
  data_ci <- data.frame()
  data_obs <- data.frame()
  for(j in groups.to.plot){
    for(i in 1:3){
      if(length(as.vector(Matrix::t(CEM$object$teams[[j]]$indv[[i]]$A)%*%CEM$object$teams[[j]]$data$Y)) != length(time)){
        if (as.vector(Matrix::t(CEM$object$teams[[j]]$indv[[i]]$A)%*%CEM$object$teams[[j]]$data$Y)[1]==time[1]) {
          time_i = time[1:(length(time)-1)]
        } else {
          time_i = time[2:length(time)]
        }
      }else{
        time_i = time
      }
      
      df <- as.data.frame(cbind(time_i,as.vector(Matrix::t(CEM$object$teams[[j]]$indv[[i]]$A)%*%CEM$object$teams[[j]]$data$Y)))
      df <- rename(df, y_obs = V2)
      df$person <- i
      df$group <- j
      
      # CEM
      m <- Smooth1$teams[[j]]$indv[[i]]$mean
      s <-  sqrt(Smooth1$teams[[j]]$indv[[i]]$var  +  exp(paramList1$error[[1]][1] + paramList1$error[[1]][2]*time_i))
      df1 <- as.data.frame(cbind(m + 2*s,m - 2*s,time_i))
      df1 <- rename(df1, ul = V1, ll = V2)
      
      
      # CEI
      m <- Smooth2$teams[[j]]$indv[[i]]$mean
      s <-  sqrt(Smooth2$teams[[j]]$indv[[i]]$var  + exp(paramList2$error[[1]]))
      df2 <- as.data.frame(cbind(m + 2*s,m - 2*s,time_i))
      df2 <- rename(df2, ul = V1, ll = V2)
      
      # GP
      m <- Smooth3$teams[[j]]$indv[[i]]$mean
      s <-  sqrt(Smooth3$teams[[j]]$indv[[i]]$var  +  exp(paramList3$error[[1]][1]))
      df3 <- as.data.frame(cbind(m + 2*s,m - 2*s,time_i))
      df3 <- rename(df3, ul = V1, ll = V2)
      
      temp_obs <- df
      data_obs <- rbind(data_obs,temp_obs)
      temp_data <-  bind_rows("CEM"=df1, "CEI"=df2, "GP"=df3, .id = "Model")
      temp_data$person <- i
      temp_data$group <- j
      data_ci <- rbind(data_ci,temp_data)
      
    }
  }
  
  # with different colors and linetypes
  ggplot(data_obs) +
    geom_point(aes(time_i,y_obs)) +
    facet_grid(rows = vars(group), cols = vars(person),labeller = label_both) +
    geom_line(data =data_ci, 
              mapping = aes(time_i, ul, color = Model,linetype = Model)) +
    geom_line(data =data_ci, 
              mapping = aes(time_i, ll, color = Model, linetype = Model)) +
    theme_bw() +
    labs(title = "Title") +
    xlab("Time") +
    ylab("Y") +
    scale_color_manual(values=c("#6BA2FF", "#F57970","#27C152"))
  
}

r.plot_old2 <- function(CEM, CEI, GP, y, group, time, type = "level2") {
  
  paramList1 <- CEM$covariances # CEM
  paramList2 <- CEI$covariances # CEI
  paramList3 <- GP$covariances # GP
  
  y <- y
  group <- group
  time <- time
  
  #type <- type # total or level2 #TODO
  ## should be model dependent, not choose for all?
  
  ###
  #'
  #' computing r for emperical and the three models
  #'
  ###
  
  r <- r.emperical(y, group, time)
  
  #
  time <- seq(min(r$time), max(r$time), length.out = 100)
  
  # CEM
  
  r_1 <- exp(paramList1$error[[1]][1] + paramList1$error[[1]][2]*time) + exp(2* paramList1$indv[[1]]) #indv term, same as pdlogchol
  r_1 <- sqrt(r_1/r_1[1])
  
  # CEI
  
  r_2 <- exp(2*paramList2$indv[[1]][2] + 2*paramList2$indv[[1]][1]*time)  + exp(paramList2$error[[1]])
  r_2 <- sqrt(r_2/r_2[1])
  
  # GP
  
  Cov <- rep(0, length(time))
  for(i in 1:length(time)){
    indv <- list(D=as.matrix(1),time=time[i])
    Cov[i] <- GP$object$indvCovs[[1]]$get_Cov(paramList3$indv[[1]], indv) + GP$object$indvCovs[[2]]$get_Cov(paramList3$indv[[2]], indv)  
  }
  r_3 <- sqrt(Cov/Cov[1])
  
  # put together in one data frame
  r_1 <- as.data.frame(cbind(time, r_1))
  r_1 <- rename(r_1, r = r_1)
  r_2 <- as.data.frame(cbind(time, r_2))
  r_2 <- rename(r_2, r = r_2)
  r_3 <- as.data.frame(cbind(time, r_3))
  r_3 <- rename(r_3, r = r_3)
  r <- rename(r, r_empirical = r)
  
  data <- bind_rows("CEM"=r_1, "CEI"=r_2, "GP"=r_3, .id = "Model")
  data <- bind_rows(data, r)
  
  # plot
  ggplot(data) +
    geom_line(mapping = aes(time, r, linetype = Model, color = Model), na.rm = T) +
    geom_point(mapping = aes(time, r_empirical),na.rm = T) +
    scale_linetype_discrete(na.translate=FALSE) +
    scale_color_manual(values=c("#6BA2FF", "#F57970","#27C152"), na.translate = F) +
    theme_bw() +
    labs(title = "r(t) plot") +
    xlab("Time") 
  
}

r.plot_old <- function(CEM, CEI, GP, y, group, time, type = "total") {
  
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


smooth.plot_old <- function(CEM, CEI, GP, y, time, type = "total", groups.to.plot = c(6,8)) {
  
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
  for(j in groups.to.plot){
    for(i in 1:3){
      if(length(as.vector(Matrix::t(CEM$object$teams[[j]]$indv[[i]]$A)%*%CEM$object$teams[[j]]$data$Y)) ==3){
        time_i = time[2:4]
      }else{
        time_i = time[1:4]
      }
      plot(time_i, as.vector(Matrix::t(CEM$object$teams[[j]]$indv[[i]]$A)%*%CEM$object$teams[[j]]$data$Y), 
           ylab='obs',ylim=c(range_y[1]-1,range_y[2]+1),cex=1, pch=19)
    
      # CEM
      m <- Smooth1$teams[[j]]$indv[[i]]$mean
      s <-  sqrt(Smooth1$teams[[j]]$indv[[i]]$var  +  exp(paramList1$error[[1]][1] + paramList1$error[[1]][2]*time_i))
      lines(time_i,m + 2*s, col='red' )
      lines(time_i,m - 2*s, col='red' )
    
      # CEI
      m <- Smooth2$teams[[j]]$indv[[i]]$mean
      s <-  sqrt(Smooth2$teams[[j]]$indv[[i]]$var  + exp(paramList2$error[[1]]))
      lines(time_i,m + 2*s, col='blue' )
      lines(time_i,m - 2*s, col='blue' )
    
      # GP
      m <- Smooth3$teams[[j]]$indv[[i]]$mean
      s <-  sqrt(Smooth3$teams[[j]]$indv[[i]]$var  +  exp(paramList3$error[[1]][1]))
      lines(time_i,m + 2*s, col='green' )
      lines(time_i,m - 2*s, col='green' )
      }
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend('bottom',legend = c("CEM", "CEI", "GP"), col = c("red","blue", "green"), lwd = 5, xpd = TRUE, horiz = TRUE, cex = 1, seg.len=1, bty = 'n')
  }



