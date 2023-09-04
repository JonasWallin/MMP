predict.indv <- function(param, Obj, Obj.orig){
  
  beta <- getbeta(param, Obj.orig)
  paramList <- paramToList(param, Obj)
  
  
  
  res.pred <- matrix(0,
                       nrow  =length(Obj$indexTeams),
                       ncol  = 3)
  
  for(i in 1:length(Obj$teams)){
    
    Team_i <- Obj$teams[[i]]
    y_i    <- Obj$teams[[i]]$data$Y
    
    if(length(beta) > 0){
      X_i    <- Obj$teams[[i]]$data$X
      y_i    <- y_i - X_i%*%beta
    }
    n_i    <- length(y_i)
    
    ###
    # constructing covariance matrix of y
    ##
    ##########################
    # noise and X part
    ##########################
    Sigma_X    <- matrix(0, nrow = Team_i$n_obs,
                         ncol = Team_i$n_obs)
    Sigma_I    <- Sigma_X
    
    for(ii in 1:Team_i$nindv){
      #meauserment error setup
      sigmaI  <- matrix(0,
                        nrow = Team_i$indv[[ii]]$n,
                        ncol = Team_i$indv[[ii]]$n)
      for(iii in 1:length(Obj$indvCovs))
        sigmaI <- sigmaI + Obj$indvCovs[[iii]]$get_AtCA(paramList$indv[[iii]], Team_i$indv[[ii]])
      
      Sigma_I    <- Sigma_I + Team_i$indv[[ii]]$A%*%sigmaI%*%SparseM::t(Team_i$indv[[ii]]$A)
      
      
    }
    
    ##
    # team components
    ##
    for(ii in 1:length(Obj$teamCovs))
      Sigma_I <- Sigma_I + as.matrix(Obj$teamCovs[[ii]]$get_AtCA(paramList$team[[ii]], Team_i))
    
    ##
    # measuerment error
    ##
    Sigma_X <- Sigma_I
    for(iii in 1:length(Obj$errorCovs))
      Sigma_X <- Sigma_X + Obj$errorCovs[[iii]]$get_AtCA(paramList$error[[iii]], Team_i)
    
    index.na = as.vector(is.na(y_i)==F)
    y_smooth = Sigma_I[,index.na]%*%solve(Sigma_X[index.na,index.na], y_i[index.na])
    if(length(beta) > 0){
      X_i    <- Obj$teams[[i]]$data$X
      y_smooth    <- y_smooth + X_i%*%beta
    }
    TEMP    <- Sigma_I[,index.na]%*%Matrix::solve(Sigma_X[index.na,index.na], Sigma_I[index.na,])
    Sigma_X_cond <- Sigma_I - Sigma_I[,index.na]%*%Matrix::solve(Sigma_X[index.na,index.na], Sigma_I[index.na,])
    Sigma_XpN_cond <- Sigma_X - Sigma_I[,index.na]%*%Matrix::solve(Sigma_X[index.na,index.na], Sigma_I[index.na,])
    res.pred[,1] <- res.pred[,1] + as.vector(Team_i$A%*%y_smooth)
    res.pred[,2] <- res.pred[,2] + as.vector(Team_i$A%*%Matrix::diag(Sigma_X_cond))
    res.pred[,3] <- res.pred[,3] + as.vector(Team_i$A%*%Matrix::diag(Sigma_XpN_cond))
  }
  return(res.pred)
}


#'
#' Takes a ce.obj and new.data and creates a prediction around the position
#' 
#' @param  ce.object
#' @param  new.data [m x k] needs to have all relevant covariates
#' @return pred.vec [m x 3] return posterior mean, variance, variance + noise for
#'                          each row in new.data
#'
predict.ce <- function(ce.object, new.data, return.all = F){
  
  #add NA.data to new
  # 2023-09-04 (YB): this also creates NAs in the group variable
  # if it is not named "group", TODO: fix this
  data <- bind_rows(new.data,ce.object$model$data)
  pred.data <- getData(ce.object$model$`Fixed effects`[[1]], 
                  ce.object$model$`Individual random effects`[[1]], 
                  ce.object$model$`Team random effects`[[1]], 
                  ce.object$model$`Emergence model`[[1]], 
                  ce.object$model$Method,
                  data,
                  ce.object$model$`Team process`,
                  ce.object$model$Time)
  object <- dataToObject(pred.data)
 
  smoothRes <- predict.indv(ce.object$res$par, object,ce.object$object)
  if(return.all)
    return(smoothRes)
  
  return(smoothRes[1:dim(new.data)[1],])
}

