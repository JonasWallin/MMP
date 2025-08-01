
#' loglikelihood new
#' @param - the parameters of the model
#' @param - object type
#'          $teams - list of team data
#'                   $data$y
#'
#' @return - the logliklihood
loglik <- function(param, Obj, REML=F){

  return(likAndBeta(param, Obj, REML)$lik)
}
#'
#' Get the fixed effect for the model
#'
#'
getbeta <- function(param, Obj){
  return(likAndBeta(param, Obj)$beta)
}

#'
#' Get the covariance matrix of fixed effect for the model
#' hold variance parameters fixed
#'
getCovbeta <- function(param, Obj){
  return(solve(likAndBeta(param, Obj)$Hbeta))
}

#'
#' Computes the log-likelihood (or Restricted)
#' for variance compoents and also the ML-estimated fixed effect
#' and resulting Hessian
#' 
likAndBeta  <- function(param, Obj, REML=FALSE){

  paramList <- paramToList(param, Obj)
  lik <- 0

  if(is.null(Obj$X) ==F){
    XtSigmainvX <- matrix(0,
                          nrow = dim(Obj$X)[2],
                          ncol = dim(Obj$X)[2]
    )
    ySigmainvX <- rep(0, dim(Obj$X)[2])

  }
  if(is.null(Obj$X) ==F)
    XtSigmainvX <- matrix(0, nrow= dim(Obj$X)[2],ncol= dim(Obj$X)[2])
  
  for(i in 1:length(Obj$teams)){
    Team_i <- Obj$teams[[i]]
    y_i    <- as.matrix(Obj$teams[[i]]$data$Y)
    if(is.null(Obj$X) ==F)
      X_i    <- Obj$teams[[i]]$data$X
    n_i    <- length(y_i)
    ###
    # constructing covariance matrix of y
    ##
    ##########################
    # noise and X part
    ##########################
    Sigma_X    <- matrix(0, nrow = Team_i$n_obs,
                         ncol = Team_i$n_obs)
    for(ii in 1:Team_i$nindv){
      #meauserment error setup
      sigmaI  <- matrix(0,
                        nrow = Team_i$indv[[ii]]$n,
                        ncol = Team_i$indv[[ii]]$n)
      meanI   <- rep(0,  Team_i$indv[[ii]]$n)
      for(iii in 1:length(Obj$indvCovs)){
        sigmaI <- sigmaI + Obj$indvCovs[[iii]]$get_AtCA(paramList$indv[[iii]], Team_i$indv[[ii]])
        meanI  <- meanI  + Obj$indvCovs[[iii]]$get_Amean(paramList$indv[[iii]], Team_i$indv[[ii]])
      }
      index_ <- Team_i$indv[[ii]]$A_list
      #Sigma_X    <- Sigma_X + Team_i$indv[[ii]]$A%*%sigmaI%*%SparseM::t(Team_i$indv[[ii]]$A)
      Sigma_X[index_,index_]    <- Sigma_X[index_,index_] + sigmaI
      y_i[index_]        <- y_i[index_] - as.vector(meanI)
    }
    #building the blockmatrix
    #Sigma_X    <- Matrix::.bdiag(Sigma_X)

    #for(iii in 1:length(Obj$errorCovs))
    #  Sigma_X <- Sigma_X + Obj$errorCovs[[iii]]$get_AtCA(paramList$error[[iii]], Team_i)
    Sigma_X_d <- rep(0, dim(Sigma_X)[1])
    for(iii in 1:length(Obj$errorCovs))
      Sigma_X_d <- Sigma_X_d + Obj$errorCovs[[iii]]$add_diagonal(paramList$error[[iii]], Team_i)
    diag(Sigma_X) <- diag(Sigma_X) + Sigma_X_d
    ##
    # team components
    ##
    for(ii in 1:length(Obj$teamCovs))
      Sigma_X <- Sigma_X + as.matrix(Obj$teamCovs[[ii]]$get_AtCA(paramList$team[[ii]], Team_i))

    Sigma_X <- as.matrix(Sigma_X)
    L = tryCatch({L <- chol(Sigma_X)}, error = function(err){return(-Inf)})
    if(is.null(dim(L)))
      return(list(lik =-Inf))
    # determinant of covariance matrix
    lik <- lik - sum(log(diag(L)))
    ##
    # computing likelihood
    # lik = -0.5 * y^T (Sigma^-1 y)
    Sigma12invY <- solve(t(L),y_i)
    lik <- lik -0.5 * as.matrix(t(Sigma12invY)%*%Sigma12invY)


    # add beta component!!!
    # and REML
    #SigmaInvX = solve(Sigma_X, X)
    if(is.null(Obj$X) ==F){
      SigmainvX   <- forwardsolve(L, backsolve(L, X_i, transpose = TRUE), upper.tri = TRUE)
      XtSigmainvX <- XtSigmainvX + t(X_i)%*%SigmainvX
      ySigmainvX  <- ySigmainvX  + t(y_i)%*%SigmainvX
    }


  }

  beta = c()
  
  #-0.5 * beta_hat^T \sigma_hat^-1 beta_hat + \sum_i y_i \sigma_hat^-1 beta_hat
  # = -0.5 (D%*%C)^T D^{-1} D%*%C + C^T * D^{-1} (D%*%C)
  # 
  if(is.null(Obj$X) ==F){
    beta = solve(XtSigmainvX, t(ySigmainvX))
    lik <- lik + 0.5 * ySigmainvX%*%beta
    
    #use restitriced maximum likelihood
    if(REML==T){
      L <- chol(XtSigmainvX)
      # determinant of covariance matrix
      lik <- lik + sum(log(diag(L)))
      
    }
  }
  return(list(lik = lik, beta=beta, Hbeta = XtSigmainvX))
}
#gives the smoothing distribution of the indivual
# in form of mean and covariance
smoothIndivual <- function(param, Obj){

 beta <- getbeta(param, Obj)
 paramList <- paramToList(param, Obj)

 smoothRes <- list()
 smoothRes$teams <- list()
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

    y_smooth = Sigma_I%*%solve(Sigma_X, y_i)
    smoothRes$teams[[i]] <- list()
    smoothRes$teams[[i]]$indv <- list()
   for(ii in 1:Team_i$nindv){
     
     sigmaI <- Matrix::t(Team_i$indv[[ii]]$A)%*%Sigma_I%*%Team_i$indv[[ii]]$A
     
     #Indiv | obs
     mu_I    <- Matrix::t(Team_i$indv[[ii]]$A)%*%y_smooth
     if(length(beta) > 0){
       X_i    <- Obj$teams[[i]]$data$X
       mu_I    <- as.vector(mu_I) + Matrix::t(Team_i$indv[[ii]]$A)%*%X_i%*%beta
     }
     sigma_I <- sigmaI - sigmaI%*%Matrix::solve(Matrix::t(Team_i$indv[[ii]]$A)%*%Sigma_X%*%Team_i$indv[[ii]]$A,sigmaI)
     smoothRes$teams[[i]]$indv[[ii]] <- data.frame(mean = as.vector(mu_I),
                                                 var  = Matrix::diag(sigma_I))
   }


 }
 return(smoothRes)
}

