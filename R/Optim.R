
safe_loglik_term <- function(y, Sigma) {
  tryCatch({
    qf <- as.numeric(crossprod(y, solve(Sigma, y)))  
    det_log <- determinant(Sigma, logarithm = TRUE)
    if (det_log$sign != 1) return(-Inf)              # non PD or det<=0
    -0.5 * (qf + as.numeric(det_log$modulus))        # (+ const if needed)
  }, error = function(e) -Inf)
}

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
    val <- tryCatch({
      Sigma12invY <- solve(t(L), y_i)
      -0.5 * sum(Sigma12invY^2)
    }, error = function(e) -Inf)
    if(val == -Inf)
      return(list(lik =-Inf))
    lik <- lik + val


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
    
    tc <- tryCatch({
      beta <- solve(XtSigmainvX, t(ySigmainvX))
      
      # Guard against NaN/Inf in beta
      if (any(!is.finite(beta))) stop("beta contains non-finite entries")
      
      val <- 0.5 * as.numeric(ySigmainvX %*% beta)
      
      # Guard against NaN/Inf in val
      if (!is.finite(val)) stop("val is non-finite")
      
      list(beta = beta, val = val)
    }, error = function(e) {
      list(beta = NULL, val = -Inf)
    })
    
    if (!is.finite(tc$val)) {
      return(list(lik = -Inf))
    }
    
    beta <- tc$beta
    lik  <- lik + tc$val

    
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

#' Estimate group effects via maximum likelihood
#'
#' This function wraps `likelihood.groupavg` and performs optimization
#' of group-level parameters given a model object.
#'
#' @param Obj A model object with fields `teamCovs`, `teams`, and
#'   methods to update parameters.
#' @param init_param Numeric vector of initial parameter guesses
#'   (first element is global mean mu, remainder are covariance params).
#' @param control Optional list of control parameters passed to `optim`.
#' @return A list with elements:
#'   * `fit`: the final `optim` result,
#'   * `param_list`: a list of updated team parameters,
#'   * `param_vector`: the concatenated parameter vector returned by `listToParam`.
estimate_group_effects <- function(param_list, Obj, init_param, control = list(maxit = 1000, trace = FALSE)) {
  # Negative log-likelihood wrapper
  nll <- function(p) -likelihood.groupavg(p, Obj$object)
  
  # First optimization
  fit1 <- optim(par     = init_param,
                fn      = nll,
                method  = "BFGS",
                control = control)
  
  # Refine around the solution
  fit2 <- optim(par     = fit1$par,
                fn      = nll,
                method  = "BFGS",
                control = control)
  
  # Update covariance parameter list from the fitted pars (excluding mu)
  raw_cov_pars <- fit2$par[-1]
  param_list   <- updateTeamParamList(param_list, Obj$object, raw_cov_pars)
  
  
  return(param_list)
}

#' Estimate error and individual-level effects via maximum likelihood
#'
#' This function wraps `likelihood.degroup` and performs optimization
#' of error and individual-level parameters given a model object.
#'
#' @param Obj A model object with fields `errorCovs`, `indvCovs`, `teams`.
#' @param init_param Numeric vector of initial parameter guesses.
#' @param control Optional list of control parameters passed to `optim`.
#' @return A list with elements:
#'   * `fit`: the final `optim` result,
#'   * `param_list`: a list of updated error and individual parameters.
estimate_error_indv_effects <- function(param_list, Obj, init_param, control = list(maxit = 1000, trace = FALSE)) {
  # Negative log-likelihood wrapper
  nll <- function(p) -likelihood.degroup(p, Obj$object)
  
  # First optimization
  fit1 <- optim(par     = init_param,
                fn      = nll,
                method  = "BFGS",
                control = control)
  
  # Refine around the solution
  fit2 <- optim(par     = fit1$par,
                fn      = nll,
                method  = "BFGS",
                control = control)
  
  # Update model covariances with the fitted parameters
  param_list <-  updateErrorIndvParamList(param_list, Obj$object, fit2$par)
  return(param_list)
}


likelihood.groupavg <- function(param, Obj){
  loglik <- 0
  paramList <- list()
  paramList$team <- list()
  mu <- param[1]
  param <- param[-1]
  for (i in 1:length(Obj$teamCovs)) {
    covObj <- Obj$teamCovs[[i]]
    nParObj <- covObj$get_param_length()
    paramList$team[[i]] <- param[1:nParObj]
    param <- param[-(1:nParObj)]
  }
  for(i in 1:length(Obj$teams)){
    y.avg.pos <- tibble(
      time = Obj$teams[[i]]$time,
      y_i  = Obj$teams[[i]]$data$Y,
      idx  = seq_along(time)        # original position
    ) %>%
      dplyr::group_by(time) %>%
      dplyr::summarise(
        first_idx = first(idx),     # grabs the smallest original index
        y_avg     = mean(y_i),      # group mean
        .groups   = "drop"
      )
    Sigma_X    <- matrix(0, nrow = Obj$teams[[i]]$n_obs,
                         ncol = Obj$teams[[i]]$n_obs)
    for(ii in 1:length(Obj$teamCovs)){
      Sigma_X <- Sigma_X + as.matrix(Obj$teamCovs[[ii]]$get_AtCA(paramList$team[[ii]], Obj$teams[[i]]))
    }
    Sigma_X <- Sigma_X[y.avg.pos$first_idx,y.avg.pos$first_idx]
    y_i <- y.avg.pos$y_avg
    # calculate the likelihood of y_i which is zero mean and covariance Sigma_X
    # check if Sigma_X si positive definite
 
    L = tryCatch({L <- chol(Sigma_X)}, error = function(err){return(Inf)})
    Ly = tryCatch({Ly <- solve(t(L),y_i-mu)}, error = function(err){return(Inf)})
    if(is.null(dim(L)))
      return(Inf)
    loglik <- loglik -0.5 * t(Ly) %*% Ly - sum(log(diag(L)))
  }
  
  return(loglik)
}


likelihood.degroup <- function(param, Obj){
  loglik <- 0
  paramList <- list()
  paramList$error <- list()
  for(i in 1:length(Obj$errorCovs)){
    covObj <- Obj$errorCovs[[i]]
    nParObj <- covObj$get_param_length()
    paramList$error[[i]] <- param[1:nParObj]
    param <- param[-(1:nParObj)]
  }
  
  if(length(Obj$indvCovs)>0){
    paramList$indv <- list()
    for(i in 1:length(Obj$indvCovs)){
      covObj <- Obj$indvCovs[[i]]
      nParObj <- covObj$get_param_length()
      paramList$indv[[i]] <- param[1:nParObj]
      param <- param[-(1:nParObj)]
    }
  }
  
  for(i in 1:length(Obj$teams)){
    Team_i <- Obj$teams[[i]]
    y.residuals <- tibble(
      time = Team_i$time,
      y_i  = Team_i$data$Y,
      idx  = seq_along(time)
    ) %>%
      dplyr::group_by(time) %>%
      dplyr::mutate(
        # subtract the group mean from each observation
        y_centered = y_i - mean(y_i)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(idx) 
    Sigma_X    <- matrix(0, nrow = Team_i$n_obs,
                         ncol = Team_i$n_obs)
    y_i <- y.residuals$y_centered
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
    Sigma_X_d <- rep(0, dim(Sigma_X)[1])
    for(iii in 1:length(Obj$errorCovs))
      Sigma_X_d <- Sigma_X_d + Obj$errorCovs[[iii]]$add_diagonal(paramList$error[[iii]], Team_i)
    diag(Sigma_X) <- diag(Sigma_X) + Sigma_X_d
    
    # calculate the likelihood of y_i which is zero mean and covariance Sigma_X
    # check if Sigma_X si positive definite
  #  print(Sigma_X)
  #  print(eigen(Sigma_X)$values)
  #  if (any(eigen(Sigma_X)$values <= 0)) {
 #     return(Inf)
  #  }
    
    
    ## use it like:
    loglik <- loglik + safe_loglik_term(y_i, Sigma_X)
 #   print(loglik)
    if (!is.finite(loglik)) return(-Inf)
  #  if(abs(loglik) == Inf)
  #    return(-Inf)
  }
  
  return(loglik)
}

#' Estimate mu, error, and team effects via maximum likelihood
#'
#' Jointly estimates the global mean \eqn{\mu}, **error**, and **team**
#' covariance parameters by maximizing the likelihood of time-averaged data.
#' The parameter vector must be ordered as:
#' \code{c(mu, error blocks..., team blocks...)}.
#'
#' @param param_list A parameter list as returned by \code{paramToList()}.
#' @param Obj A wrapper with field \code{object} holding \code{errorCovs},
#'   \code{teamCovs}, and \code{teams}.
#' @param init_param Numeric vector of initial guesses concatenated as
#'   \code{c(mu, error blocks..., team blocks...)}.
#' @param control List passed to \code{optim}; default \code{list(maxit=1000, trace=FALSE)}.
#'
#' @return The updated \code{param_list} with new \code{error} and \code{team}
#'   entries (other components unchanged). The MLE of \code{mu} is not stored in
#'   \code{param_list}; retrieve it from the optimizer if needed.
#'
#' @examples
#' # init <- c(mu0, err_pars..., team_pars...)
#' # param_list <- estimate_error_team_effects(param_list, Obj, init)
#' @export
estimate_error_team_effects <- function(param_list,
                                        Obj,
                                        init_param,
                                        control = list(maxit = 1000, trace = FALSE)) {
  nll <- function(p) -likelihood.errteam_avg(p, Obj$object)
  
  fit1 <- optim(par = init_param, fn = nll, control = control)
  fit2 <- optim(par = fit1$par, fn = nll, method = "BFGS", control = control)
  
  # Write back ONLY error+team (exclude the first element mu)
  fitted_mu <- fit2$par[1]
  fitted_rest <- fit2$par[-1]
  
  param_list <- updateErrorTeamParamList(param_list, Obj$object, fitted_rest)
  
  # If you want to keep mu accessible without changing structures elsewhere:
  # attr(param_list, "mu_hat") <- fitted_mu
  
  return(param_list)
}

# Internal: joint log-likelihood for c(mu, error..., team...) on time-averaged data
likelihood.errteam_avg <- function(param, Obj) {
  if (length(Obj$errorCovs) == 0L || length(Obj$teamCovs) == 0L) {
    stop("likelihood.errteam_avg requires both errorCovs and teamCovs.")
  }
  
  # Split param -> mu | error blocks | team blocks
  mu <- param[1]
  rest <- param[-1]
  
  # Extract error
  paramList <- list(error = vector("list", length(Obj$errorCovs)),
                    team  = vector("list", length(Obj$teamCovs)))
  idx <- 1L
  for (i in seq_along(Obj$errorCovs)) {
    n <- as.integer(Obj$errorCovs[[i]]$get_param_length())
    paramList$error[[i]] <- rest[idx:(idx + n - 1L)]
    idx <- idx + n
  }
  # Extract team
  for (i in seq_along(Obj$teamCovs)) {
    n <- as.integer(Obj$teamCovs[[i]]$get_param_length())
    paramList$team[[i]] <- rest[idx:(idx + n - 1L)]
    idx <- idx + n
  }
  
  loglik <- 0
  
  for (i in seq_along(Obj$teams)) {
    Team_i <- Obj$teams[[i]]
    n_obs  <- Team_i$n_obs
    y_full <- Team_i$data$Y
    times  <- Team_i$time
    
    # Build index sets per unique time, preserving first-occurrence order
    idx_list_raw <- split(seq_len(n_obs), times)
    first_pos <- vapply(idx_list_raw, function(id) min(id), numeric(1))
    ord <- order(first_pos)
    idx_list <- idx_list_raw[ord]
    m <- length(idx_list)
    
    # Time-averaged observations and mean vector
    y_bar <- vapply(idx_list, function(id) mean(y_full[id]), numeric(1))
    mu_vec <- rep(mu, m)
    
    # Aggregation matrix M (m x n_obs) with rows summing to 1 within time
    M <- matrix(0, nrow = m, ncol = n_obs)
    for (t in seq_len(m)) {
      id <- idx_list[[t]]
      M[t, id] <- 1 / length(id)
    }
    
    # Team covariance on full grid, then aggregate: M * Sigma_team * M^T
    Sigma_team_full <- matrix(0, nrow = n_obs, ncol = n_obs)
    for (k in seq_along(Obj$teamCovs)) {
      Sigma_team_full <- Sigma_team_full +
        as.matrix(Obj$teamCovs[[k]]$get_AtCA(paramList$team[[k]], Team_i))
    }
    Sigma_team_bar <- M %*% Sigma_team_full %*% t(M)
    
    # Error variances on full grid, then aggregate diagonals with 1/n_t^2
    d_full <- numeric(n_obs)
    for (k in seq_along(Obj$errorCovs)) {
      d_full <- d_full + Obj$errorCovs[[k]]$add_diagonal(paramList$error[[k]], Team_i)
    }
    d_bar <- vapply(idx_list, function(id) sum(d_full[id]) / (length(id)^2), numeric(1))
    
    # Final covariance for y_bar
    Sigma <- Sigma_team_bar
    diag(Sigma) <- diag(Sigma) + d_bar
    
    # PD check and log-likelihood via Cholesky (like your other functions)
    L <- tryCatch(chol(Sigma), error = function(e) NULL)
    if (is.null(L)) return(-Inf)  # => nll = Inf
    
    Ly <- tryCatch(solve(t(L), y_bar - mu_vec), error = function(e) NULL)
    if (is.null(Ly)) return(-Inf)
    
    loglik <- loglik - 0.5 * drop(t(Ly) %*% Ly) - sum(log(diag(L)))
  }
  
  return(loglik)
}
