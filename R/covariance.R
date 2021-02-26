#
# covaraince processes for teamObject
#
#
library(rSPDE)

#'OUcov
#'
#' Ornstein Uhlenbeck process
#'
#' dx_t = \theta (\mu - x_t) dt + \sigma dW_t
#'
#' @param d      - (n x n) distance matrix
#' @param param  - (2 x 1) \sigma, \theta,
#'                         \sigma - standard devation
#'                         \theta - 1/length_scale
OUcov <-function(d, param){

  sigma <- exp(param[1])
  theta <- exp(param[2])
  Cov   <- (sigma^2 /(2 * theta) ) * exp(-theta * d)
  return(Cov)
}

#'
#' Matern Covariance function
#'
#' @param d      - (n x n) distance matrix
#' @param param  - (3 x 1) \sigma, \kappa, \nu
#'                         \sigma - standard devation
#'                         \kappa - 1/length_scale
#'                         \nu    - differntiablity
#'                         
Materncov <- function(d, param){
  sigma <- exp(param[1])
  kappa <- exp(param[2])
  nu    <- exp(param[3])
  return(rSPDE::matern.covariance(d, 
                                  kappa = kappa, 
                                  nu = nu, 
                                  sigma = sigma))
}
#'Xcov
#'
#' Xcov covariance structure obtanied through random effect i.e
#' The covariance of Y from the model
#' Y = X*Z where Z ~ N(0,\Sigma)
#'
Xcov <- R6::R6Class("Xcov", list(
  d = NULL,
  initialize = function(d) {self$d <- d},
  get_name = function(){return('XCov')},
  get_param_length = function(){return(self$d*(self$d+1)/2 )},
  get_Cov = function(param,obj) {

    Sigma <- as.matrix(nlme::pdLogChol(param[1:length(param)]))
    return(Sigma)
  },
  get_mean = function(param, obj){
    m <- rep(0, self$d)
    return(m)
  },
  get_Amean = function(param, obj){
    m <- self$get_mean(param, obj)
    X <- obj$X
    return(X%*%m)
  },
  get_AtCA  = function(param, obj){
    X <- obj$X
    Sigma <- self$get_Cov(param, obj)
    A     <- self$get_A(param, obj)
    return(A%*%Sigma%*%t(A))
  },
  get_A  = function(param, obj){
    X <- obj$X
    return(X)
  }
))

#'XcovSmooth
#'
#' A X covaraince with expontial smoothing
#' The covariance of Y from the model
#' Y = diag(exp(W %*% delta)) %*% X %*% Z where Z ~ N(0,\Sigma)
#'
#' @param obj - needs to contain X,d which is list of times
#'

XcovSmooth <- R6::R6Class("XcovSmooth", list(
  d   = NULL,
  d_w = NULL,
  initialize = function(d, d_w) {self$d <- d;
                                 self$d_w <- d_w},
  get_name = function(){return('XCovSmooth')},
  get_param_length = function(){return(self$d*(self$d+1)/2 +  self$d_w)  },
  get_mean = function(param, obj){
    m <- rep(0, self$d)
    return(m)
  },
  get_Amean = function(param, obj){
    m <- self$get_mean(param, obj)
    X <- obj$X
    return(X%*%m)
  },
  get_Cov = function(param,obj) {

    Sigma <- as.matrix(nlme::pdLogChol(param[(self$d_w + 1):length(param)]))
    return(Sigma)
  },
  get_AtCA  = function(param, obj, cov_name='W'){
    X <- obj$X
    Sigma <- self$get_Cov(param, obj)
    A     <- self$get_A(param, obj, cov_name= cov_name)
    return(A%*%Sigma%*%t(A))
  },

  get_A  = function(param, obj, cov_name='W'){
    X <- obj$X
    DX <- diag(c(exp(obj[[cov_name]]%*%param[1:self$d_w])))%*%X
    return(DX)
  }
))

#'expWeightDiag
#'
#'
#' weighted covariance
#' The covariance of Y from the model
#' Z where Z ~ N(0, exp(E%*% delta )) where
#' W is n x d                                   ## *E is n x d?
#'
expWeightDiag <- R6::R6Class("expWeightDiag", list(
  d = NULL,
  initialize = function(d) {self$d <- d},
  get_name = function(){return('diagonal exp')},
  get_param_length = function(){return(self$d)},
  get_Cov = function(param,obj, cov_name = 'E') {

    Sigma <- diag(c(exp(obj[[cov_name]]%*%param)))
    return(Sigma)
  },
  get_AtCA  = function(param, obj, cov_name = 'E'){
    Sigma <- self$get_Cov(param, obj, cov_name)
    return(Sigma)
  },
  
  get_mean = function(param, obj){
    m <- rep(0, dim(obj[[cov_name]])[1])
    return(m)
  },
  get_A  = function(param, obj){
    return(NULL)
  }
))

#'
#' Covariance of an OU processes conditioning on
#' X(exp(D%*%\beta_delta))  = 0 where D%*%beta_delta is a parameter and
#' X(t) = 0 for t>delta=tmin + exp(D%*%beta_delta)
#' @param \delta 
#' @param \sigma
#' @param \range


OUbridge <- R6::R6Class("OUbridge", list(
  d = NULL,
  tmin = NULL,
  tmax = NULL,
  initialize = function(tmin, tmax,d_D) {self$tmin <- tmin;
                                     self$tmax <- tmax;
                                     self$d<- 2 + d_D },
  get_name = function(){return('OU bridge')},
  get_param_length = function(){return(self$d)},
  get_Cov = function(param, obj, cov_name = 'D', time = 'time') {

    delta <- self$tmin + c(exp(obj[[cov_name]]%*%as.vector(param[1:(self$d-2)])))

    Sigma <- matrix(0, 
                    nrow = length(obj[[time]]),
                    ncol = length(obj[[time]]))
    # time < delta
    
    less_delta <- obj[[time]] < delta
    time <- c(obj[[time]][less_delta], delta) 
    Dist <- as.matrix(dist(time))    
    n_ <- length(time)
   
    Sigma_p <- OUcov(Dist, param[(self$d-1):self$d]) # the two last param
    # conditional distribution
    Sigma[less_delta,less_delta] <- Sigma_p[1:(n_-1),1:(n_-1)] -
                                    Sigma_p[1:(n_-1),n_, drop = FALSE]%*%t(Sigma_p[1:(n_-1),n_, drop = FALSE]/Sigma_p[n_,n_])
  return(Sigma)
  },
  get_AtCA  = function(param, obj,cov_name ='D', time = 'time'){
    Sigma <- self$get_Cov(param, obj, cov_name = cov_name, time = time)
    return(Sigma)
  },
  
  get_mean = function(param, obj, time = 'time'){
    #TODO
    m <- rep(0, length(obj[[time]]))
    return(m)
  },
  
  get_Amean = function(param, obj, time = 'time'){
    return(self$get_mean(param, obj, time = 'time'))
  },
  get_A  = function(param, obj){
    return(NULL)
  }
))

MaternBridge <- R6::R6Class("Matern bridge", list(
  d = NULL,
  tmin = NULL,
  tmax = NULL,
  initialize = function(tmin, tmax,d_D) {self$tmin <- tmin;
  self$tmax <- tmax;
  self$d<- 3 + d_D },
  get_name = function(){return('Matern bridge')},
  get_param_length = function(){return(self$d)},
  get_Cov = function(param, obj, cov_name = 'D', time = 'time') {
    
    delta <- self$tmin + c(exp(obj[[cov_name]]%*%as.vector(param[1:(self$d-3)])))
    Sigma <- matrix(0, 
                    nrow = length(obj[[time]]),
                    ncol = length(obj[[time]]))
    # time < delta
    
    less_delta <- obj[[time]] < delta
    time <- c(obj[[time]][less_delta], delta)
    D <- as.matrix(dist(time))
    n_ <- length(time)
    
    Sigma_p <- Materncov(D, param[(self$d-2):self$d])
    # conditional distribution
    Sigma[less_delta,less_delta] <- Sigma_p[1:(n_-1),1:(n_-1)] -
      Sigma_p[1:(n_-1),n_, drop = FALSE]%*%t(Sigma_p[1:(n_-1),n_, drop = FALSE]/Sigma_p[n_,n_])
    return(Sigma)
  },
  get_AtCA  = function(param, obj,cov_name ='D', time = 'time'){
    Sigma <- self$get_Cov(param, obj, cov_name = cov_name, time = time)
    return(Sigma)
  },
  
  get_mean = function(param, obj, time = 'time'){
    #TODO
    m <- rep(0, length(obj[[time]]))
    return(m)
  },
  
  get_Amean = function(param, obj, time = 'time'){
    return(self$get_mean(param, obj, time = 'time'))
  },
  get_A  = function(param, obj){
    return(NULL)
  }
))

