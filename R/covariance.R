#
# covaraince processes for teamObject
#
#


#'OUcov
#'
#' Ornstein Uhlenbeck process
#'
#' dx_t = \theta (\mu - x_t) dt + \sigma dW_t
#'
#' @param d      - (n x n) distance matrix
#' @param param  - (2 x 1) \sigma, \theta,
OUcov <-function(d, param){

  sigma <- exp(param[1])
  theta <- exp(param[2])
  Cov   <- (sigma^2 /(2 * theta) ) * exp(-theta * d)
  return(Cov)
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
  get_name = function(){return('XCovSmooth')},
  get_param_length = function(){return(self$d*(self$d+1)/2 )},
  get_Cov = function(param,obj) {

    Sigma <- as.matrix(nlme::pdLogChol(param[1:length(param)]))
    return(Sigma)
  },
  get_mean = function(param, obj){
    m <- rep(0, self$d)
    return(m)
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
#' Y = diag(exp(-delta * d)) %*% X %*% Z where Z ~ N(0,\Sigma)
#'
#' @param obj - needs to contain X,d which is list of times
#'

XcovSmooth <- R6::R6Class("XcovSmooth", list(
  d = NULL,
  initialize = function(d) {self$d <- d},
  get_name = function(){return('XCovSmooth')},
  get_param_length = function(){return(self$d*(self$d+1)/2 + 1)},
  get_mean = function(param, obj){
    m <- rep(0, self$d)
    return(m)
  },
  get_Cov = function(param,obj) {

    Sigma <- as.matrix(nlme::pdLogChol(param[2:length(param)]))
    return(Sigma)
  },
  get_AtCA  = function(param, obj){
    X <- obj$X
    d <- obj$w
    Sigma <- self$get_Cov(param, obj)
    A     <- self$get_A(param, obj)
    return(A%*%Sigma%*%t(A))
  },

  get_A  = function(param, obj){
    X <- obj$X
    d <- obj$w
    DX <- diag(c(exp(-d*param[1])))%*%X
    return(DX)
  }
))

#'expWeightDiag
#'
#'
#' weighted covariance
#' The covariance of Y from the model
#' Z where Z ~ N(0, exp(E%*% delta )) where
#' W is n x d
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
#' X(\delta)  = 0 where \delta is a parameter and
#' X(t) = 0 for t>\delta
#' @param \delta in [tmin, \inf)
#' @param \sigma
#' @param \range


OUbridge <- R6::R6Class("OUbridge", list(
  d = 3,
  tmin = NULL,
  tmax = NULL,
  initialize = function(tmin, tmax) {self$tmin <- tmin
                                     self$tmax <- tmax},
  get_name = function(){return('OU bridge')},
  get_param_length = function(){return(self$d)},
  get_Cov = function(param, obj, cov_name = 'time') {
    
    delta <- self$tmin + exp(param[1])
    Sigma <- matrix(0, 
                    nrow = length(obj[[cov_name]]),
                    ncol = length(obj[[cov_name]]))
    # time < delta
    less_delta <- obj[[cov_name]] < delta
    time <- c(obj[[cov_name]][less_delta], delta)
    
    D <- as.matrix(dist(time))
    n_ <- length(time)
    Sigma_p <- OUcov(D, param[2:3])
    # conditional distribution
    Sigma[less_delta,less_delta] <- Sigma_p[1:(n_-1),1:(n_-1)] -
                                    Sigma_p[1:(n_-1),n_, drop = FALSE]%*%t(Sigma_p[1:(n_-1),n_, drop = FALSE]/Sigma_p[n_,n_])
    return(Sigma)
  },
  get_AtCA  = function(param, obj, cov_name = 'time'){
    Sigma <- self$get_Cov(param, obj, cov_name)
    return(Sigma)
  },
  
  get_mean = function(param, obj, cov_name = 'time'){
    #TODO
    m <- rep(0, length(obj[[cov_name]]))
    return(m)
  },
  get_A  = function(param, obj){
    return(NULL)
  }
))
