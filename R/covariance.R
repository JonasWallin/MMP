#'
#' covaraince processes for teamObject
#'
#'
#'
#'
library(R6)
library(nlme)
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
#'
#'
#' Xcov covariance structure obtanied through random effect i.e
#' The covariance of Y from the model
#' Y = X*Z where Z ~ N(0,\Sigma)
#'
Xcov <- R6Class("Xcov", list(
  d = NULL,
  initialize = function(d) {self$d <- d},
  get_name = function(){return('XCovSmooth')},
  get_param_length = function(){return(self$d*(self$d+1)/2 )},
  get_Cov = function(param,obj) {

    Sigma <- as.matrix(nlme::pdLogChol(param[1:length(param)]))
    return(Sigma)
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

#'
#' A X covaraince with expontial smoothing
#' The covariance of Y from the model
#' Y = diag(exp(-delta * d)) %*% X %*% Z where Z ~ N(0,\Sigma)
#'
#' @param obj - needs to contain X,d which is list of times
#'

XcovSmooth <- R6Class("XcovSmooth", list(
  d = NULL,
  initialize = function(d) {self$d <- d},
  get_name = function(){return('XCovSmooth')},
  get_param_length = function(){return(self$d*(self$d+1)/2 + 1)},
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

#'
#'
#' weighted covariance
#' The covariance of Y from the model
#' Z where Z ~ N(0, exp(E%*% delta )) where
#' W is n x d
#'
expWeightDiag <- R6Class("expWeightDiag", list(
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
  get_A  = function(param, obj){
    return(NULL)
  }
))

