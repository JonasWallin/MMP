### wrapping function for analysis


ce <- function(formula1, 
               formula2, 
               formula3, 
               emergence, 
               method,
               method.team = NULL,
               time = NULL,
               data,
               REML = F) {
  
  data <- getData(formula1, 
                  formula2, 
                  formula3, 
                  emergence, 
                  method,
                  data,
                  method.team,
                  time)
  
  n <- length(unlist(data[1]))
  
  model <- list('Fixed effects'             = list(formula1,data$names[1]), 
                'Individual random effects' = list(formula2,data$names[2]), 
                'Team random effects'       = list(formula3,data$names[3]), 
                'Emergence model'           = list(emergence,data$names[4]),
                'Method'                    = method,
                'Team process'              = method.team,
                'Time'                      = time)
  
  
  object <- dataToObject(data)
  
  param <- param0(object)
  
  if (model$Method == "GP") {
  
  res <-optim(param, function(x){-loglik(x, object,REML ) })
  res <-optim(res$par, function(x){-loglik(x, object, REML) })
  
    if (res$convergence != 0) {
      res <-optim(res$par, function(x){-loglik(x, object, REML) },control = list(maxit=5000))
      if (res$convergence != 0) {
      res <-optim(res$par, function(x){-loglik(x, object, REML) })
      }
    }
  
  } else {
    res <-optim(param, function(x){-loglik(x, object,REML) })
    res <-optim(res$par, function(x){-loglik(x, object, REML) },control = list(maxit=5000))
    if (res$convergence != 0) {
      res <-optim(res$par, function(x){-loglik(x, object, REML) })
      if (res$convergence != 0) {
        res <-optim(res$par, function(x){-loglik(x, object, REML) })
      }
    }
  }
  
  betas <-getbeta(res$par, object)
  
  cov_beta <- getCovbeta(res$par, object)
  
  paramList <- paramToList(res$par, object)
  
  paramPlot <- res$par
  
  loglik <- loglik(res$par, object)
  
  # AIC: 2*k - 2*loglik, k = number of estimated parameters 
  k <- sum(length(paramPlot),
           length(betas))
  
  aic <- (2*k)-(2*loglik)
  
  ce_model <- list(model, 
                   loglik,
                   betas, 
                   cov_beta, 
                   paramList,
                   n,
                   object,
                   paramPlot,
                   res,
                   aic)
  
  names(ce_model) <- c("model",
                       "loglik",
                       "betas", 
                       "cov_beta", 
                       "covariances",
                       "n",
                       "object",
                       "unlisted_covariances",
                       "res",
                       "AIC")
  
  # create class
    class(ce_model) <- "ce"

    
  return(ce_model)
}


