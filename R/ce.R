### wrapping function for analysis


ce <- function(formula1, 
               formula2, 
               formula3, 
               emergence, 
               method,
               method.team = NULL,
               time = NULL,
               data) {
  
  data <- getData(formula1, 
                  formula2, 
                  formula3, 
                  emergence, 
                  method,
                  data,
                  method.team,
                  time)
  
  n <- length(unlist(data[1]))
  
  model <- list(list(formula1,data$names[1]), 
                list(formula2,data$names[2]), 
                list(formula3,data$names[3]), 
                list(emergence,data$names[4]),
                method,
                method.team,
                time)
  
  names(model) <- c("Fixed effects", 
                    "Individual random effects", 
                    "Team random effects",
                    "Emergence model",
                    "Method",
                    "Team process",
                    "Time")
  
  
  
  object <- dataToObject(data)
  
  param <- param0(object)
  
  if (model$Method == "GP") {
  
  res <-optim(param, function(x){-loglik(x, object) })
  res <-optim(res$par, function(x){-loglik(x, object) })
  
    if (res$convergence != 0) {
      res <-optim(res$par, function(x){-loglik(x, object) }, method = "BFGS")
      if (res$convergence != 0) {
      res <-optim(res$par, function(x){-loglik(x, object) }, method = "BFGS")
      }
    }
  
  } else {
    res <-optim(param, function(x){-loglik(x, object) })
    if (res$convergence != 0) {
      res <-optim(res$par, function(x){-loglik(x, object) }, method = "BFGS")
      if (res$convergence != 0) {
        res <-optim(res$par, function(x){-loglik(x, object) }, method = "BFGS")
      }
    }
  }
  
  betas <-getbeta(res$par, object)
  
  cov_beta <- getCovbeta(res$par, object)
  
  paramList <- paramToList(res$par, object)
  
  paramPlot <- res$par
  
  loglik <- -res$value
  
  ce_model <- list(model, 
                   loglik,
                   betas, 
                   cov_beta, 
                   paramList,
                   n,
                   object,
                   paramPlot,
                   res)
  
  names(ce_model) <- c("model",
                       "loglik",
                       "betas", 
                       "cov_beta", 
                       "covariances",
                       "n",
                       "object",
                       "unlisted_covariances",
                       "res")
  
  # create class
    class(ce_model) <- "ce"

    
  return(ce_model)
}


