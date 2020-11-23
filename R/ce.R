### wrapping function for analysis


ce <- function(formula1, 
               formula2, 
               formula3, 
               emergence, 
               method,
               time = NULL,
               data) {
  
  data <- getData(formula1, 
                  formula2, 
                  formula3, 
                  emergence, 
                  method,
                  time,
                  data)
  
  n <- length(unlist(data[1]))
  
  model <- list(list(formula1,data$names[1]), 
                list(formula2,data$names[2]), 
                list(formula3,data$names[3]), 
                list(emergence,data$names[4]),
                method,
                time)
  
  names(model) <- c("Fixed effects", 
                    "Individual random effects", 
                    "Team random effects",
                    "Emergence model",
                    "Method",
                    "Time")
  
  
  
  object <- dataToObject(data)
  
  param <- param0(object)
  
  res <-optim(param, function(x){-loglik(x, object) })
  
  betas <-getbeta(res$par, object)
  
  cov_beta <- getCovbeta(res$par, object)
  
  paramList <- paramToList(res$par, object)
  
  loglik <- -res$value
  
  ce_model <- list(model, 
                   loglik,
                   betas, 
                   cov_beta, 
                   paramList,
                   n)
  
  names(ce_model) <- c("model",
                       "loglik",
                       "betas", 
                       "cov_beta", 
                       "covariances",
                       "n")
  
  # create class
    class(ce_model) <- "ce"

    
  return(ce_model)
}


