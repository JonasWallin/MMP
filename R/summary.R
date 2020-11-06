### Summary function


summary.ce <- function(object) {
  
  model <- as.character(object$model[5])
  
  if (model == "CEM") {
  ## fixed effects
  fe <- as.character(object$model$`Fixed effects`[1]) # formula
  fe_names <- unlist(object$model$`Fixed effects`[2]) # names of covariates
  fe_param <- object$betas
  
  fe_mat <- cbind(fe_param,diag(t(object$cov_beta)))
  dimnames(fe_mat) <- list(fe_names,c("Estimate", "Variance"))
  
  ## random effects
  # individual
  
  ire <- as.character(object$model$`Individual random effects`[1]) # formula
  ire_names <- unlist(object$model$`Individual random effects`[2]) # names of covariates
  
  
  ire_param <- unlist(object$covariances$indv)
  ire_mat <- as.matrix(nlme::pdLogChol(ire_param))
  dimnames(ire_mat) <- list(ire_names,ire_names)
  

  # team
  
  tre <- as.character(object$model$`Team random effects`[1]) # formula
  tre_names <- unlist(object$model$`Team random effects`[2]) # names of covariates
  
  tre_param <- unlist(object$covariances$team)
  tre_mat <- as.matrix(nlme::pdLogChol(tre_param))
  dimnames(tre_mat) <- list(tre_names,tre_names)
  
  
  # emergence
  em <- as.character(object$model$`Emergence model`[1]) # formula
  em_names <- unlist(object$model$`Emergence model`[2]) # names of covariates
  
  em_param <- unlist(object$covariances$error)
  
  sigma2 <- exp(em_param[1])
  delta_param <- em_param[-1]/2
    
  em_mat <- matrix(c(sigma2, delta_param), 
                     dimnames = list(c("sigma^2",em_names[-1]),"Estimate"))
  }
  
  else if (model == "CEI") {
    
    ## fixed effects
    fe <- as.character(object$model$`Fixed effects`[1]) # formula
    fe_names <- unlist(object$model$`Fixed effects`[2]) # names of covariates
    fe_param <- object$betas
    
    fe_mat <- cbind(fe_param,diag(t(object$cov_beta)))
    dimnames(fe_mat) <- list(fe_names,c("Estimate", "Variance"))
    
    fe_measurement_error <- exp(unlist(object$covariances$error)) # do we want to print this?
    
    ## random effects
    
    # individual (exclude for now, does not work to have covariates here)
    
    #ire <- as.character(object$model$`Individual random effects`[1]) # formula
    #ire_names <- unlist(object$model$`Individual random effects`[2]) # names of covariates
    
    
    #ire_param <- unlist(object$covariances$indv)
    #ire_mat <- as.matrix(nlme::pdLogChol(ire_param))
    #dimnames(ire_mat) <- list(ire_names,ire_names)
    
    # team
    
    tre <- as.character(object$model$`Team random effects`[1]) # formula
    tre_names <- unlist(object$model$`Team random effects`[2]) # names of covariates
    
    tre_param <- unlist(object$covariances$team)
    tre_mat <- as.matrix(nlme::pdLogChol(tre_param))
    dimnames(tre_mat) <- list(tre_names,tre_names)
    
    
    # emergence (on individual level for CEI)
    em <- as.character(object$model[4])
    
    em <- as.character(object$model$`Emergence model`[1]) # formula
    em_names <- unlist(object$model$`Emergence model`[2]) # names of covariates
    
    em_param <- unlist(object$covariances$indv)
    
    sigma2 <- exp(em_param[1])
    delta_param <- em_param[-1]/2
    
    em_mat <- matrix(c(sigma2, delta_param), 
                     dimnames = list(c("sigma^2",em_names[-1]),"Estimate"))
  }
  
  else if (model == "GP") { # not sure what we want to show for this one
    
    ## fixed effects
    fe <- as.character(object$model$`Fixed effects`[1]) # formula
    fe_names <- unlist(object$model$`Fixed effects`[2]) # names of covariates
    fe_param <- object$betas
    
    fe_mat <- cbind(fe_param,diag(t(object$cov_beta)))
    dimnames(fe_mat) <- list(fe_names,c("Value", "Variance"))
    
    fe_measurement_error <- exp(unlist(object$covariances$error)) # do we want to print this?
    
    ## random effects
    
    # individual (exclude for now, does not work to have covariates here)
    
    #ire <- as.character(object$model$`Individual random effects`[1]) # formula
    #ire_names <- unlist(object$model$`Individual random effects`[2]) # names of covariates
    
    
    #ire_param <- unlist(object$covariances$indv)
    #ire_mat <- as.matrix(nlme::pdLogChol(ire_param))
    #dimnames(ire_mat) <- list(ire_names,ire_names)
    
    # team
    
    tre <- as.character(object$model$`Team random effects`[1]) # formula
    tre_names <- unlist(object$model$`Team random effects`[2]) # names of covariates
    
    tre_param <- unlist(object$covariances$team)
    tre_mat <- as.matrix(nlme::pdLogChol(tre_param))
    dimnames(tre_mat) <- list(tre_names,tre_names)
    
    
    # emergence (on individual level for CEI)
    em <- as.character(object$model[4])
    
    em <- as.character(object$model$`Emergence model`[1]) # formula
    em_names <- unlist(object$model$`Emergence model`[2]) # names of covariates
    
    em_param <- unlist(object$covariances$indv)
    
    sigma2 <- exp(em_param[1])
    delta_param <- em_param[-1]/2
    
    em_mat <- matrix(c(sigma2, delta_param), 
                     dimnames = list(list("sigma^2",em_names[-1]),"Variance"))
  }
  
  
  cat(" Type of model: ", model, "\n", "\n",
      
      "# Fixed effects: \n",
      "Formula: ", fe, "\n", "\n")
  
  print(fe_mat)
  
  cat("\n # Random effects: \n \n", 
      "Individual effects: \n",
      "Formula: ", ire, "\n", 
      "Covariance matrix \n")
  
  print(ire_mat)
  
  cat("\n Team effects: \n",
      "Formula: ", tre, "\n",
      "Covariance matrix \n")
  
  print(tre_mat)
  
  cat("\n # Emergence: \n",
      "Formula: ", em, "\n  \n")
  
  print(em_mat)
  
}

