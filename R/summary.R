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
  
  tre_param <- unlist(object$covariances$team[1])
  tre_mat <- as.matrix(nlme::pdLogChol(tre_param))
  dimnames(tre_mat) <- list(tre_names,tre_names)
  
  # team GP
  if (!is.null(object$model$`Team process`)) {
    
    tgp_param <- unlist(object$covariances$team[[2]])
    
    sigma <- exp(tgp_param[1])
    theta <- exp(tgp_param[2])
    
    
    tgp_mat <- matrix(c(sigma,theta),
                      dimnames = list(c("sigma","theta"),"Estimate"))
  
  }
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
    # individual
    
    
    # emergence (on individual level for CEI)
    
    em <- as.character(object$model$`Emergence model`[1]) # formula
    em_names <- unlist(object$model$`Emergence model`[2]) # names of covariates
    
    em_param <- unlist(object$covariances$indv)[1:length(em_names)]
    
    if (em_names[1]=="(Intercept)"){    
      indv_baseline <- exp(em_param[1]) # inte rätt
      
      delta_param <- em_param[-1]/2 # inte rätt
      
      
      em_mat <- matrix(c(indv_baseline, delta_param), 
                       dimnames = list(c("individual baseline variance", em_names[-1]),"Estimate"))
    } else {
      delta_param <- em_param/2 # ska inte delas med två
      
      
      em_mat <- matrix(delta_param, 
                       dimnames = list(em_names,"Estimate"))
    }
    
    # individual
    
    ire <- as.character(object$model$`Individual random effects`[1]) # formula
    ire_names <- unlist(object$model$`Individual random effects`[2]) # names of covariates
    
    n_indv_param <- length(unlist(object$covariances$indv))
    
    ire_param <- unlist(object$covariances$indv)[(length(em_names)+1):n_indv_param]
    ire_mat <- as.matrix(nlme::pdLogChol(ire_param))
    dimnames(ire_mat) <- list(ire_names,ire_names)
    
    # sigma_u 
    
    
    # team
      # random effects
    
    tre <- as.character(object$model$`Team random effects`[1]) # formula
    tre_names <- unlist(object$model$`Team random effects`[2]) # names of covariates
    
    tre_param <- unlist(object$covariances$team[1])
    tre_mat <- as.matrix(nlme::pdLogChol(tre_param))
    
    dimnames(tre_mat) <- list(tre_names,tre_names)
    
      # team GP
    if (!is.null(object$model$`Team process`)) {
      
    tgp_param <- unlist(object$covariances$team[[2]])
    
    sigma <- exp(tgp_param[1])
    theta <- exp(tgp_param[2])
    
 
    tgp_mat <- matrix(c(sigma,theta),
                       dimnames = list(c("sigma","theta"),"Estimate"))

    }
    
    
    
  } else if (model == "GP") {
    
    ## fixed effects
    fe <- as.character(object$model$`Fixed effects`[1]) # formula
    fe_names <- unlist(object$model$`Fixed effects`[2]) # names of covariates
    fe_param <- object$betas
    
    fe_mat <- cbind(fe_param,diag(t(object$cov_beta)))
    dimnames(fe_mat) <- list(fe_names,c("Estimate", "Variance"))
    
    fe_measurement_error <- exp(unlist(object$covariances$error)) # do we want to print this?
    
    ## random effects
    
    # individual
    
    ire <- as.character(object$model$`Individual random effects`[1]) # formula
    ire_names <- unlist(object$model$`Individual random effects`[2]) # names of covariates
    
    ire_param <- unlist(object$covariances$indv[[1]])
    ire_mat <- as.matrix(nlme::pdLogChol(ire_param))
    dimnames(ire_mat) <- list(ire_names,ire_names)
    
    # sigma_u 
    
    
    # team
    
    tre <- as.character(object$model$`Team random effects`[1]) # formula
    tre_names <- unlist(object$model$`Team random effects`[2]) # names of covariates
    
    tre_param <- unlist(object$covariances$team[1])
    tre_mat <- as.matrix(nlme::pdLogChol(tre_param))
    dimnames(tre_mat) <- list(tre_names,tre_names)
    
    # team GP
    if (!is.null(object$model$`Team process`)) {
      
      tgp_param <- unlist(object$covariances$team[[2]])
      
      sigma <- exp(tgp_param[1])
      theta <- exp(tgp_param[2])
      
      
      tgp_mat <- matrix(c(sigma,theta),
                        dimnames = list(c("sigma","theta"),"Estimate"))
    }
    # emergence GP
    
    
    em <- as.character(object$model$`Emergence model`[1]) # formula
    em_names <- unlist(object$model$`Emergence model`[2])[-1] # names of covariates
    
  
    em_param <- unlist(object$covariances$indv[[2]])
    n_param <- length(em_param)
    
    beta_delta <- em_param[1]
    sigma <- exp(em_param[n_param-1])
    theta <- exp(em_param[n_param])
    
    if (n_param > 3) {
      
      additional_param <- em_param[2:(n_param-2)]
      
      em_mat <- matrix(c(beta_delta,sigma,theta, additional_param),
                       dimnames = list(c("beta_delta","sigma","theta",em_names),"Estimate"))
    } else {
      em_mat <- matrix(c(beta_delta,sigma,theta),
                       dimnames = list(c("beta_delta","sigma","theta"),"Estimate"))
    }
  }
    
  
  # model fit
  loglik <- object$loglik
  
  
  k <- sum(length(unlist(object$covariances)),
           length(fe_param))
  
  n <- object$n
  
  ## 
  # AIC: 2*k - 2*loglik, k = number of estimated parameters  
  aic <- (2*k)-(2*loglik)
  
  # BIC: k*ln(n)-2*loglik, n=nr of observations
  bic <- (k*log(n))-(2*loglik)
  
  # fit measures
  fit <- c(loglik, aic, bic, round(k))
  names(fit) <- list("logLik", "AIC", "BIC", "# of parameters")
  
  cat(" Type of model: ", model, "\n", "\n")
      
      print(fit) 
      
  cat("\n", "\n","# Fixed effects: \n",
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
  
  if (!is.null(object$model$`Team process`)) {
    cat("\n Team GP: \n")
    
    print(tgp_mat)
  }
  
  cat("\n # Emergence: \n",
      "Formula: ", em, "\n  \n")
  
  print(em_mat)
  
}

