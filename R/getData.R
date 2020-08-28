



#' Title
#'
#' @param formula1 y ~ fixed
#' @param formula2 ~ indv.random | indv
#' @param formula3 ~ team.random | team
#' @param emergence ~ TIME + additional for methods CEM and CEI, ~ covariates (without time) for GP
#' @param method "CEM" Lang et al. model, "CEI" our altered model, "GP" Gaussian Process
#' @param time name of time variable in quote, "TIME"
#' @param data 
#'
#' @return 
#' team      - (n x 1) list of team y belongs to
#' indv      - (n x 1) list of indv y belongs to
#' Xf        - fixed effects
#' XI        - indivual effects (random)
#' XT        - team effects (random)
#' WEI       - expontial weighting indivual effect
#' WNOISE    - covariates to the measurement error

#' @export
#'
#' @examples
#' 
getData <- function(formula1, 
                    formula2, 
                    formula3, 
                    emergence, 
                    method,
                    time = NULL,
                    data) {
  
  ### Xf and y
  
  level1 <- formula1
  
  # extract info from formula
  terms1 <- terms(level1)
  term1_names <- attr(terms1,"term.labels")
  
  
  vars1 <- as.character(attr(terms1,"variables")) # all variables in model in order they appear in formula
  vars1 <- vars1[-1] # remove "list" thing
  response <- vars1[1]
  
  # model frame
  #mf1 <- model.frame(level1, data, na.action = na.pass) 
  mf1 <- as.data.frame(data[,vars1])
  
  # create Xf and y

  if (ncol(mf1) == 1) { 
     
     y <- as.data.frame(mf1)
     colnames(y) <- response
     
     # add intercept to Xf
     
     n <- nrow(mf1)
     
     if (attr(terms1,"intercept") == 1) {
       intercept <- rep(1, n)
       Xf <- as.data.frame(intercept)
     }
     
  } else {
    Xf <- as.data.frame(mf1[,term1_names])
    colnames(Xf) <- term1_names
    
    y <- as.data.frame(mf1[,response]) 
    colnames(y) <- response
    
    # add intercept to Xf
    
    n <- nrow(mf1)
    
    if (attr(terms1,"intercept") == 1) {
      intercept <- rep(1, n)
      Xf <- data.frame(intercept,Xf)
    }
  }
  

  
  ## TODO: testa om full rank
  
   
  
  
  ### XI and indv
  
  
  
  l2 <- formula2
  l2_ch <- deparse1(l2) # converts it to string so the bar can be removed
  
  l2_list <- unlist(strsplit(l2_ch, "\\|")) # remove vertical bar, creates a list. unlist
  level2 <- as.formula(l2_list[1])
  
  
  # extract info from formula
  terms2 <- terms(level2)
  term2_names <- attr(terms2,"term.labels")
  
  
  vars2 <- as.character(attr(terms2,"variables")) # all variables in model in order they appear in formula
  vars2 <- vars2[-1] # remove "list" thing
  
  
  # model frame
  #mf2 <- model.frame(level2, data, na.action = na.pass) 
   mf2 <- as.data.frame(data[,vars2])
  
  # create XI
  if (ncol(mf2) == 1) {
    XI <- mf2
  } else {
    XI <- mf2[,term2_names]
  }
  
  colnames(XI) <- term2_names 
  
  # add intercept to XI
  
  if (attr(terms2,"intercept") == 1) {
    intercept <- rep(1, n) # will this be a problem if all intercepts are called "intercept"?
    XI <- data.frame(intercept,XI)
  }
  
  
  ## create indv (does it need to be unique within teams or overall?)
  groupvar2 <- l2_list[2]
  
  groupvar2 <- sub(" ", "", groupvar2)
  
  indv <- as.data.frame(data[,groupvar2])
  
  colnames(indv) <- groupvar2
  
  
  ### XT and team
  
  l3 <- formula3
  l3_ch <- deparse1(l3) # converts it to string so the bar can be removed
  
  l3_list <- unlist(strsplit(l3_ch, "\\|")) # remove vertical bar, creates a list. unlist
  level3 <- as.formula(l3_list[1])
  
  
  # extract info from formula
  terms3 <- terms(level3)
  term3_names <- attr(terms3,"term.labels")
  
  
  vars3 <- as.character(attr(terms3,"variables")) # all variables in model in order they appear in formula
  vars3 <- vars3[-1] # remove "list" thing
  
  
  # model frame
  #mf3 <- model.frame(level3, data, na.action = na.pass) 
  mf3 <- as.data.frame(data[,vars3])
  
  # create XT
  if (ncol(mf3) == 1) {
    XT <- mf3
  } else {
    XT <- mf3[,term3_names]
  }
  
  colnames(XT) <- term3_names 
  
  # add intercept to XT
  
  if (attr(terms3,"intercept") == 1) {
    intercept <- rep(1, n) 
    XT <- data.frame(intercept,XT)
  }
  
  
  ## create team 
  groupvar3 <- l3_list[2]
  
  groupvar3 <- sub(" ", "", groupvar3)
  
  team <- as.data.frame(data[,groupvar3])
  
  colnames(team) <- groupvar3
  
  
  ### consensus method
  # possible values of method: "CEM", "CEI", "GP"
  ceFormula <- emergence
  
  # extract info from formula
  ceTerms <- terms(ceFormula)
  ceTerms_names <- attr(ceTerms,"term.labels")
  
  
  ceVars <- as.character(attr(ceTerms,"variables")) # all variables in model in order they appear in formula
  ceVars <- ceVars[-1] 
  
  # model frame
  #ceMf <- model.frame(ceFormula, data, na.action = na.pass)
  ceMf <- as.data.frame(data[,ceVars])
  
  
  # create wNOISE/WEI/TI/(WET, TODO)
  
  ## initiate all extra wanted matrices
  WEI <- NULL
  WET <- NULL 
  TI <- NULL
  wNOISE <- NULL
  
  
  if (ncol(ceMf) == 1) {
    
    if (method == "CEM") {
      
      wNOISE <- ceMf
      colnames(wNOISE) <- ceTerms_names
      
      
    } else if (method == "CEI") {
      
      WEI <- ceMf
      colnames(WEI) <- ceTerms_names
      
    } else if (method == "GP") {
      
      TI <- ceMf
      colnames(TI) <- ceTerms_names
      
    } else {
      
    }
  } else {
    
    if (method == "CEM") {
      
      wNOISE <- ceMf[,ceTerms_names]
      colnames(wNOISE) <- ceTerms_names
 
      
    } else if (method == "CEI") {
      
      WEI <- ceMf[,ceTerms_names]
      colnames(WEI) <- ceTerms_names
      
      
      
    } else if (method == "GP"){
      TI <- ceMf[,ceTerms_names]
      colnames(TI) <- ceTerms_names
      
      
    } else {
 
    }
  }
  
  # create time
  if (method == "GP" && is.null(time) == TRUE) {
    
    warning("Time variable not specified")
    
  } else if (method == "GP" && exists(time,TI) == TRUE) {
    
    # TODO: find a way to do this without specifying the name as
    #       a character in argument
    
    warning("Emergence variables cannot include time when method = GP")
  
  }
    else if (method == "GP") {
    timeName <- time
    time <- data[,timeName]
    
    # TODO: add name to time variable? 
    
    
  } else {
    time <- NULL 
  }
  
  
  matrices <- list(y, 
                   Xf, 
                   XI, 
                   indv, 
                   XT, 
                   team, 
                   WEI, 
                   WET, 
                   TI, 
                   time, 
                   wNOISE)
  
  names(matrices) <- c("y", 
                       "Xf", 
                       "XI", 
                       "indv", 
                       "XT", 
                       "team", 
                       "WEI", 
                       "WET", 
                       "TI", 
                       "time", 
                       "wNOISE")
  
  return(matrices)
  
}
