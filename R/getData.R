



#' Title
#'
#' @param formula1 y ~ fixed
#' @param formula2 ~ indv.random | indv
#' @param formula3 ~ team.random | team
#' @param emergence ~ TIME + additional for methods CEM and CEI, ~ covariates (without time) for GP
#' @param method "CEM" Lang et al. model, "CEI" our altered model, "GP" Gaussian Process
#' @param method "OU - Guassian process for teams"
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
                    data,
                    method.team = NULL,
                    time = NULL) {
  
  ### Xf and y
  
  level1 <- formula1
  
  # extract info from formula
  terms1 <- terms(level1)
  
  vars1 <- as.character(attr(terms1,"variables")) # all variables in model in order they appear in formula
  vars1 <- vars1[-1] # remove "list" thing
  response <- vars1[1]
  
  # names
  f1_names <- attr(terms1,"term.labels")
  f1_int <- attr(terms1,"intercept")
  
  if (f1_int == 1) {
    f1_names <- c("(Intercept)", f1_names)
  } 
  
  # model frame
  mf1 <- model.frame(level1, data, na.action = na.pass) 
  
  y <- as.data.frame(model.extract(mf1, "response"))
  colnames(y) <- response
  
  Xf <- as.data.frame(model.matrix(level1, data=mf1))
  
  ## TODO: testa om full rank
  
  
  
  
  ### XI and indv
  
  l2 <- formula2
  l2_ch <- deparse1(l2) # converts it to string so the bar can be removed
  
  l2_list <- unlist(strsplit(l2_ch, "\\|")) # remove vertical bar, creates a list. unlist
  level2 <- as.formula(l2_list[1])
  
  # names
  terms2 <- terms(level2)
  f2_names <- attr(terms2,"term.labels")
  f2_int <- attr(terms2,"intercept")
  
  if (f2_int == 1) {
    f2_names <- c("(Intercept)", f2_names)
  }
 
  
  # model frame
  mf2 <- model.frame(level2, data, na.action = na.pass) 
  
  # XI
  XI <- as.data.frame(model.matrix(level2, data=mf2))
  
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
  
  # names
  terms3 <- terms(level3)
  f3_names <- attr(terms3,"term.labels")
  f3_int <- attr(terms3,"intercept")
  
  if (f3_int == 1) {
    f3_names <- c("(Intercept)", f3_names)
  } 
 
  
  # model frame
  mf3 <- model.frame(level3, data, na.action = na.pass) 
  
  # XT
  XT <- as.data.frame(model.matrix(level3, data=mf3))
  
  
  ## create team 
  groupvar3 <- l3_list[2]
  
  groupvar3 <- sub(" ", "", groupvar3)
  
  team <- as.data.frame(data[,groupvar3])
  
  colnames(team) <- groupvar3
  
  
  ### consensus method
  # possible values of method: "CEM", "CEI", "GP"
  ceFormula <- emergence
  
  
  # model/data frame
  ceMf <- model.frame(ceFormula, data, na.action = na.pass) 
  
  ceMf <- as.data.frame(model.matrix(ceFormula, data=ceMf))
  
  # names
  termsE <- terms(ceFormula)
  e_names <- attr(termsE,"term.labels")
  e_int <- attr(termsE,"intercept")
  
  if (e_int == 1) {
    e_names <- c("(Intercept)", e_names)
  }
  
  # create wNOISE/WEI/TI/(WET, TODO)
  
  ## initiate all extra wanted matrices
  WEI <- NULL
  WET <- NULL 
  TI <- NULL
  wNOISE <- NULL
  
  
  if (method == "CEM") {
    
    wNOISE <- ceMf
    
  } else if (method == "CEI") {
    
    WEI <- ceMf
    
  } else if (method == "GP") {
    
    TI <- ceMf
    
  } else {
    
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
    if(is.null(method.team)==F){
      if(method.team == "OU" || method.team == "OU.homeostasis"){
        timeName <- time
        time <- data[,timeName]
      }
    }else{
    time <- NULL 
    }
  }
  
  names <- list(f1_names, f2_names, f3_names, e_names)
  names(names) <- c("Fixed effects",
                    "Individual random effects", 
                    "Team random effects", 
                    "Emergence")
  
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
                   wNOISE,
                   names,
                   method,
                   method.team)
  
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
                       "wNOISE",
                       "names",
                       "method",
                       "method.team")
  
  return(matrices)
  
}
