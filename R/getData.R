



#' Title
#'
#' @param formula1 y ~ fixed
#' @param formula2 ~ indv.random | indv
#' @param formula3 ~ team.random | team
#' @param emergence ~ TIME + additional 
#' @param model "CEM" Lang et al. model, "CEI" our altered model, "GP"
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
                    model, 
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
  #mf <- model.frame(level1, data) # removes NA's
  mf1 <- data[,vars1]
  
  # create Xf and y
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
  #mf <- model.frame(level2, univbct) # removes NA's
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
  #mf <- model.frame(level3, univbct) # removes NA's
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
    intercept <- rep(1, n) # will this be a problem if all intercepts are called "intercept"?
    XT <- data.frame(intercept,XT)
  }
  
  
  ## create team 
  groupvar3 <- l3_list[2]
  
  groupvar3 <- sub(" ", "", groupvar3)
  
  team <- as.data.frame(data[,groupvar3])
  
  colnames(team) <- groupvar3
  
  
  ### consensus model
  # possible values of model: "CEM", "CEI", "GP"
  ceFormula <- emergence
  
  # extract info from formula
  ceTerms <- terms(ceFormula)
  ceTerms_names <- attr(ceTerms,"term.labels")
  
  
  ceVars <- as.character(attr(ceTerms,"variables")) # all variables in model in order they appear in formula
  ceVars <- ceVars[-1] 
  
  # model frame
  ceMf <- as.data.frame(data[,ceVars])
  
  
  # create wNOISE/WEI/(GP, TODO)
  
  
  if (ncol(ceMf) == 1) {
    
    if (model == "CEM") {
      
      wNOISE <- ceMf
      colnames(wNOISE) <- ceTerms_names
      
      matrices <- list(y, Xf, XI, indv, XT, team, wNOISE)
      names(matrices) <- c("y", "Xf", "XI", "indv", "XT", "team", "wNOISE")
      
    } else if (model == "CEI") {
      
      WEI <- ceMf
      colnames(WEI) <- ceTerms_names
      
      matrices <- list(y, Xf, XI, indv, XT, team, WEI)
      names(matrices) <- c("y", "Xf", "XI", "indv", "XT", "team", "WEI")
      
   # } else if (model == "GP") {
      # TODO
    } else {
      matrices <- list(y, Xf, XI, indv, XT, team)
      names(matrices) <- c("y", "Xf", "XI", "indv", "XT", "team")
    }
  } else {
    
    if (model == "CEM") {
      
      wNOISE <- ceMf[,ceTerms_names]
      colnames(wNOISE) <- ceTerms_names
      
      matrices <- list(y, Xf, XI, indv, XT, team, wNOISE)
      names(matrices) <- c("y", "Xf", "XI", "indv", "XT", "team", "wNOISE")
      
    } else if (model == "CEI") {
      
      WEI <- ceMf[,ceTerms_names]
      colnames(WEI) <- ceTerms_names
      
      matrices <- list(y, Xf, XI, indv, XT, team, WEI)
      names(matrices) <- c("y", "Xf", "XI", "indv", "XT", "team", "WEI")
      
   # } else if (model == "GP"){
      # TODO
    } else {
      matrices <- list(y, Xf, XI, indv, XT, team)
      names(matrices) <- c("y", "Xf", "XI", "indv", "XT", "team")
    }
  }
  
  
  return(matrices)
  
}
