#'
#' creates optimization object, need to be ordered within team, indv
#'
#'
#'
#' team      - (n x 1) vector of team y belongs to
#' indv      - (n x 1) vector of indv y belongs to
#' Xf        - fixed effects
#' XI        - individual effects (random)
#' XT        - team effects (random)
#' WLI       - linear weighting individual effect 
#' WLT       - linear weighting team effect
#' WEI       - exponential weighting individual effect
#' WII       - independent weighting of inidvual effect 
#' WET       - exponential weighting team effect
#' TI        -  n_indv x p (time individual effect (covariance model) requires time)
#' time      - (n x 1) time point of observations
#' WNOISE    - covariates to the measurement error
#' dataOrder - for filtering
#'
#'
#' @param data_list
#'        $method.team - method for time series team (allow null)
#'@export
dataToObject <- function(data_list, Teams = NULL){
  
  # data_list - output list from function getData()
  data <- data_list
  
  y <- as.matrix(data[["y"]])
  team <- as.matrix(data[["team"]])
  indv <- as.matrix(data[["indv"]])
  Xf <- as.matrix(data[["Xf"]])
  XI <- as.matrix(data[["XI"]])
  XT <- as.matrix(data[["XT"]])
  
  method <- as.character(data[["method"]])
  
  if (method == "CEM") {
    WLI <- T
  } else if (method == "CEM2"){
    WLI <- T
  } else if (method == "CEI2"){
    WLI <- T
  } else if (method == "GP"){
    WLI <- T
  } else {
    WLI <- F
  }
  
  WLT <- T  
  
  if (is.null(data[["WEI"]]) == F) {
    WEI <- as.matrix(data[["WEI"]])
  } else {
    WEI <- data[["WEI"]]
  }
  
  if (is.null(data[["WII"]]) == F) {
    WII <- as.matrix(data[["WII"]])
  } else {
    WII <- data[["WII"]]
  }
  if (is.null(data[["WET"]]) == F) {
    WET <- as.matrix(data[["WET"]])
  } else {
    WET <- data[["WET"]]
  }
  
  if (is.null(data[["TI"]]) == F) {
  TI <- as.matrix(data[["TI"]])
  } else {
    TI <- data[["TI"]]
  }
  
  if (is.null(data[["time"]]) == F) {
  time <- as.matrix(data[["time"]])
  } else {
    time <- data[["time"]]
  }
  
  if (is.null(data[["wNOISE"]]) == F) {
    wNOISE <- as.matrix(data[["wNOISE"]])
  } else {
    wNOISE <- (data[["wNOISE"]])
  }
  
  dataOrder <- NULL

  TeamObj <- list()
  if(is.null(Teams))
    Teams <- factor(team)
  uTeams <- unique(Teams)
  TeamObj$indexTeams <- Teams
  if(is.null(Xf)==F){
    TeamObj$X          <- as.matrix(Xf)
  }else{
    TeamObj$X <- Xf
  }
  TeamObj$nTeams     <- length(uTeams)
  TeamObj$teams      <- list()
  TeamObj$teamCovs   <- list()
  TeamObj$indvCovs   <- list()
  TeamObj$errorCovs  <- list()
  if(is.null(wNOISE)){
    wNOISE = as.matrix(rep(1,length(y)))
  }
  TeamObj$errorCovs[[1]] <- expWeightDiag$new(dim(wNOISE)[2])
  
  
  ######
  # setting up team random effects
  ######
  count <- 1
  
  if(is.null(XT)==F){
    if(WLT){
      d <-ncol(as.matrix(XT))
      TeamObj$teamCovs[[count]] <- Xcov$new(d)
      count <- count + 1
    }
  
    if(is.null(WET)==F){
      d <- ncol(as.matrix(XT))
      TeamObj$teamCovs[[count]] <- XcovSmooth$new(d, dim(WET)[2])
    }
  }
  if(is.null(data[["method.team"]])==F){
    method.team <- as.character(data[["method.team"]])
    if(method.team == "OU"){
      TeamObj$teamCovs[[count]] <- OU$new()
      count <- count + 1
    }else if(method.team == "OU.homeostasis"){
      TeamObj$teamCovs[[count]] <- OU.homeostasis$new(1) # two parameters time and delta
      count <- count + 1
      
    }
  }
  
  
  ###  
  # setting up the indiuval random effects
  ###
  count <- 1
  if(is.null(XI)==F & length(XI)>0){
    if(WLI){
      d <-ncol(as.matrix(XI))
      TeamObj$indvCovs[[count]] <- Xcov$new(d)
      count <- count + 1
    }
  
    if(is.null(WEI)==F){
      d <- ncol(as.matrix(XI))
      TeamObj$indvCovs[[count]] <- XcovSmooth$new(d, dim(WEI)[2])
      count <- count + 1
    }
  }
  if(is.null(WII)==F){
    TeamObj$indvCovs[[count]] <- expWeightDiag$new(dim(WII)[2])
    #OUbridge$new(min(time), max(time), dim(TI)[2])
    count <- count + 1
  }
  if(is.null(TI)==F){
    if(data[["GP.type"]] == "OU.homeostasis")
      TeamObj$indvCovs[[count]] <- OU.homeostasis$new(1) # two parameters time and delta 
    if(data[["GP.type"]] == "BM.homeostasis")
      TeamObj$indvCovs[[count]] <- BM.homeostasis$new(1) # two parameters time and delta 
    #OUbridge$new(min(time), max(time), dim(TI)[2])
    count <- count + 1
  }

  for(i in 1:TeamObj$nTeams){
    TeamObj$teams[[i]] <-list()
    TeamObj$teams[[i]]$name <- uTeams[i]
    n_obs <- sum(Teams == uTeams[i])
    TeamObj$teams[[i]]$n_obs <-  n_obs
    TeamObj$teams[[i]]$data <- list(Y     = y[Teams == uTeams[i]],
                                    indv  = factor(indv[Teams == uTeams[i]]))
    TeamObj$teams[[i]]$E <- wNOISE[Teams == uTeams[i],, drop=FALSE]
    if(is.null(Xf)==F)
      TeamObj$teams[[i]]$data$X <- TeamObj$X[Teams== uTeams[i],, drop=FALSE]
    TeamObj$teams[[i]]$nindv   <- length(levels(TeamObj$teams[[i]]$data$indv))

    TeamObj$teams[[i]]$A      <-  Matrix::sparseMatrix(j = 1:n_obs,
                                               i = which(Teams == uTeams[i]),
                                               x = rep(1, n_obs),
                                               dims = c(length(y),
                                                        n_obs))
    if(is.null(XT)==F)
      TeamObj$teams[[i]]$X = as.matrix(XT)[Teams == uTeams[i],, drop = FALSE]
    if(is.null(WET)==F)
      TeamObj$teams[[i]]$W = as.matrix(WET)[Teams == uTeams[i],, drop = FALSE]

    if(is.null(data[["method.team"]])==F){
      if(method.team == "OU" || method.team == "OU.homeostasis"){
        timeTeam                <- time[Teams == uTeams[i]]
        TeamObj$teams[[i]]$time  <- timeTeam
        TeamObj$teams[[i]]$D    <- as.matrix(timeTeam)
        
      }
    }
    
    TeamObj$teams[[i]]$indv <- list()
    indv_i <-  levels(TeamObj$teams[[i]]$data$indv)
    TeamObj$teams[[i]]$indvNames <- indv_i  
  
    
    ##
    # setting up indivuals data
    ##
    for(ii in 1:TeamObj$teams[[i]]$nindv){
      TeamObj$teams[[i]]$indv[[ii]] <- list()
      index <- TeamObj$teams[[i]]$data$indv == indv_i[ii]

      TeamObj$teams[[i]]$indv[[ii]]$n = sum(index)
      TeamObj$teams[[i]]$indv[[ii]]$A_list = which(index) 
      TeamObj$teams[[i]]$indv[[ii]]$A      <-  Matrix::sparseMatrix(j = 1:sum(index),
                                                 i = which(index),
                                                 x = rep(1, sum(index)),
                                                 dims = c(n_obs,
                                                          sum(index)))
      TeamObj$teams[[i]]$indv[[ii]]$data <- list()
      #is there a random linear effect
      if(is.null(XI)==F){
        XIt <-  as.matrix(XI)[Teams == uTeams[i],, drop = FALSE]
        TeamObj$teams[[i]]$indv[[ii]]$X = as.matrix(XIt)[index,, drop = FALSE]
      }
      #is there a expontialy weighted effect
      if(is.null(WEI)==F){
        wIt <- as.matrix(WEI)[Teams == uTeams[i],, drop = FALSE]
        TeamObj$teams[[i]]$indv[[ii]]$W =as.matrix(wIt)[index,, drop = FALSE]
      }
      if(is.null(WII)==F){
        WIIt <-  as.matrix(WII)[Teams == uTeams[i],, drop = FALSE]
        TeamObj$teams[[i]]$indv[[ii]]$E = as.matrix(WIIt)[index,, drop = FALSE]
        
      }
      #is there is a time series component
      if(is.null(TI)==F){
        timeTeam                <- time[Teams == uTeams[i]]
        COVTI      <- TI[Teams == uTeams[i], , drop=FALSE]
        COVTI_indv <- COVTI[index, , drop=FALSE]
        TeamObj$teams[[i]]$indv[[ii]]$time <- timeTeam[index]
        #TeamObj$teams[[i]]$indv[[ii]]$D    <- COVTI_indv[1,, drop=FALSE]
        TeamObj$teams[[i]]$indv[[ii]]$D     <- as.matrix(timeTeam[index])
      }
    }
  }
  return(TeamObj)
}
