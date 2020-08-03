

#'
#' creates optimization object, need to be orderd within team, indv
#'
#'
#'
#' team      - (n x 1) list of team y belongs to
#' indv      - (n x 1) list of indv y belongs to
#' Xf        - fixed effects
#' XI        - indivual effects (random)
#' XT        - team effects (random)
#' WLI       - Linear weighting indivual effect
#' WLT       - expontial weighting team effect
#' WEI       - expontial weighting indivual effect
#' WET       - expontial weighting team effect
#' WNOISE    - covariates to the measurement error
#' dataOrder - for filtering
dataToObject <- function(y,
                         team,
                         indv,
                         Xf = NULL,
                         XI = NULL,
                         XT = NULL,
                         WLI = TRUE,
                         WLT = TRUE,
                         wEI = NULL,
                         wET = NULL,
                         wNOISE = NULL,
                         dataOrder = NULL){

  TeamObj <- list()
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
  count <- 1
  if(is.null(XT)==F){
    if(WLT){
      d <-ncol(as.matrix(XT))
      TeamObj$teamCovs[[count]] <- Xcov$new(d)
      count <- count + 1
    }
  }
  if(is.null(wET)==F){
    d <- ncol(as.matrix(XT))
    TeamObj$teamCovs[[count]] <- XcovSmooth$new(d)
  }
  count <- 1
  if(is.null(XI)==F){
    if(WLI){
      d <-ncol(as.matrix(XI))
      TeamObj$indvCovs[[count]] <- Xcov$new(d)
      count <- count + 1
    }
  }
  if(is.null(wEI)==F){
    d <- ncol(as.matrix(XI))
    TeamObj$indvCovs[[count]] <- XcovSmooth$new(d)
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
    if(is.null(wET)==F)
      TeamObj$teams[[i]]$w = as.matrix(wET)[Teams == uTeams[i],, drop = FALSE]

    TeamObj$teams[[i]]$indv <- list()
    indv_i <-  levels(TeamObj$teams[[i]]$data$indv)
    for(ii in 1:TeamObj$teams[[i]]$nindv){
      TeamObj$teams[[i]]$indv[[ii]] <- list()
      index <- TeamObj$teams[[i]]$data$indv == indv_i[ii]

      TeamObj$teams[[i]]$indv[[ii]]$n = sum(index)

      TeamObj$teams[[i]]$indv[[ii]]$A      <-  Matrix::sparseMatrix(j = 1:sum(index),
                                                 i = which(index),
                                                 x = rep(1, sum(index)),
                                                 dims = c(n_obs,
                                                          sum(index)))
      TeamObj$teams[[i]]$indv[[ii]]$data <- list()
      if(is.null(XI)==F){
        XIt <-  as.matrix(XI)[Teams == uTeams[i],, drop = FALSE]
        TeamObj$teams[[i]]$indv[[ii]]$X = as.matrix(XIt)[index,, drop = FALSE]
      }
      if(is.null(wEI)==F){
        wIt <- as.matrix(wEI)[Teams == uTeams[i],, drop = FALSE]
        TeamObj$teams[[i]]$indv[[ii]]$w =as.matrix(wIt)[index,, drop = FALSE]
      }
    }
  }
  return(TeamObj)
}
