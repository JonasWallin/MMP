library(Matrix)
library(MASS)



#'
#' takes a dataset and creates an object suitable for estimation
#'
#'
#' @param data - data.frame
#' @param time - time variable        name in data
#' @param y    - observation variable name in data
#' @param team - team variable        name in data
#'
dataToObject_old <- function(data,
                             time  ='time',
                             y     ='y',
                             team  ='team',
                              indv = 'indv'){

  TeamObj <- list()
  Teams <- factor(data[, team])
  uTeams <- unique(Teams)
  TeamObj$nTeams  <- length(uTeams)
  TeamObj$teams <- list()
  for(i in 1:TeamObj$nTeams){
    TeamObj$teams[[i]] <-list()
    TeamObj$teams[[i]]$name <- uTeams[i]

    TeamObj$teams[[i]]$data <- data.frame(Y     = data[Teams == uTeams[i],y],
                                          time  = data[Teams == uTeams[i],time],
                                          indv  = factor(data[Teams == uTeams[i],indv]))
    TeamObj$teams[[i]]$utime   <- unique(TeamObj$teams[[i]]$data$time)
    TeamObj$teams[[i]]$dMAtrix <- as.matrix(dist(TeamObj$teams[[i]]$utime))
    TeamObj$teams[[i]]$nindv   <- length(levels(TeamObj$teams[[i]]$data$indv))
    Aj = match(TeamObj$teams[[i]]$data$time,
               TeamObj$teams[[i]]$utime)
    n_ut <- length(TeamObj$teams[[i]]$utime)
    TeamObj$teams[[i]]$AO      <-  sparseMatrix(i = 1:length(Aj),
                                                j = Aj,
                                                x = rep(1, length(Aj)),
                                                dims = c(length(Aj),
                                                         n_ut))
    indv_i <-  levels(TeamObj$teams[[i]]$data$indv)
    TeamObj$teams[[i]]$n_indv <- length(indv_i)
    TeamObj$teams[[i]]$AI <- list()
    for(ii in 1:TeamObj$teams[[i]]$n_indv){
      index <- TeamObj$teams[[i]]$data$indv == indv_i[ii]
      Aii <- match(TeamObj$teams[[i]]$data$time[index],
                   TeamObj$teams[[i]]$utime)
      TeamObj$teams[[i]]$AI[[ii]] <- sparseMatrix(i = 1:length(Aii),
                                                  j = Aii,
                                                  x = rep(1, length(Aii)),
                                                  dims = c(length(Aii),
                                                           n_ut))
    }



  }
  return(TeamObj)
}

#'
#' setting up process type of processes for team and
#' indiviual
#' @param TeamObj - list
#' @param indvP   - which process type should be used for indv
#'                  OU - ornstein uhlenbeck type process
#'                  R1 - rank one 'process'
#' @param teamP   - which process type should be used for team
SetObject <- function(TeamObj,
                      indvP = 'OU',
                      teamP = 'OU'){

  TeamObj$IndvProcess <- indvP
  TeamObj$TeamProcess <- teamP
  if(TeamObj$IndvProcess == 'OU'){
    TeamObj$Param_indv <- c(0,0,0)
  }else if(TeamObj$IndvProcess == 'R1'){
    TeamObj$Param_indv <- c(0,0)
    TeamObj$Indv_CovObj <- RankOneCov
  }else{
    TeamObj$Param_indv <- c(0)
    TeamObj$IndvProcess = 'constant'
  }

  if(TeamObj$TeamProcess == 'OU'){
    TeamObj$Param_team <- c(0,0,0,0)
  }else if(TeamObj$TeamProcess == 'R1'){
    TeamObj$Param_team <- c(0,0,0)
    TeamObj$Team_CovObj <- RankOneCov
  }else{
    TeamObj$Param_team <- c(0,0)
    TeamObj$TeamProcess =  'constant'
  }
  TeamObj$Param_noise <- c(0)
  return(TeamObj)
}

#'
#' Takes a team Object and simulate data
#'
#'
#' TODO: change the structure here so TeamProcess only
#'       has a covariance object to call
#'       Same for Idnv
simulate <- function(TeamObj){

  Data <- list(data=c(), latentVariables=list())
  Data$latentVariables$O <- c()
  Data$latentVariables$X <- c()
  Data$latentVariables$E <- c()
  data <- c()
  for(i in 1:length(TeamObj$teams)){
    Team_i <- TeamObj$teams[[i]]
    ##################################
    # simulating team processes
    ##################################
    O <- TeamObj$Param_team[1] + exp(TeamObj$Param_team[2]) * rep(rnorm(n=1), length(Team_i$utime))
    if(TeamObj$TeamProcess == 'OU'){
      sigmaO <- OUcov(Team_i$dMAtrix, TeamObj$Param_team[3:4])
      O <- mvrnorm(1, O, sigmaO)
    }else if(TeamObj$TeamProcess == 'R1'){
      sigmaO <- TeamObj$Team_CovObj$get_Cov(loc= Team_i$utime,
                                            param = TeamObj$Param_team[3:4])
      O <- mvrnorm(1, O, sigmaO)
    }
    Data$latentVariables$O <- rbind(Data$latentVariables$O,
                                    cbind(i,Team_i$utime,O ))
    #################################
    #  Simulating indivual processes
    #################################
    if(TeamObj$IndvProcess == 'OU'){
      sigmaI <- OUcov(Team_i$dMAtrix, TeamObj$Param_indv[2:3])
    }else if(TeamObj$IndvProcess == 'R1'){

      sigmaI <- TeamObj$Indv_CovObj$get_Cov(loc   = Team_i$utime,
                                            param = TeamObj$Param_indv[1:2])
    }
    Xi     <- c()
    indv   <- c()
    time_i <- c()
    for(ii in 1:Team_i$n_indv){
      Aindv <- Team_i$AI[[ii]]
      Xii <- exp(TeamObj$Param_indv[1]) * rep( rnorm(1),dim(Aindv)[1])
      if(TeamObj$IndvProcess == 'OU'){
        sigmaIii <- Aindv%*% sigmaI%*% Matrix::t(Aindv)
        Xii <- mvrnorm(1, Xii, sigmaIii)
      }else if(TeamObj$IndvProcess == 'R1'){
        sigmaIii <- Aindv%*% sigmaI%*% Matrix::t(Aindv)
        Xii <- mvrnorm(1, mu=rep(0,dim(Aindv)[1]), sigmaIii)
      }
      time_ii <- as.vector(Aindv%*%Team_i$utime)
      n_ii    <- length(time_ii)
      indv   <- c(indv, rep(ii, n_ii))
      Xi     <- c(Xi, Xii)
      time_i <- c(time_i, time_ii)

      Data$latentVariables$X <- rbind(Data$latentVariables$X,
                                      cbind(i, ii, time_ii, Xii ))

    }
    #################################
    #  Simulating noise processes
    #################################
    E_i <- exp(TeamObj$Param_noise) * rnorm( n = length(time_i))
    Data$latentVariables$E <- c(Data$latentVariables$E, E_i)

    Y_i <- Xi + E_i + O
    data <- rbind(data,
                  cbind(i, indv, time_i, Y_i))
  }
  Data$data <- data.frame(time = data[,3],
                          team = data[,1],
                          indv = data[,2],
                          y    = data[,4])
  return(Data)
}

#'
#' estimate parameters
#'
#'
#' @param TeamObj
#'
#'
estimate <- function(TeamObj){
  param <- c(TeamObj$Param_noise)
  param <- c(param, TeamObj$Param_team)
  param <- c(param, TeamObj$Param_indv)
  neglog <- function(param){ -loglik(param, TeamObj)}
  res <- optim(param, neglog)
  n_param <- 1
  TeamObj$Param_noise <- res$par[1]
  TeamObj$Param_team <- res$par[n_param + 1:length(TeamObj$Param_team)]
  n_param <- n_param + length(TeamObj$Param_team)
  TeamObj$Param_indv <- res$par[n_param + 1:length(TeamObj$Param_indv)]
  TeamObj$param <- res$param
  return(TeamObj)
}

#'
#' Creates a vector of parameters into a list where each parameter
#' param     - (p x 1) parameter vector
#' paramList - (List)
#'             $error - (1 x 1) log(sigma) Gaussian measurement error
#'             $indv  -  (list) one for each indiviuals covariance function
#'             $team  -  (list) one for each team       covariance function
#' D: 2020-03-16
#'
#'
paramToList <- function(param, obj){

  paramList <- list()
  paramList$team  <- NULL
  paramList$indv  <- NULL
  paramList$error <- NULL


  paramList$error <- list()
  for(i in 1:length(obj$errorCovs)){
    covObj <- obj$errorCovs[[i]]
    nParObj <- covObj$get_param_length()
    paramList$error[[i]] <- param[1:nParObj]
    param <- param[-(1:nParObj)]
  }

  if(length(obj$indvCovs)>0){
    paramList$indv <- list()
    for(i in 1:length(obj$indvCovs)){
      covObj <- obj$indvCovs[[i]]
      nParObj <- covObj$get_param_length()
      paramList$indv[[i]] <- param[1:nParObj]
      param <- param[-(1:nParObj)]
    }
  }

  if(length(obj$teamCovs)>0){
    paramList$team <- list()
    for(i in 1:length(obj$teamCovs)){
      covObj <- obj$teamCovs[[i]]
      nParObj <- covObj$get_param_length()
      paramList$team[[i]] <- param[1:nParObj]
      param <- param[-(1:nParObj)]
    }
  }
  return(paramList)
}
#'
#'
#' creates a zero guess of the parameters list 
#' 
#'
#'
param0 <- function(obj){

  param <- c()
  for(i in 1:length(obj$errorCov))
    param <- c(param, rep(0, obj$errorCovs[[i]]$get_param_length())) 

  if(length(obj$indvCovs)>0){
    for(i in 1:length(obj$indvCovs))
      param <- c(param, rep(0, obj$indvCovs[[i]]$get_param_length())) 
  }

  if(length(obj$teamCovs)>0){
    for(i in 1:length(obj$teamCovs)){
      param <- c(param, rep(0, obj$teamCovs[[i]]$get_param_length())) 
      
    }
  }
  return(param)
}

