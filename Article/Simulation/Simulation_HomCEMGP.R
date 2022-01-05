### Simulation HomCEM GP
# 2021-06-22

# Note: running this simulation may take several weeks, 
# depending on your computer. 

## Rename:
# datatypes:  HoCEM = HomCEM,
#             HeCEM = HomCEM,
#             HoCEM+GP = HomCEM+GP,
#             HeCEM+GP = HetCEM+GP

library(parallel)
library(nlme)
library(robustbase)

# setup

N <- 1000 # # of reps
no_cores <- detectCores() - 1 
cl<-makeCluster(no_cores)

# Specify a path for saving all output
mypath <- "/Users/baurne/Documents/Consensus emergence/SimulationResults/HomCEM_GP"

######

## HomCEM GP
# delta1 = 0

clusterExport(cl,ls())
clusterEvalQ(cl, {
  library(MMP) 
  library(nlme)
  library(dplyr)})


system.time(res <- parLapply(cl,1:N, function(i,...) {
  set.seed(123+i); 
  c(ce_sim(l3n=20,l2n=7,l1n=7,
           mu0=3,mu1=0.01,
           sp0=0.1,
           sp=0.1,
           se=0.1,
           sgp=0.05, 
           corr=0.05, 
           delta1=0,
           delta2=0.06, 
           datatype = "HoCEM+GP",
           REML=F,
           index = i,
           path = mypath),
    seed =123+i)}  ))

stopCluster(cl)

output <- sapply(1:8,function(i) sapply(res, `[[`, i))
names(output) <- c("HeCEM",
                   "HoCEM",
                   "GP",
                   "HeCEMGP",
                   "HoCEMGP",
                   "GPGP",
                   "CEM",
                   "seed")

for (i in 1:7) {
  output[[i]] <- cbind(rbind(output[[i]],seed =output[[8]]),
                       median = c(rowMedians(output[[i]][1:3,1:1000],na.rm = T),NA,NA),
                       qdev = c((quantile(output[[i]][1,1:1000], .66,na.rm = T)- 
                                   quantile(output[[i]][1,1:1000], .33,na.rm = T))/2,
                                (quantile(output[[i]][2,1:1000], .66,na.rm = T)- 
                                   quantile(output[[i]][2,1:1000], .33,na.rm = T))/2,
                                (quantile(output[[i]][3,1:1000], .66,na.rm = T)- 
                                   quantile(output[[i]][3,1:1000], .33,na.rm = T))/2,NA,NA),
                       NAs= c(rowMeans(is.na(output[[i]][1:3,])),NA,NA))
}



sumres <- lapply(1:7, function(i) {
  output[[i]][1:3,c(N+1,N+2,N+3)]
})

names(sumres) <- c("HeCEM",
                   "HoCEM",
                   "GP",
                   "HeCEMGP",
                   "HoCEMGP",
                   "GPGP",
                   "CEM")

result_list <- list(ce_sim_output = res,
                    summary_output = output,
                    mean_results = sumres)

saveRDS(result_list, file = "HomCEGP_delta0.RData")

###########

# delta1 = -0.03

clusterExport(cl,ls())
clusterEvalQ(cl, {
  library(MMP) 
  library(nlme)
  library(dplyr)})


system.time(res <- parLapply(cl,1:N, function(i,...) {
  set.seed(123+i); 
  c(ce_sim(l3n=20,l2n=7,l1n=7,
           mu0=3,mu1=0.01,
           sp0=0.1,
           sp=0.1,
           se=0.1,
           sgp=0.05, 
           corr=0.05, 
           delta1=-0.03,
           delta2=0.06, 
           datatype = "HoCEM+GP",
           REML=F,
           index = i,
           path = mypath),
    seed =123+i)}  ))

stopCluster(cl)

output <- sapply(1:8,function(i) sapply(res, `[[`, i))
names(output) <- c("HeCEM",
                   "HoCEM",
                   "GP",
                   "HeCEMGP",
                   "HoCEMGP",
                   "GPGP",
                   "CEM",
                   "seed")

for (i in 1:7) {
  output[[i]] <- cbind(rbind(output[[i]],seed =output[[8]]),
                       median = c(rowMedians(output[[i]][1:3,1:1000],na.rm = T),NA,NA),
                       qdev = c((quantile(output[[i]][1,1:1000], .66,na.rm = T)- 
                                   quantile(output[[i]][1,1:1000], .33,na.rm = T))/2,
                                (quantile(output[[i]][2,1:1000], .66,na.rm = T)- 
                                   quantile(output[[i]][2,1:1000], .33,na.rm = T))/2,
                                (quantile(output[[i]][3,1:1000], .66,na.rm = T)- 
                                   quantile(output[[i]][3,1:1000], .33,na.rm = T))/2,NA,NA),
                       NAs= c(rowMeans(is.na(output[[i]][1:3,])),NA,NA))
}



sumres <- lapply(1:7, function(i) {
  output[[i]][1:3,c(N+1,N+2,N+3)]
})

names(sumres) <- c("HeCEM",
                   "HoCEM",
                   "GP",
                   "HeCEMGP",
                   "HoCEMGP",
                   "GPGP",
                   "CEM")

result_list <- list(ce_sim_output = res,
                    summary_output = output,
                    mean_results = sumres)

saveRDS(result_list, file = "HomCEGP_delta-0.03.RData")

###########

# delta1 = -0.05

clusterExport(cl,ls())
clusterEvalQ(cl, {
  library(MMP) 
  library(nlme)
  library(dplyr)})


system.time(res <- parLapply(cl,1:N, function(i,...) {
  set.seed(123+i); 
  c(ce_sim(l3n=20,l2n=7,l1n=7,
           mu0=3,mu1=0.01,
           sp0=0.1,
           sp=0.1,
           se=0.1,
           sgp=0.05, 
           corr=0.05, 
           delta1=-0.05,
           delta2=0.06, 
           datatype = "HoCEM+GP",
           REML=F,
           index = i,
           path = mypath),
    seed =123+i)}  ))

stopCluster(cl)

output <- sapply(1:8,function(i) sapply(res, `[[`, i))
names(output) <- c("HeCEM",
                   "HoCEM",
                   "GP",
                   "HeCEMGP",
                   "HoCEMGP",
                   "GPGP",
                   "CEM",
                   "seed")

for (i in 1:7) {
  output[[i]] <- cbind(rbind(output[[i]],seed =output[[8]]),
                       median = c(rowMedians(output[[i]][1:3,1:1000],na.rm = T),NA,NA),
                       qdev = c((quantile(output[[i]][1,1:1000], .66,na.rm = T)- 
                                   quantile(output[[i]][1,1:1000], .33,na.rm = T))/2,
                                (quantile(output[[i]][2,1:1000], .66,na.rm = T)- 
                                   quantile(output[[i]][2,1:1000], .33,na.rm = T))/2,
                                (quantile(output[[i]][3,1:1000], .66,na.rm = T)- 
                                   quantile(output[[i]][3,1:1000], .33,na.rm = T))/2,NA,NA),
                       NAs= c(rowMeans(is.na(output[[i]][1:3,])),NA,NA))
}



sumres <- lapply(1:7, function(i) {
  output[[i]][1:3,c(N+1,N+2,N+3)]
})

names(sumres) <- c("HeCEM",
                   "HoCEM",
                   "GP",
                   "HeCEMGP",
                   "HoCEMGP",
                   "GPGP",
                   "CEM")

result_list <- list(ce_sim_output = res,
                    summary_output = output,
                    mean_results = sumres)

saveRDS(result_list, file = "HomCEGP_delta-0.05.RData")

###########

# delta1 = -0.08

clusterExport(cl,ls())
clusterEvalQ(cl, {
  library(MMP) 
  library(nlme)
  library(dplyr)})


system.time(res <- parLapply(cl,1:N, function(i,...) {
  set.seed(123+i); 
  c(ce_sim(l3n=20,l2n=7,l1n=7,
           mu0=3,mu1=0.01,
           sp0=0.1,
           sp=0.1,
           se=0.1,
           sgp=0.05, 
           corr=0.05, 
           delta1=-0.08,
           delta2=0.06, 
           datatype = "HoCEM+GP",
           REML=F,
           index = i,
           path = mypath),
    seed =123+i)}  ))

stopCluster(cl)

output <- sapply(1:8,function(i) sapply(res, `[[`, i))
names(output) <- c("HeCEM",
                   "HoCEM",
                   "GP",
                   "HeCEMGP",
                   "HoCEMGP",
                   "GPGP",
                   "CEM",
                   "seed")

for (i in 1:7) {
  output[[i]] <- cbind(rbind(output[[i]],seed =output[[8]]),
                       median = c(rowMedians(output[[i]][1:3,1:1000],na.rm = T),NA,NA),
                       qdev = c((quantile(output[[i]][1,1:1000], .66,na.rm = T)- 
                                   quantile(output[[i]][1,1:1000], .33,na.rm = T))/2,
                                (quantile(output[[i]][2,1:1000], .66,na.rm = T)- 
                                   quantile(output[[i]][2,1:1000], .33,na.rm = T))/2,
                                (quantile(output[[i]][3,1:1000], .66,na.rm = T)- 
                                   quantile(output[[i]][3,1:1000], .33,na.rm = T))/2,NA,NA),
                       NAs= c(rowMeans(is.na(output[[i]][1:3,])),NA,NA))
}



sumres <- lapply(1:7, function(i) {
  output[[i]][1:3,c(N+1,N+2,N+3)]
})

names(sumres) <- c("HeCEM",
                   "HoCEM",
                   "GP",
                   "HeCEMGP",
                   "HoCEMGP",
                   "GPGP",
                   "CEM")

result_list <- list(ce_sim_output = res,
                    summary_output = output,
                    mean_results = sumres)

saveRDS(result_list, file = "HomCEGP_delta-0.08.RData")

###########

