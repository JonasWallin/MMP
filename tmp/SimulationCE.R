### Simulation
# 2021-06-22

library(parallel)
library(nlme)
library(robustbase)

## run this for delta1 = 0, -0.03, -0.05, -0.08,
## for each datatype = "HeCEM","HoCEM","HeCEM+GP","HoCEM+GP"



# setup

N <- 1000 # # of reps

no_cores <- detectCores() - 1 
cl<-makeCluster(no_cores)
clusterExport(cl,ls())
clusterEvalQ(cl, {
  library(MMP) 
  library(nlme)
  library(dplyr)})


system.time(res <- parLapply(cl,1:N, function(i,...) {
    set.seed(123+i); 
    c(ce_sim(l3n=20,l2n=7,l1n=7,
             mu0=3,mu1=0.01,
             #sg00=0.05,sg11=0.005,sg01=-0.005, # for HeCEM and HoCEM
             sp0=0.1,
             sp=0.1,
             se=0.1,
             sgp=0.05, # for "HeCEM+GP" and "HoCEM+GP"
             corr=0.05, # for "HeCEM+GP" and "HoCEM+GP"
             delta1=0,
             delta2=0.06, # for "HeCEM+GP" and "HoCEM+GP"
             datatype = "HeCEM+GP",
             REML=F,
             index = i,
             path = "/Users/baurne/Documents/Consensus emergence/SimulationResults/"),
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

saveRDS(result_list, file = "HeCEGP_0.RData")


