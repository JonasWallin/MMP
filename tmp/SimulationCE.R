### Simulation
# 2021-06-22

library(parallel)


# # of reps
N <- 5

## test
system.time(res <- lapply(1:N, function(i,...) {
  set.seed(123+i); 
  c(ce_sim(l3n=5,l2n=3,l1n=7,
           mu0=3,mu1=0.01,
           sg00=0.05,sg11=0.005,sg01=-0.005,
           sp0=0.2,
           sp=0.2,
           se=0.2,
           delta1=-0.08,
           datatype = "HeCEM"),
    seed =123+i)}  ))

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
                       mean = c(rowMeans(output[[i]][1:3,],na.rm = T),NA),
                       NAs= c(rowMeans(is.na(output[[i]][1:3,])),NA))
}



## funkar inte att köra nedan med CEM
# setup
no_cores <- detectCores() - 1 
cl<-makeCluster(no_cores)
clusterExport(cl,ls())
clusterEvalQ(cl, library("MMP"))
#library(nlme)

system.time(res <- parLapply(cl,1:N, function(i,...) {
    set.seed(123+i); 
    c(ce_sim(l3n=3,l2n=3,l1n=7,
             mu0=3,mu1=0.01,
             sg00=0.05,sg11=0.005,sg01=-0.005,
             sp0=0.2,
             sp=0.2,
             se=0.2,
             #sigma = 0.01,
             #corr = 0.05,
             delta1=-0.08,
             #delta2 = 0.25,
             datatype = "HeCEM"),
      seed =123+i)}  ))

simresults <- cbind(delta1 = rowMeans(res[1:6,],na.rm=T), NAs=rowMeans(is.na(res[1:6,])))

stopCluster(cl)

