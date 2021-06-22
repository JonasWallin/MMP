### Simulation
# 2021-06-22

library(parallel)



# # of reps
N <- 3

## test
system.time(res <- sapply(1:N, function(i,...) {
  set.seed(123+i); 
  c(ce_sim(l3n=3,l2n=3,l1n=3,
           mu0=3,mu1=0.01,
           sg00=0.05,sg11=0.005,sg01=-0.005,
           sp0=0.01,
           sp=0.2,
           se=0.01,
           delta1=-0.8,
           datatype = "HeCEM"),
    seed =123+i)}  ))

simresults <- cbind(delta1 = rowMeans(res[1:6,],na.rm=T), NAs=rowMeans(is.na(res[1:6,])))


# setup
no_cores <- detectCores() - 1 
cl<-makeCluster(no_cores)
clusterExport(cl,ls())
clusterEvalQ(cl, library("MMP"))


system.time(res <- parSapply(cl,1:N, function(i,...) {
    set.seed(123+i); 
    c(ce_sim(l3n=3,l2n=3,l1n=3,
             mu0=3,mu1=0.01,
             sg00=0.05,sg11=0.005,sg01=-0.005,
             sp0=0.01,
             sp=0.2,
             se=0.01,
             #sigma = 0.01,
             #corr = 0.05,
             delta1=-0.8,
             #delta2 = 0.25,
             datatype = "HeCEM"),
      seed =123+i)}  ))

simresults <- cbind(delta1 = rowMeans(res[1:6,],na.rm=T), NAs=rowMeans(is.na(res[1:6,])))

stopCluster(cl)

