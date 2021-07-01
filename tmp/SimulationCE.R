### Simulation
# 2021-06-22

library(parallel)
library(nlme)


# # of reps
N <- 5

## test

# REML
system.time(res <- lapply(1:N, function(i,...) {
  set.seed(123+i); 
  c(ce_sim(l3n=20,l2n=7,l1n=7,
           mu0=3,mu1=0.01,
           sg00=0.05,sg11=0.005,sg01=-0.005,
           sp0=0.2,
           sp=0.2,
           se=0.2,
           delta1=-0.1,
           datatype = "HeCEM",
           REML=F),
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
                       mean = c(rowMeans(output[[i]][1:3,],na.rm = T),NA,NA),
                       NAs= c(rowMeans(is.na(output[[i]][1:3,])),NA,NA))
}

save(res, file = "res_5 <- .RData")
res_he <- res
output_he <- output

# ML
system.time(res <- lapply(1:N, function(i,...) {
  set.seed(123+i); 
  c(ce_sim(l3n=5,l2n=3,l1n=7,
           mu0=3,mu1=0.01,
           sg00=0.05,sg11=0.005,sg01=-0.005,
           sp0=0.2,
           sp=0.2,
           se=0.1,
           delta1=-0.1,
           datatype = "HoCEM",
           REML=F),
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
                       mean = c(rowMeans(output[[i]][1:3,],na.rm = T),NA,NA),
                       NAs= c(rowMeans(is.na(output[[i]][1:3,])),NA,NA))
}

## test
system.time(res <- lapply(1:N, function(i,...) {
  set.seed(123+i); 
  c(ce_sim(l3n=5,l2n=3,l1n=7,
           mu0=3,mu1=0.01,
           #sg00=0.05,sg11=0.005,sg01=-0.005,
           sp0=0.2,
           sp=0.2,
           se=0.2,
           sigma=sqrt(0.2),
           corr=0.05,
           delta1=-0.08,
           delta2 = -0.15,
           datatype = "HeCEM+GP"),
    seed =123+i)}  ))


# setup

N <- 100 # # of reps

no_cores <- detectCores() - 1 
cl<-makeCluster(no_cores)
clusterExport(cl,ls())
clusterEvalQ(cl, {
  library(MMP) 
  library(nlme)
  library(dplyr)})


system.time(res <- parLapply(cl,1:N, function(i,...) {
    set.seed(123+i); 
    c(ce_sim(l3n=8,l2n=3,l1n=3,
             mu0=3,mu1=0.01,
             sg00=0.1,sg11=0.005,sg01=-0.005,
             sp0=0.1,
             sp=0.1,
             se=0.1,
             delta1=-0.03,
             datatype = "HeCEM",
             REML=F),
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
                       mean = c(rowMeans(output[[i]][1:3,],na.rm = T),NA,NA),
                       NAs= c(rowMeans(is.na(output[[i]][1:3,])),NA,NA))
}




stopCluster(cl)


summary <- lapply(1:7, function(i) {
  output[[i]][1:3,c(N+1,N+2)]
})

names(summary) <- c("HeCEM",
                   "HoCEM",
                   "GP",
                   "HeCEMGP",
                   "HoCEMGP",
                   "GPGP",
                   "CEM")


saveRDS(output, file = "output_HeCE_d03.RData")
saveRDS(res, file = "res_HeCE_d03.RData")
saveRDS(summary, file = "summary_HeCE_d03.RData")

