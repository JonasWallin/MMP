### Simulation


library(MASS)
library(dplyr)
library(MMP)

set.seed(2025)

########################
## Cohesion
########################
data("cohesiondat")

cohesiondat$y <- rowMeans(scale(cohesiondat[,4:6], center = FALSE))

# GP + GP
par.GPNL <- c(-4.2035533, -5.2085661 , 0.1726854 ,-2.9251160,
              -1.5413470 ,-4.3169783 , 0.0823050 ,-2.2834835 ,-1.1748932)
GPNL <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 | group, 
           emergence = ~ 1 , 
           method = "GP",
           method.team = "OU.homeostasis",
           time = "time",
           data =cohesiondat,
           param=par.GPNL)

# Esimated params
Dist <- as.matrix(dist(unique(GPNL$model$data$time)))
Table.coh <- c(GPNL$betas,  #betas
           exp(2 * GPNL$covariances$indv[[1]]), #sigma^2_v0
           exp(2 * GPNL$covariances$indv[[2]][2] - GPNL$covariances$indv[[2]][3])/2, #sigma^2_v1
           GPNL$covariances$indv[[2]][1], #delta_v
           exp( -GPNL$covariances$indv[[2]][3]), #kappa_v
           exp(-max(Dist)/exp( -GPNL$covariances$indv[[2]][3])), # corr = exp(-max(Dist)/k)
           exp(2 * GPNL$covariances$team[[1]]), # sigma^2_tau0
           exp(2 * GPNL$covariances$team[[2]][2]- GPNL$covariances$team[[2]][3])/2, # sigma^2_tau1
           GPNL$covariances$team[[2]][1], #delta_tau
           exp( -GPNL$covariances$team[[2]][3]), #kappa_tau
           exp(-max(Dist)/exp( -GPNL$covariances$team[[2]][3])), # corr = exp(-max(Dist)/k)
           exp(2 * GPNL$covariances$error[[1]]) #sigma2_eps
) 

param_names <- c("beta0", "beta1", 
                 "sigma^2_v0", "sigma^2_v1", "delta_v","kappa_v", "corr_v", 
                 "sigma^2_tau0", "sigma^2_tau1",  "delta_tau", "kappa_tau","corr_tau",
                 "sigma2_eps")
names(Table.coh) <- param_names
print(round(Table.coh,3)) 



# Simulation parameters
J <- 40       # groups
mj <- 6       # persons per group
Ti <- 4       # time points

mu0 <- 1.001573e+00     # intercept
mu1 <- -2.190427e-02    # slope

delta_group <- 8.047156e-02  # group GP scale change over time
delta_indiv <- 1.698632e-01  # individual GP scale change
sgp_group <- 1.700019e-02     # group GP variance at time 0
sgp_indiv <- 6.914247e-03     # individual GP variance at time 0
corr_group <- 3.963931e-01          # correlation at max distance    
corr_indiv <- 5.077924e-01          # correlation at max distance   
sigma2_e <- 2.173637e-04    # residual variance
sigma2_v0 <- 3.907732e-07     # indv random intercept variance
sigma2_tau0 <- 2.123306e-06     # group random intercept variance

  

# corr = exp(-maxDist/k)
# Time and distance
time <- 0:(Ti - 1)
Dist <- as.matrix(dist(time))
k_group <- max(Dist) / (-log(corr_group))  # length-scale parameter
k_indiv <- max(Dist) / (-log(corr_indiv))  # length-scale parameter

# --------------------------------
# 1. GROUP LEVEL GP
# --------------------------------
param_group <- c(log(sqrt((2 * 1/k_group) *sgp_group)), -log(k_group)) # exp(param[2] ) = 1/k OBS! verify this is correct
Sigma_group <- OUcov(Dist, param_group)
covs_group <- diag(exp(delta_group * time)) %*% Sigma_group %*% diag(exp(delta_group * time)) # sigma in mvrnorm
mean_group <- mu0 + mu1 * time

g_group <- MASS::mvrnorm(n = J, mu = mean_group, Sigma = covs_group)

gpdat_group <- data.frame(
  group = rep(1:J, each = Ti),
  time = rep(time, J),
  g_group = as.vector(t(g_group))
)

# --------------------------------
# 2. INDIVIDUAL LEVEL GP
# --------------------------------
total_indiv <- J * mj
param_indiv <- c(log(sqrt((2 * 1/k_indiv) *sgp_indiv)), -log(k_indiv)) # exp(param[2] ) = 1/k OBS! verify this is correct
Sigma_indiv <- OUcov(Dist, param_indiv)
covs_indiv <- diag(exp(delta_indiv * time)) %*% Sigma_indiv %*% diag(exp(delta_indiv * time))
mean_indiv <- rep(0, Ti)  # zero-centered GP

g_indiv <- MASS::mvrnorm(n = total_indiv, mu = mean_indiv, Sigma = covs_indiv)

gpdat_indiv <- data.frame(
  id_full = rep(1:total_indiv, each = Ti),
  time = rep(time, total_indiv),
  g_indiv = as.vector(t(g_indiv))
)

# --------------------------------
# 3. Design Matrix and Merging
# --------------------------------
dat <- expand.grid(time = time, person = 1:mj, group = 1:J) %>%
  arrange(group, person, time) %>%
  mutate(id_full = (group - 1) * mj + person)

dat <- dat %>%
  left_join(gpdat_group, by = c("group", "time")) %>%
  left_join(gpdat_indiv, by = c("id_full", "time"))

# --------------------------------
# 4. Add Random Intercepts and Residuals
# --------------------------------
v0ij_vals <- rnorm(n = total_indiv, mean = 0, sd = sqrt(sigma2_v0))
tau0j_vals <- rnorm(n = J, mean = 0, sd = sqrt(sigma2_tau0))

dat <- dat %>%
  mutate(
    v0ij = v0ij_vals[id_full],
    tau0j = tau0j_vals[group],
    eps = rnorm(n(), 0, sqrt(sigma2_e)),
    y = g_group + tau0j + v0ij + g_indiv + eps
  )

# Preview
#head(dat)

# GP
GP.h <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 | group, 
           emergence = ~ 1, 
           method = "GP",
           method.team = "OU.homeostasis", # byt namn till GP?
           time = "time",
           data = dat)

# Esimated params
Dist <- as.matrix(dist(unique(GP.h$model$data$time)))
Table.sim.coh <- c(GP.h$betas,  #betas
           exp(2 * GP.h$covariances$indv[[1]]), #sigma^2_v0
           exp(2 * GP.h$covariances$indv[[2]][2] - GP.h$covariances$indv[[2]][3])/2, #sigma^2_v1
           GP.h$covariances$indv[[2]][1], #delta_v
           exp( -GP.h$covariances$indv[[2]][3]), #kappa_v
           exp(-max(Dist)/exp( -GP.h$covariances$indv[[2]][3])), # corr = exp(-max(Dist)/k)
           exp(2 * GP.h$covariances$team[[1]]), # sigma^2_tau0
           exp(2 * GP.h$covariances$team[[2]][2]- GP.h$covariances$team[[2]][3])/2, # sigma^2_tau1
           GP.h$covariances$team[[2]][1], #delta_tau
           exp( -GP.h$covariances$team[[2]][3]), #kappa_tau
           exp(-max(Dist)/exp( -GP.h$covariances$team[[2]][3])), # corr = exp(-max(Dist)/k)
           exp(2 * GP.h$covariances$error[[1]]) #sigma2_eps
) 

param_names <- c("beta0", "beta1", 
                 "sigma^2_v0", "sigma^2_v1", "delta_v","kappa_v", "corr_v", 
                 "sigma^2_tau0", "sigma^2_tau1",  "delta_tau", "kappa_tau","corr_tau",
                 "sigma2_eps")
names(Table.sim.coh) <- param_names

print(round(cbind(Table.coh,Table.sim.coh),3)) 


#########################
##### Sherif
#########################
data("sherifdat")
sherifdat$time <- sherifdat$time + 1

# excluding last time point (individual measurement)
sherifdat <- subset(sherifdat, time <= 3)

# GP
GP.h.sherif <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 | group, 
           emergence = ~ 1, 
           method = "GP",
           method.team = "OU.homeostasis",
           time = "time",
           data = sherifdat)


# Esimated params
Dist <- as.matrix(dist(unique(GP.h.sherif$model$data$time)))
Table.sherif <- c(GP.h.sherif$betas,  #betas
           exp(2 * GP.h.sherif$covariances$indv[[1]]), #sigma^2_v0
           exp(2 * GP.h.sherif$covariances$indv[[2]][2] - GP.h.sherif$covariances$indv[[2]][3])/2, #sigma^2_v1
           GP.h.sherif$covariances$indv[[2]][1], #delta_v
           exp( -GP.h.sherif$covariances$indv[[2]][3]), #kappa_v
           exp(-max(Dist)/exp( -GP.h.sherif$covariances$indv[[2]][3])), # corr = exp(-max(Dist)/k)
           exp(2 * GP.h.sherif$covariances$team[[1]]), # sigma^2_tau0
           exp(2 * GP.h.sherif$covariances$team[[2]][2]- GP.h.sherif$covariances$team[[2]][3])/2, # sigma^2_tau1
           GP.h.sherif$covariances$team[[2]][1], #delta_tau
           exp( -GP.h.sherif$covariances$team[[2]][3]), #kappa_tau
           exp(-max(Dist)/exp( -GP.h.sherif$covariances$team[[2]][3])), # corr = exp(-max(Dist)/k)
           exp(2 * GP.h.sherif$covariances$error[[1]]) #sigma2_eps
) 

param_names <- c("beta0", "beta1", 
                 "sigma^2_v0", "sigma^2_v1", "delta_v","kappa_v", "corr_v", 
                 "sigma^2_tau0", "sigma^2_tau1",  "delta_tau", "kappa_tau","corr_tau",
                 "sigma2_eps")
names(Table.sherif) <- param_names
print(Table.sherif) 

set.seed(2025)


# Simulation parameters
J <- 8       # groups
mj <- 3       # persons per group
Ti <- 4       # time points

mu0 <- 3.321091e+00     # intercept
mu1 <- -2.089883e-01    # slope

delta_group <- -6.933332e-01  # group GP scale change over time
delta_indiv <- -1.082545e+00  # individual GP scale change
sgp_group <- 3.295941e+00     # group GP variance at time 0
sgp_indiv <- 3.199697e+00     # individual GP variance at time 0
corr_group <- 0.00000000000000001          # correlation at max distance    
corr_indiv <- 8.191625e-01         # correlation at max distance   
sigma2_e <- 2.505542e-03    # residual variance
sigma2_v0 <- 5.528073e-05     # indv random intercept variance
sigma2_tau0 <- 1.366200e+00     # group random intercept variance



# corr = exp(-maxDist/k)
# Time and distance
time <- 0:(Ti - 1)
Dist <- as.matrix(dist(time))
k_group <- max(Dist) / (-log(corr_group))  # length-scale parameter
k_indiv <- max(Dist) / (-log(corr_indiv))  # length-scale parameter

# --------------------------------
# 1. GROUP LEVEL GP
# --------------------------------
param_group <- c(log(sqrt((2 * 1/k_group) *sgp_group)), -log(k_group)) # exp(param[2] ) = 1/k OBS! verify this is correct
Sigma_group <- OUcov(Dist, param_group)
covs_group <- diag(exp(delta_group * time)) %*% Sigma_group %*% diag(exp(delta_group * time)) # sigma in mvrnorm
mean_group <- mu0 + mu1 * time

g_group <- MASS::mvrnorm(n = J, mu = mean_group, Sigma = covs_group)

gpdat_group <- data.frame(
  group = rep(1:J, each = Ti),
  time = rep(time, J),
  g_group = as.vector(t(g_group))
)

# --------------------------------
# 2. INDIVIDUAL LEVEL GP
# --------------------------------
total_indiv <- J * mj
param_indiv <- c(log(sqrt((2 * 1/k_indiv) *sgp_indiv)), -log(k_indiv)) # exp(param[2] ) = 1/k OBS! verify this is correct
Sigma_indiv <- OUcov(Dist, param_indiv)
covs_indiv <- diag(exp(delta_indiv * time)) %*% Sigma_indiv %*% diag(exp(delta_indiv * time))
mean_indiv <- rep(0, Ti)  # zero-centered GP

g_indiv <- MASS::mvrnorm(n = total_indiv, mu = mean_indiv, Sigma = covs_indiv)

gpdat_indiv <- data.frame(
  id_full = rep(1:total_indiv, each = Ti),
  time = rep(time, total_indiv),
  g_indiv = as.vector(t(g_indiv))
)

# --------------------------------
# 3. Design Matrix and Merging
# --------------------------------
dat <- expand.grid(time = time, person = 1:mj, group = 1:J) %>%
  arrange(group, person, time) %>%
  mutate(id_full = (group - 1) * mj + person)

dat <- dat %>%
  left_join(gpdat_group, by = c("group", "time")) %>%
  left_join(gpdat_indiv, by = c("id_full", "time"))

# --------------------------------
# 4. Add Random Intercepts and Residuals
# --------------------------------
v0ij_vals <- rnorm(n = total_indiv, mean = 0, sd = sqrt(sigma2_v0))
tau0j_vals <- rnorm(n = J, mean = 0, sd = sqrt(sigma2_tau0))

dat <- dat %>%
  mutate(
    v0ij = v0ij_vals[id_full],
    tau0j = tau0j_vals[group],
    eps = rnorm(n(), 0, sqrt(sigma2_e)),
    y = g_group + tau0j + v0ij + g_indiv + eps
  )

# Preview
#head(dat)

# GP
GP.h.sim <- ce(y ~ 1+time, 
           ~ 1 | person, 
           ~ 1 | group, 
           emergence = ~ 1, 
           method = "GP",
           method.team = "OU.homeostasis", # byt namn till GP?
           time = "time",
           data = dat)

# Esimated params
Dist <- as.matrix(dist(unique(GP.h.sim$model$data$time)))
Table.sim <- c(GP.h.sim$betas,  #betas
           exp(2 * GP.h.sim$covariances$indv[[1]]), #sigma^2_v0
           exp(2 * GP.h.sim$covariances$indv[[2]][2] - GP.h.sim$covariances$indv[[2]][3])/2, #sigma^2_v1
           GP.h.sim$covariances$indv[[2]][1], #delta_v
           exp( -GP.h.sim$covariances$indv[[2]][3]), #kappa_v
           exp(-max(Dist)/exp( -GP.h.sim$covariances$indv[[2]][3])), # corr = exp(-max(Dist)/k)
           exp(2 * GP.h.sim$covariances$team[[1]]), # sigma^2_tau0
           exp(2 * GP.h.sim$covariances$team[[2]][2]- GP.h.sim$covariances$team[[2]][3])/2, # sigma^2_tau1
           GP.h.sim$covariances$team[[2]][1], #delta_tau
           exp( -GP.h.sim$covariances$team[[2]][3]), #kappa_tau
           exp(-max(Dist)/exp( -GP.h.sim$covariances$team[[2]][3])), # corr = exp(-max(Dist)/k)
           exp(2 * GP.h.sim$covariances$error[[1]]) #sigma2_eps
) 

param_names <- c("beta0", "beta1", 
                 "sigma^2_v0", "sigma^2_v1", "delta_v","kappa_v", "corr_v", 
                 "sigma^2_tau0", "sigma^2_tau1",  "delta_tau", "kappa_tau","corr_tau",
                 "sigma2_eps")
names(Table.sim) <- param_names
print(round(cbind(Table,Table.sim),3)) 
