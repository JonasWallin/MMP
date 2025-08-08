# install.packages(c("MASS","dplyr","tidyr","ggplot2","viridis","patchwork"))
library(MASS)       # for mvrnorm
library(dplyr)      # for data wrangling
library(tidyr)      # for pivoting
library(ggplot2)    # for plotting
library(viridis)    # for color scales
library(patchwork)  # for arranging plots
set.seed(4)
# Helper to sample from (possibly low-rank) covariance Σ
sample_gp_paths <- function(Sigma, n.samps = 5, tol = 1e-8) {
  ev  <- eigen(Sigma, symmetric = TRUE)
  pos <- which(ev$values > tol)
  L   <- ev$vectors[,pos] %*% diag(sqrt(ev$values[pos]), length(pos))
  Z   <- matrix(rnorm(length(pos) * n.samps), nrow = length(pos))
  t(L %*% Z)
}

# Covariance kernelsm
cov_intBM <- function(s, t) {
  smin <- outer(s, t, pmin); smax <- outer(s, t, pmax)
  (smin^2 * smax)/2 - (smin^3)/6
}
cov_intBM_scaled <- function(s, t) {
  K0      <- cov_intBM(s, t)
  scaling <- 3 / sqrt(outer(s, t, "*"))
  K0 * scaling + outer(rep(1,length(s)),rep(1,length(t)),"*")
}
cov_rank1 <- function(s, t) { outer(s, t, "*") + outer(rep(1,length(s)),rep(1,length(t)),"*")}

# Time grids + observations
times_disc <- 1:5
times_cont <- seq(1, 5, length.out = 100)
obs_df     <- tibble(time = 1:4, y = c(0.5, -1.0, 1.5, 0.0))

# Function to build all four plots for a given kernel
make_plots <- function(kernel, name,
                       times_disc, times_cont, obs_df, n.samps = 5,
                       tau_name = expression(tau),
                       base_size =26) {
  # 1) Covariances
  Sigma_cont <- kernel(times_cont, times_cont)
  Sigma_disc <- kernel(times_disc, times_disc)
  K_xX       <- kernel(times_cont, obs_df$time)           # 100×4
  K_XX       <- kernel(obs_df$time, obs_df$time) +
    diag(1e-1, nrow(obs_df))                  # 4×4
  K_inv      <- solve(K_XX)
  
  # 2) PRIOR
  prior_paths_cont <- sample_gp_paths(Sigma_cont, n.samps)
  # sample‐cloud at discrete times
  idxs <- sapply(times_disc, function(t) which.min(abs(times_cont - t)))
  prior_paths_disc <- prior_paths_cont[, idxs]
  prior_mv_disc <- tibble(
    time  = times_disc,
    mean  = 0,
    sd    = sqrt(diag(Sigma_disc))
  ) %>% mutate(lower = mean - 1.96*sd, upper = mean + 1.96*sd)
  prior_mv_cont <- tibble(
    time  = times_cont,
    mean  = 0,
    sd    = sqrt(diag(Sigma_cont))
  ) %>% mutate(lower = mean - 1.96*sd, upper = mean + 1.96*sd)
  
  # 3) POSTERIOR (continuous)
  post_mean_cont <- K_xX %*% (K_inv %*% obs_df$y)                      # 100×1
  post_cov_cont  <- Sigma_cont - K_xX %*% K_inv %*% t(K_xX)            # 100×100
  post_mv_cont   <- tibble(
    time  = times_cont,
    mean  = as.vector(post_mean_cont),
    sd    = sqrt(diag(post_cov_cont))
  ) %>% mutate(lower = mean - 1.96*sd, upper = mean + 1.96*sd)
  post_paths_cont <- sweep(
    sample_gp_paths(post_cov_cont, n.samps),
    2, post_mean_cont, FUN = "+"
  )
  
  # 4) POSTERIOR (discrete)
  K_dX            <- kernel(times_disc, obs_df$time)                  # 5×4
  K_Xd            <- kernel(obs_df$time, times_disc)                  # 4×5
  post_mean_disc  <- as.vector(K_dX %*% (K_inv %*% obs_df$y))         # 5×1
  post_cov_disc   <- Sigma_disc - K_dX %*% K_inv %*% K_Xd             # 5×5
  post_mv_disc <- tibble(
    time  = times_disc,
    mean  = post_mean_disc,
    sd    = sqrt(diag(post_cov_disc))
  ) %>% mutate(lower = mean - 1.96*sd, upper = mean + 1.96*sd)
  post_paths_disc <- post_paths_cont[, idxs]
  
  # reshape clouds
  df_prior_cloud_disc <- as_tibble(prior_paths_disc) %>%
    mutate(sample = row_number()) %>%
    pivot_longer(-sample, names_to="i", values_to="y") %>%
    mutate(time = times_disc[as.integer(sub("V","",i))])
  df_post_cloud_disc  <- as_tibble(post_paths_disc)  %>%
    mutate(sample = row_number()) %>%
    pivot_longer(-sample, names_to="i", values_to="y") %>%
    mutate(time = times_disc[as.integer(sub("V","",i))])
  
  df_prior_paths_cont <- as_tibble(prior_paths_cont) %>%
    mutate(sample = row_number()) %>%
    pivot_longer(-sample, names_to="i", values_to="y") %>%
    mutate(time = times_cont[as.integer(sub("V","",i))])
  df_post_paths_cont  <- as_tibble(post_paths_cont)  %>%
    mutate(sample = row_number()) %>%
    pivot_longer(-sample, names_to="i", values_to="y") %>%
    mutate(time = times_cont[as.integer(sub("V","",i))])
  
  # 5) Build ggplots
  p_disc_prior <- ggplot() +
    geom_errorbar(data = prior_mv_disc,
                  aes(x = time, ymin = lower, ymax = upper),
                  width = 0.2, color = "grey40") +
    geom_point(data = df_prior_cloud_disc,
               aes(x = time, y = y, color = factor(sample)),
               alpha = 0.7, size = 2) +
    geom_point(data = prior_mv_disc,
               aes(x = time, y = mean),
               shape = 4, size = 3, stroke = 1.2, color = "black") +
    labs(
      x     = "Time",
      y     = tau_name,
      color = NULL      # no legend title
    ) +
    scale_color_viridis_d(guide = "none") +
    scale_x_continuous(breaks = times_disc) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title.y = element_text(
        size  = base_size, 
      ),
      axis.title.x = element_text(
        size  = base_size, 
      )
    )
  
  p_disc_post <- ggplot() +
    geom_errorbar(data = post_mv_disc,
                  aes(x = time, ymin = lower, ymax = upper),
                  width = 0.2, color = "grey40") +
    geom_point(data = df_post_cloud_disc,
               aes(x = time, y = y, color = factor(sample)),
               alpha = 0.7, size = 2) +
    geom_point(data = post_mv_disc,
               aes(x = time, y = mean),
               shape = 4, size = 3, stroke = 1.2, color = "black") +
    geom_point(data = obs_df,
               aes(x = time, y = y),
               shape = 21, fill = "black", color = "white", size = 4) +
    labs(
      x     = "Time",
      y     = tau_name,
      color = NULL
    ) +
    scale_color_viridis_d(guide = "none") +
    scale_x_continuous(breaks = times_disc) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title.y = element_text(
        size  = base_size, 
      ),
      axis.title.x = element_text(
        size  = base_size, 
      )
    )
  
  p_cont_prior <- ggplot() +
    geom_ribbon(data = prior_mv_cont,
                aes(time, ymin = lower, ymax = upper),
                fill = "grey80") +
    geom_line(data = df_prior_paths_cont,
              aes(time, y = y, color = factor(sample)),
              alpha = 0.8) +
    labs(
      x     = "t",
      y     = "f(t)",
      color = NULL
    ) +
    scale_color_viridis_d(guide = "none") +
    theme_minimal(base_size = 14) +
    theme(
      axis.title.y = element_text(
        size  = base_size, 
      ),
      axis.title.x = element_text(
        size  = base_size, 
      )
    )
  
  p_cont_post <- ggplot() +
    geom_ribbon(data = post_mv_cont,
                aes(time, ymin = lower, ymax = upper),
                fill = "grey80") +
    geom_line(data = df_post_paths_cont,
              aes(time, y = y, color = factor(sample)),
              alpha = 0.8) +
    geom_line(data = post_mv_cont,
              aes(time, y = mean),
              size = 1, color = "black") +
    geom_point(data = obs_df,
               aes(x = time, y = y),
               shape = 21, fill = "black", color = "white", size = 4) +
    labs(
      x     = "t",
      y     = "f(t)",
      color = NULL
    ) +
    scale_color_viridis_d(guide = "none") +
    theme_minimal(base_size = 14) +
    theme(
      axis.title.y = element_text(
        size  = base_size, 
      ),
      axis.title.x = element_text(
        size  = base_size, 
      )
    )
  
  list(
    disc_prior     = p_disc_prior,
    disc_posterior = p_disc_post,
    cont_prior     = p_cont_prior,
    cont_posterior = p_cont_post
  )
}

# 6) Generate and arrange for both kernels
intbm_plots  <- make_plots(cov_intBM_scaled, "Scaled IntBM",
                           times_disc, times_cont, obs_df)
rank1_plots <- make_plots(cov_rank1,       "Rank-1 Kernel",
                          times_disc, times_cont, obs_df,tau_name = expression(tau[0] + tau[1] %.% t))

# Arrange in a 4×2 grid: each kernel’s discrete and continuous plots side by side
layout <- (intbm_plots$disc_prior    + intbm_plots$disc_posterior) /
  (intbm_plots$cont_prior    + intbm_plots$cont_posterior) /
  (rank1_plots$disc_prior   + rank1_plots$disc_posterior) /
  (rank1_plots$cont_prior   + rank1_plots$cont_posterior) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print((intbm_plots$disc_prior  + intbm_plots$cont_prior))
print((rank1_plots$disc_prior  + rank1_plots$cont_prior))
print((intbm_plots$disc_posterior  + intbm_plots$cont_posterior))
print((rank1_plots$disc_posterior  + rank1_plots$cont_posterior))

ggsave("disc_prior_bm.pdf", intbm_plots$disc_prior, width = 12, height = 12)
ggsave("disc_prior_r1.pdf", rank1_plots$disc_prior, width = 12, height = 12)
ggsave("disc_post_bm.pdf", intbm_plots$disc_posterior, width = 12, height = 12)
ggsave("disc_post_r1.pdf", rank1_plots$disc_posterior, width = 12, height = 12)

ggsave("cont_prior_bm.pdf", intbm_plots$cont_prior, width = 12, height = 12)
ggsave("cont_prior_r1.pdf", rank1_plots$cont_prior, width = 12, height = 12)
ggsave("cont_post_bm.pdf", intbm_plots$cont_posterior, width = 12, height = 12)
ggsave("cont_post_r1.pdf", rank1_plots$cont_posterior, width = 12, height = 12)
