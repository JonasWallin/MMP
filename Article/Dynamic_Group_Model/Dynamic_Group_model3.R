# --- Packages -----------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(mvtnorm)

# --- Your sampling routine ----------------------------------------------------
sample_gp_paths <- function(Sigma, n.samps = 5, tol = 1e-8) {
  ev  <- eigen(Sigma, symmetric = TRUE)
  pos <- which(ev$values > tol)
  L   <- ev$vectors[,pos] %*% diag(sqrt(ev$values[pos]), length(pos))
  Z   <- matrix(rnorm(length(pos) * n.samps), nrow = length(pos))
  t(L %*% Z)
}

# --- Your covariance kernels --------------------------------------------------
cov_intBM <- function(s, t) {
  smin <- outer(s, t, pmin); smax <- outer(s, t, pmax)
  (smin^2 * smax)/2 - (smin^3)/6
}
cov_intBM_scaled <- function(s, t) {
  K0      <- cov_intBM(s, t)
  scaling <- 3 / sqrt(outer(s, t, "*"))
  K0 * scaling + outer(rep(1,length(s)),rep(1,length(t)),"*")
}
cov_rank1 <- function(s, t) { outer(s, t, "*") + outer(rep(1,length(s)),rep(1,length(t)),"*") }

# --- Helper: GP conditioning (works with any of the kernels above) ------------
gp_condition <- function(t_obs, y_obs, t_pred, kernel_fun, sigma_eps) {
  K_oo <- kernel_fun(t_obs, t_obs)
  K_op <- kernel_fun(t_obs, t_pred)
  K_pp <- kernel_fun(t_pred, t_pred)
  
  Sigma_y <- K_oo + diag(sigma_eps^2, length(t_obs))
  # mean = K_po * Sigma_y^{-1} * y
  alpha <- solve(Sigma_y, y_obs)
  mean_pred <- as.vector(t(K_op) %*% alpha)
  
  # cov = K_pp - K_po * Sigma_y^{-1} * K_op
  V <- solve(Sigma_y, K_op)
  cov_pred <- K_pp - t(K_op) %*% V
  list(mean = mean_pred, cov = cov_pred)
}
## ---------------------------------------------------------------------------
## 1) SIMULATION ONLY: generate one dataset (no fitting)
## ---------------------------------------------------------------------------
simulate_gp_dataset <- function(groups = c("Group 1","Group 2","Group 3"),
                                t_obs = 1:5, t_pred = 6:10,
                                kernel_gen = cov_intBM_scaled,  # kernel used for data gen
                                sigma_eps = 1.0,
                                seed = 42) {
  set.seed(seed)
  t_all <- c(t_obs, t_pred)
  G <- length(groups)
  
  obs_df   <- vector("list", G)
  truth_df <- vector("list", G)  # optional: latent f and noisy y on all times
  
  for (g in seq_len(G)) {
    K_all <- kernel_gen(t_all, t_all)
    f_all <- as.numeric(sample_gp_paths(K_all, n.samps = 1))  # latent path
    y_all <- f_all + rnorm(length(t_all), 0, sigma_eps)
    
    idx_obs <- match(t_obs, t_all)
    y_obs   <- y_all[idx_obs]
    
    obs_df[[g]] <- data.frame(group = groups[g], time = t_obs, y = y_obs)
    truth_df[[g]] <- data.frame(group = groups[g], time = t_all,
                                f = f_all, y_noisy = y_all)
  }
  
  list(
    obs   = do.call(rbind, obs_df),
    truth = do.call(rbind, truth_df),
    t_obs = t_obs,
    t_pred = t_pred
  )
}

## ---------------------------------------------------------------------------
## 2) FITTING ONLY: given observed data, produce preds and trajectories
## ---------------------------------------------------------------------------
fit_gp_by_group <- function(obs_df,
                            t_pred,
                            t_grid,
                            kernel_fun = cov_intBM_scaled,
                            sigma_eps = 1.0) {
  groups <- sort(unique(obs_df$group))
  pred_df <- vector("list", length(groups))
  grid_df <- vector("list", length(groups))
  traj_df <- vector("list", length(groups))
  traj_pred_df <- vector("list", length(groups))  # same draw as traj_df
  
  for (i in seq_along(groups)) {
    g <- groups[i]
    dat_g <- obs_df[obs_df$group == g, , drop = FALSE]
    dat_g <- dat_g[order(dat_g$time), ]
    t_obs <- dat_g$time
    y_obs <- dat_g$y
    
    ## (A) Predictions at discrete times (for means/CI)
    fit_pred  <- gp_condition(t_obs, y_obs, t_pred, kernel_fun, sigma_eps)
    pred_mean <- fit_pred$mean
    pred_sd   <- sqrt(pmax(diag(fit_pred$cov), 0))
    
    ## (B) Fine-grid posterior (for ribbons/mean)
    fit_grid  <- gp_condition(t_obs, y_obs, t_grid, kernel_fun, sigma_eps)
    grid_mean <- fit_grid$mean
    grid_sd   <- sqrt(pmax(diag(fit_grid$cov), 0))
    
    ## (C) ONE joint posterior draw over union of times, then split
    t_joint <- unique(c(t_grid, t_pred))  # avoids duplicate times -> stable numerics
    fit_joint <- gp_condition(t_obs, y_obs, t_joint, kernel_fun, sigma_eps)
    mu_joint  <- fit_joint$mean
    cov_joint <- fit_joint$cov
    draw_joint <- as.numeric(sample_gp_paths(cov_joint, n.samps = 1) + matrix(mu_joint, nrow = 1))
    
    # map joint draw back to grid and pred times
    idx_grid <- match(t_grid, t_joint)
    idx_pred <- match(t_pred, t_joint)
    traj_draw_grid <- draw_joint[idx_grid]
    traj_draw_pred <- draw_joint[idx_pred]
    
    pred_df[[i]] <- data.frame(
      group = g, time = t_pred,
      mean = pred_mean, sd = pred_sd,
      lower = pred_mean - 1.96 * pred_sd,
      upper = pred_mean + 1.96 * pred_sd
    )
    grid_df[[i]] <- data.frame(
      group = g, time = t_grid,
      mean = grid_mean, sd = grid_sd,
      lower = grid_mean - 1.96 * grid_sd,
      upper = grid_mean + 1.96 * grid_sd
    )
    traj_df[[i]]      <- data.frame(group = g, time = t_grid, y = traj_draw_grid)
    traj_pred_df[[i]] <- data.frame(group = g, time = t_pred, y = traj_draw_pred)
  }
  
  list(
    pred       = do.call(rbind, pred_df),
    grid       = do.call(rbind, grid_df),
    traj_grid  = do.call(rbind, traj_df),
    traj_pred  = do.call(rbind, traj_pred_df)  # crosses at 6,7,8 in your discrete plot
  )
}

## ---------------------------------------------------------------------------
## Example usage
## ---------------------------------------------------------------------------
# Simulate once (generation kernel can be chosen independently of fitting kernel)
sim <- simulate_gp_dataset(kernel_gen = cov_intBM_scaled,
                           t_obs = 1:5, t_pred = 6:8,
                           sigma_eps = .5, seed = 40)

# Choose fine grid for the WITH-trajectory plot
t_grid <- seq(min(sim$t_obs), max(sim$t_pred), by = 0.1)

# Fit (can try cov_intBM_scaled or cov_rank1 here)
fit_scaled <- fit_gp_by_group(sim$obs, t_pred = sim$t_pred, t_grid = t_grid,
                              kernel_fun = cov_intBM_scaled, sigma_eps = .5)
fit_rank1 <- fit_gp_by_group(sim$obs, t_pred = sim$t_pred, t_grid = t_grid,
                              kernel_fun = cov_rank1, sigma_eps = .5)


# --- Build plots for a given fit object ---------------------------------------
make_plots <- function(sim, fit, kernel_label,
                       group_cols   = c("Group 1"="#e41a1c","Group 2"="#4daf4a","Group 3"="#377eb8"),
                       group_shapes = c("Group 1"=16, "Group 2"=17, "Group 3"=15),
                       times_show   = c(6,7,8)) {
  
  # ----- NO trajectories (discrete at times_show) -----
  pred_show  <- subset(fit$pred, time %in% times_show)
  trajp_show <- subset(fit$traj_pred, time %in% times_show)
  # ----- WITH trajectories (fine grid) -----
  p_with_traj <-
    ggplot() +
    geom_point(data = sim$obs,
               aes(x = time, y = y, shape = group),
               color = "black", size = 2, show.legend = TRUE) +
    geom_ribbon(data = fit$grid,
                aes(x = time, ymin = lower, ymax = upper, fill = group),
                alpha = 0.20, show.legend = TRUE) +
    geom_path(data = fit$traj_grid,
              aes(x = time, y = y, color = group, group = group),
              linewidth = 0.5, alpha = 0.85, show.legend = TRUE) +
    
    geom_point(data = trajp_show,
               aes(x = time, y = y, color = group, group = group),
               shape = 4, size = .7, stroke = 1.1, show.legend = FALSE) +
    scale_color_manual(values = group_cols, name = "group") +
    scale_fill_manual(values = group_cols,  name = "group") +
    scale_shape_manual(values = group_shapes, name = "group") +
    # keep the old single "group" legend, but show black shapes clearly
    guides(shape = guide_legend(override.aes = list(color = "black", size = 3))) +
    labs(x = "Time", y = expression(y/f(t))) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right", panel.grid.minor = element_blank())
  
  dodge <- position_dodge2(width = 0.5, preserve = "single")
  p_no_traj <-
    ggplot() +
    geom_point(data = sim$obs,
               aes(x = time, y = y, shape = group),
               color = "black", size = 2, show.legend = TRUE) +
    # Error bars with clearly visible caps
    geom_errorbar(data = pred_show,
                  aes(x = time, ymin = lower, ymax = upper, color = group, group = group),
                  width = .7,  # cap width (make bigger if needed)
                  size  = .5,   # line thickness (use size for compatibility)
                  #position = dodge,
                  show.legend = TRUE) +
    # Single posterior draw shown as X crosses at 6,7,8
    geom_point(data = trajp_show,
               aes(x = time, y = y, color = group, group = group),
               shape = 4, size = .7, stroke = 1.1, show.legend = FALSE) +
    scale_color_manual(values = group_cols, name = "group") +
    scale_fill_manual(values = group_cols,  name = "group") +
    scale_shape_manual(values = group_shapes, name = "group") +
    guides(shape = guide_legend(override.aes = list(color = "black", size = 3))) +
    labs(x = "Time", y = expression(y/tau)) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right", panel.grid.minor = element_blank())
  
  list(with_traj = p_with_traj, no_traj = p_no_traj)
}
# --- Create figures for BOTH kernels ------------------------------------------
plots_scaled <- make_plots(sim, fit_scaled, kernel_label = "scaled kernel")
plots_rank1  <- make_plots(sim, fit_rank1,  kernel_label = "rank-1 kernel")

# Show all four
print(plots_scaled$with_traj)
ggsave("scaled_with_traj.png", plots_scaled$with_traj, width = 8, height = 6)
print(plots_scaled$no_traj)
ggsave("scaled_no_traj.png", plots_scaled$no_traj, width = 8, height = 6)
print(plots_rank1$with_traj)
ggsave("rank1_with_traj.png", plots_rank1$with_traj, width = 8, height = 6)
print(plots_rank1$no_traj)
ggsave("rank1_no_traj.png", plots_rank1$no_traj, width = 8, height = 6)
