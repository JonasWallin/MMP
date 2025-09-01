# Packages used (no library() required if you keep full namespaces):
# dplyr, tidyr, rlang, stats

lag1_rho <- function(S, times) {
  o <- order(times); t <- times[o]; S <- S[o, o]
  if (length(t) < 2) return(0.8)
  i <- 1:(length(t) - 1); j <- 2:length(t)
  num <- S[cbind(i, j)]
  den <- sqrt(diag(S)[i] * diag(S)[j])
  r  <- num / pmax(den, 1e-12)
  r <- r[is.finite(r)]
  if (!length(r)) return(0.8)
  stats::median(pmin(pmax(r, 1e-4), 0.999))
}

empirical_indiv_cov <- function(dat, time_col = "time", group_col = "group",
                                id_col = "id_full", y_col = "y") {
  df <- dat |>
    dplyr::rename(time = !!time_col, group = !!group_col, id = !!id_col, y = !!y_col)
  
  gbar <- df |>
    dplyr::group_by(time, group) |>
    dplyr::summarise(ybar = mean(y, na.rm = TRUE), .groups = "drop")
  
  resid_df <- df |>
    dplyr::left_join(gbar, by = c("time", "group")) |>
    dplyr::mutate(u = y - ybar)
  
  W <- resid_df |>
    dplyr::select(id, time, u) |>
    tidyr::pivot_wider(names_from = time, values_from = u) |>
    dplyr::arrange(id)
  
  M <- as.matrix(W[, -1, drop = FALSE]); storage.mode(M) <- "double"
  S <- stats::cov(M, use = "pairwise.complete.obs")
  times <- as.numeric(colnames(W)[-1])
  list(S = S, times = times)
}

## -- Empirical team covariance from group means --------------------------------
empirical_team_cov <- function(dat, time_col = "time", group_col = "group", y_col = "y") {
  df <- dat |>
    dplyr::rename(time = !!time_col, group = !!group_col, y = !!y_col)
  
  gbar <- df |>
    dplyr::group_by(time, group) |>
    dplyr::summarise(ybar = mean(y, na.rm = TRUE), .groups = "drop")
  
  W <- gbar |>
    tidyr::pivot_wider(names_from = time, values_from = ybar) |>
    dplyr::arrange(group)
  
  M <- as.matrix(W[, -1, drop = FALSE]); storage.mode(M) <- "double"
  S <- stats::cov(M, use = "pairwise.complete.obs")
  times <- as.numeric(colnames(W)[-1])
  list(S = S, times = times)
}

empirical_total_cov <- function(dat, time_col = "time", id_col = "id_full", y_col = "y") {
  df <- dat |>
    dplyr::rename(time = !!time_col, id = !!id_col, y = !!y_col)
  
  W <- df |>
    dplyr::select(id, time, y) |>
    tidyr::pivot_wider(names_from = time, values_from = y) |>
    dplyr::arrange(id)
  
  M <- as.matrix(W[, -1, drop = FALSE]); storage.mode(M) <- "double"
  S <- stats::cov(M, use = "pairwise.complete.obs")
  times <- as.numeric(colnames(W)[-1])
  list(S = S, times = times)
}

ou_start <- function(S, times) {
  r1 <- lag1_rho(S, times)
  dt <- stats::median(diff(sort(times)))
  theta0 <- max(-log(r1) / max(dt, 1e-8), 1e-4)
  var0   <- mean(diag(S), na.rm = TRUE)
  sigma0 <- sqrt(max(2 * theta0 * var0, 1e-8))
  log(c(sigma0, theta0))  # your OUcov expects exp()
}

## -- Wrapper that matches your obj$object$teamCovs ----------------------------
get_team_mom_starts <- function(obj, group_col = "group", y_col = "y") {
  emp  <- empirical_team_cov(obj$model$dat, obj$model$Time, group_col, y_col)
  S    <- emp$S
  time <- emp$times
  
  starts <- vector("list", length(obj$object$teamCovs))
  for (k in seq_along(obj$object$teamCovs)) {
    obj_k <- obj$object$teamCovs[[k]]
    nm  <- obj_k$get_name()
    
    if (identical(nm, "OU.homeostasis")) {
      par <- ou_start(S, time)
      starts[[k]] <- c(0, par)
      
    } else if (identical(nm, "XCov")) {
      X1 <- obj$object$teams[[1]]$X  # adjust if your slot is different
      if (!is.null(X1) && ncol(X1) == 1 && colnames(X1)[1] == "(Intercept)") {
        # pdLogChol uses log-Cholesky params. start with half the sd from there
        starts[[k]] <- 0.2 * 0.5 * log(mean(diag(S), na.rm = TRUE))
      } else {
        starts[[k]] <- rep(0, obj_k$d)
      }
      
    } else {
      starts[[k]] <- rep(0, obj_k$d)
    }
  }
  starts
}

get_indv_mom_start <- function(obj, group_col = "group", id_col = "id_full", y_col = "y") {
  emp_i <- empirical_indiv_cov(obj$model$dat, obj$model$Time, group_col, id_col, y_col)
  S <- emp_i$S; times <- emp_i$times
  
  starts <- vector("list", length(obj$object$indvCovs))
  for (k in seq_along(obj$object$indvCovs)) {
    obj_k <- obj$object$indvCovs[[k]]
    nm  <- obj_k$get_name()
    
    if (identical(nm, "OU.homeostasis")) {
      par <- ou_start(S, times)
      starts[[k]] <- c(0, par)
      
    } else if (identical(nm, "XCov")) {
      X1 <- obj$object$teams[[1]]$indv[[1]]$X  # adjust if your slot is different
      if (!is.null(X1) && ncol(X1) == 1 && colnames(X1)[1] == "(Intercept)") {
        # pdLogChol uses log-Cholesky params. start with 0.2 the sd from there
        starts[[k]] <- 0.2 * 0.5 * log(mean(diag(S), na.rm = TRUE))
      } else {
        starts[[k]] <- rep(0, obj_k$d)
      }
      
    } else {
      starts[[k]] <- rep(0, obj_k$d)
    }
  }
  starts
}

# Random error (expWeightDiag with intercept-only): ensure floor at 10% of total variance
get_error_mom_start <- function(obj, group_col, id_col, y_col) {
  joint_id_col <- "id_joint"
  g <- rlang::sym(group_col)
  i <- rlang::sym(id_col)
  
  obj$model$dat <- obj$model$dat |>
    dplyr::mutate(
      !!joint_id_col := base::as.character(
        base::interaction(!!g, !!i, drop = TRUE, sep = ":")
      )
    )
  
  emp_tot  <- empirical_total_cov(obj$model$dat, obj$model$Time, joint_id_col, y_col)
  emp_team <- empirical_team_cov(obj$model$dat, obj$model$Time, group_col, y_col)
  emp_indv <- empirical_indiv_cov(obj$model$dat, obj$model$Time, group_col, joint_id_col, y_col)
  
  var_tot  <- diag(emp_tot$S)
  var_team <- diag(emp_team$S)
  var_indv <- diag(emp_indv$S)
  
  # residual (per time) for epsilon
  var_eps  <- var_tot - var_team - var_indv
  # global total variance for the floor
  tot_mean_var <- mean(var_tot, na.rm = TRUE)
  floor_val <- 0.10 * tot_mean_var
  var_eps  <- pmax(var_eps, floor_val, 1e-10)
  
  # expWeightDiag with E = ~1 => diag = exp(beta)
  beta0 <- log(mean(var_eps, na.rm = TRUE))
  beta0  # length-1 parameter vector
}

param_list_mom <- function(obj, group_col, id_col, y_col) {
  joint_id_col <- "id_joint"
  g <- rlang::sym(group_col)
  i <- rlang::sym(id_col)
  
  obj$model$dat <- obj$model$dat |>
    dplyr::mutate(
      !!joint_id_col := base::as.character(
        base::interaction(!!g, !!i, drop = TRUE, sep = ":")
      )
    )
  
  list(
    error = get_error_mom_start(obj, group_col = group_col, id_col = joint_id_col, y_col = y_col),
    indv  = get_indv_mom_start(obj, group_col = group_col, id_col = joint_id_col, y_col = y_col),
    team  = get_team_mom_starts(obj, group_col = group_col, y_col = y_col)
  )
}
