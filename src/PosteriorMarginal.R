# -----------------------------
# Dependencies
# -----------------------------
library(parallel)    # For potential parallelization if desired
library(magrittr)    # For %||% operator if not already defined
source("src/ParticleFilter.R")  # Your BootstrapPF function
# -----------------------------
# Utilities
# -----------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b
# -----------------------------
# PosteriorMarginal
# Always saves X and X_obs, and the per-theta incremental log-likelihood
# matrix (Ntheta_eff x T)
# -----------------------------
PosteriorMarginal <- function(theta, SSM, opts) {
  
  opts$N <- opts$N %||% 4000
  N <- opts$N
  Ns <- opts$Ns %||% 200
  Nsub <- opts$Nsub %||% 100
  opts$PredictedObservation <- TRUE
  
  T <- opts$T
  h <- opts$forecastingHorizon %||% 0
  state_dim <- opts$state_dim
  obs_dim <- opts$obs_dim
  
  thetaSamples <- theta[sample(nrow(theta), Nsub), , drop = FALSE]
  Ntheta_eff <- nrow(thetaSamples)
  
  # Pre-sample indices once
  selected_idx <- matrix(
    sample(1:N, Ns * Ntheta_eff, replace = TRUE),
    nrow = Ntheta_eff,
    ncol = Ns
  )
  
  # Pre-allocate as matrices
  X_all      <- matrix(NA_real_, nrow = Ntheta_eff * Ns, ncol = (T + h) * state_dim)
  X_obs_all  <- matrix(NA_real_, nrow = Ntheta_eff * Ns, ncol = (T + h) * obs_dim)
  loglik_all <- matrix(NA_real_, nrow = Ntheta_eff, ncol = T)   # incremental loglik, t = 1..T
  
  for (i in 1:Ntheta_eff) {
    theta_i <- thetaSamples[i, , drop = TRUE]
    bf_res <- BootstrapPF(SSM, theta_i, opts)
    
    idx <- ((i - 1) * Ns + 1):(i * Ns)
    s_idx <- selected_idx[i, ]
    
    X_all[idx, ]     <- matrix(bf_res$X[s_idx, , ], nrow = Ns)
    X_obs_all[idx, ] <- matrix(bf_res$X_obs[s_idx, , ], nrow = Ns)
    loglik_all[i, ]  <- bf_res$incloglik
  }
  
  # Reshape to arrays
  X_all     <- array(X_all,     dim = c(Ntheta_eff * Ns, T + h, state_dim))
  X_obs_all <- array(X_obs_all, dim = c(Ntheta_eff * Ns, T + h, obs_dim))
  
  list(X = X_all, X_obs = X_obs_all, loglik = loglik_all)
}



# -----------------------------------------------------------------------------
# Helper function 
# -----------------------------------------------------------------------------

# Ribbon plot of an Ntheta x T matrix of trajectories, with optional
# overlaid true/observed points. true_points may be shorter than
# state_mat's number of columns (e.g. when state_mat includes a forecast
# period) -- only the first length(true_points) time steps get an overlay.
plot_marginal_ribbon <- function(state_mat, true_points = NULL,
                                 xlab = "Days", ylab = "",
                                 color = "steelblue") {
  Tn <- ncol(state_mat)
  qdf <- do.call(rbind, lapply(seq_len(Tn), function(t) {
    v <- state_mat[, t]
    data.frame(time = t,
               q05 = quantile(v, 0.05, na.rm = TRUE),
               q25 = quantile(v, 0.25, na.rm = TRUE),
               median = median(v, na.rm = TRUE),
               q75 = quantile(v, 0.75, na.rm = TRUE),
               q95 = quantile(v, 0.95, na.rm = TRUE))
  }))
  p <- ggplot(qdf, aes(x = time)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color, alpha = 0.15) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color, alpha = 0.35) +
    geom_line(aes(y = median), color = color, linewidth = 1.2) +
    labs(x = xlab, y = ylab) +
    theme_bw() +
    theme(
      axis.title = element_text(face = "bold", size = 14),
      axis.text  = element_text(face = "bold", size = 12)
    )
  if (!is.null(true_points)) {
    p <- p + geom_line(data = data.frame(time = seq_len(length(true_points)),
                                         value = true_points),
                       aes(x = time, y = value), color = "red", linewidth = 1)
  }
  p
}