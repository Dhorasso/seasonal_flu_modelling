# ===========================================================================
# EpiSSM_SEIHR.R
# ===========================================================================
# SEIHR state-space model -- (stochastic/deterministic)

#
# State columns (state_dim = 8):
#   1 S   2 E   3 I   4 H0   5 H1   6 R   7 beta   8 H1new
#
#   H1new = new I->H1 entries over the day. This is the observation mean
#           (scaled by rf) -- no delay convolution is needed since the
#           data are already "new hospitalizations".
#
# beta_t process:
#   - For t <= opts$T (the observed/fit period): beta_t follows a latent
#     log-normal random walk, exp(rnorm(log(beta_prev), sigma)).
#   - For t >  opts$T (a forecast period, i.e. opts$forecastingHorizon > 0):
#       * if opts$beta_path is supplied, beta_t is fixed to
#         opts$beta_path[t - opts$T] for every particle (a beta scenario)
#       * otherwise beta_t is held at its last fitted value
#
# Initial-state priors (t = 1):
#   Every compartment except S can be given a sampler in opts, e.g.:
#     opts$pE0  = list(function(N) rep(0, N))
#     opts$pI0  = list(function(N) runif(N, 1, 100))
#     opts$pH00 = list(function(N) rep(0, N))
#     opts$pH10 = list(function(N) rep(0, N))
#     opts$pR0  = list(function(N) rep(0, N))
#   Any of these left unset falls back to the default below:
#     E0 = H00 = H10 = R0 = 0,   I0 = 1
#   S0 is never sampled directly -- it is always whatever is left of the
#   population: S0 = Npop - E0 - I0 - H00 - H10 - R0.
# ===========================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

# ------------------------------------------------------------
# 1. Initial state sampler  (t = 1)
# ------------------------------------------------------------
InitState <- function(N, theta, opts, Y_hist = NULL, ...) {
  
  Npop <- opts$model_param$Npop
  
  # Draw a compartment's initial value from its user-supplied prior if one
  # was given in opts, otherwise use the stated default.
  draw_state <- function(prior_list, default_fun) {
    if (!is.null(prior_list)) prior_list[[1]](N) else default_fun(N)
  }
  
  E0  <- draw_state(opts$pE0,  function(N) rep(0, N))
  I0  <- draw_state(opts$pI0,  function(N) rep(1, N))
  H00 <- draw_state(opts$pH00, function(N) rep(0, N))
  H10 <- draw_state(opts$pH10, function(N) rep(0, N))
  R0  <- draw_state(opts$pR0,  function(N) rep(0, N))
  
  # S always absorbs whatever is left of the population
  S0 <- Npop - E0 - I0 - H00 - H10 - R0
  
  beta0  <- opts$pBeta0[[1]](N)
  H1new0 <- rep(0, N)
  
  cbind(S = S0, E = E0, I = I0, H0 = H00, H1 = H10, R = R0,
        beta = beta0, H1new = H1new0)
}

# ------------------------------------------------------------
# 2. State transition process (t >= 2), deterministic ODE only
#    X_hist: [N x t x state_dim], t = current real time step.
# ------------------------------------------------------------
StateProcess <- function(X_hist, t, theta, opts, ...) {
  
  N <- dim(X_hist)[1]
  
  S_prev  <- X_hist[, t - 1L, 1L]
  E_prev  <- X_hist[, t - 1L, 2L]
  I_prev  <- X_hist[, t - 1L, 3L]
  H0_prev <- X_hist[, t - 1L, 4L]
  H1_prev <- X_hist[, t - 1L, 5L]
  R_prev  <- X_hist[, t - 1L, 6L]
  
  mp       <- opts$model_param
  eps_rate <- theta[1] #mp$epsilon
  delta    <- mp$delta
  h        <- theta[2] #mp$h
  gamma0   <- mp$gamma0
  gamma1   <- theta[3] #mp$gamma1
  sigma <- theta[4]
  # --- beta_t process ---
  # Observed/fit period (t <= opts$T): latent log-normal random walk.
  # Forecast period (t > opts$T): either a fixed scenario (opts$beta_path)
  # or held at the last fitted value if no scenario is supplied.
  stochastic <-opts$stochastic %||% TRUE
  
  if (t <= opts$T) {
    beta_prev <- X_hist[, t - 1L, 7L]
    
    beta_t <- exp(rnorm(N, mean = log(beta_prev), sd = sigma))
  } else if (!is.null(opts$beta_scenario)) {
    beta_mean <- rep(opts$beta_scenario[t - opts$T], N)
    beta_t <- exp(rnorm(N, mean = log(beta_mean), sd = 0.1))
  } else {
    beta_t <- X_hist[, t - 1L, 7L]
  }
  
  n_sub <- mp$n_sub %||% 10
  dt <- 1 / n_sub
  
  S <- S_prev; E <- E_prev; I <- I_prev
  H0 <- H0_prev; H1 <- H1_prev; R <- R_prev
  H1new_acc <- rep(0, N)
  
  for (s in seq_len(n_sub)) {
    
    if (stochastic) {
      # ---- chain-binomial transitions ----
      Si  <- pmax(round(S),  0); Ei  <- pmax(round(E),  0)
      Ii  <- pmax(round(I),  0); H0i <- pmax(round(H0), 0)
      H1i <- pmax(round(H1), 0); Ri  <- pmax(round(R),  0)
      
      Npop_t <- pmax(Si + Ei + Ii + H0i + H1i + Ri, 1)
      lambda <- beta_t * (Ii + H0i) / Npop_t
      
      X1 <- rbinom(N, Si,  1 - exp(-lambda   * dt))
      X2 <- rbinom(N, Ei,  1 - exp(-eps_rate * dt))
      X3 <- rbinom(N, Ii,  1 - exp(-delta    * dt))
      X4 <- rbinom(N, X3,  h)
      X5 <- rbinom(N, H0i, 1 - exp(-gamma0   * dt))
      X6 <- rbinom(N, H1i, 1 - exp(-gamma1   * dt))
      
      S  <- Si - X1
      E  <- Ei + X1 - X2
      I  <- Ii + X2 - X3
      H0 <- H0i + X3 - X4 - X5
      H1 <- H1i + X4 - X6
      R  <- Ri + X5 + X6
      
      H1new_acc <- H1new_acc + X4
      
    } else {
      # ---- deterministic ODE, Euler substep ----
      Npop_t <- pmax(S + E + I + H0 + H1 + R, 1)
      lambda <- beta_t * (I + H0) / Npop_t
      
      dS  <- -lambda * S
      dE  <- lambda * S - eps_rate * E
      dI  <- eps_rate * E - delta * I
      dH0 <- delta * (1 - h) * I - gamma0 * H0
      dH1 <- delta * h * I       - gamma1 * H1
      dR  <- gamma0 * H0 + gamma1 * H1
      
      H1new_acc <- H1new_acc + delta * h * I * dt
      
      S  <- S  + dS  * dt
      E  <- E  + dE  * dt
      I  <- I  + dI  * dt
      H0 <- H0 + dH0 * dt
      H1 <- H1 + dH1 * dt
      R  <- R  + dR  * dt
    }
  }
  
  
  # keep populations non-negative under discretisation noise
  S <- pmax(S, 0); E <- pmax(E, 0); I <- pmax(I, 0)
  H0 <- pmax(H0, 0); H1 <- pmax(H1, 0); R <- pmax(R, 0)
  
  cbind(S = S, E = E, I = I, H0 = H0, H1 = H1, R = R,
        beta = beta_t, H1new = H1new_acc)
}

# ------------------------------------------------------------
# 3. Observation process
# ------------------------------------------------------------
# obs_type:
#   "raw"         : observed raw daily case counts
#   "ma7"         : observed 7-day moving-average case counts
#
# For ma7, use a Gaussian approximation to the distribution of
# the moving average of daily Poisson reporting counts.
#
# Optional:
#   opts$obs_type   = "raw" or "ma7"
#   opts$ma_window  = 7
#   opts$obs_sd     = additional observation-error SD (default 0)
#
# Assumptions:
#   - H1new is the model's daily incidence.
#   - rf is the reporting probability/fraction.
#   - X_hist[, t, 8L] contains daily incidence for each particle.
# ------------------------------------------------------------
ObsProcess <- function(X_hist, Y_hist, t, theta, opts, ...) {
  
  N  <- dim(X_hist)[1]
  rf <- opts$model_param$rf
  
  # Select observation type (default: 7-day moving average)
  obs_type <- opts$obs_type %||% "ma7"
  
  H1new <- X_hist[, t, 8L]
  
  logw  <- rep(0, N)
  C_obs <- rep(NA_real_, N)   # defined even when no weighting happens
  
  # Do not weight particles during a forecast period
  if (t >= 1L && t <= opts$T) {
    Y_t <- Y_hist[t, 2L]
    
    if (obs_type == "raw") {
      # CASE 1: Raw daily case counts
      mu_today <- pmax(rf * H1new, 1e-6)      
      if (!is.na(Y_t)){
        logw     <- dpois(x = Y_t, lambda = mu_today, log = TRUE)
      }
      C_obs    <- rpois(N, lambda = mu_today)
      
    } else if (obs_type == "ma7") {
      # CASE 2: 7-day moving-average case counts
      idx <- max(1L, t - 7L + 1L):t
      mu_ma7 <- pmax(
        rf / 7 * rowSums(X_hist[, idx, 8L, drop = FALSE]),
        1e-6
      )
      
      # Gaussian observation likelihood
      est_sd = (max(2,  mean(mu_ma7)))^0.5
      if (!is.na(Y_t)){
        logw  <- dnorm(x = Y_t, mean = mu_ma7, sd = est_sd , log = TRUE)
      }
      C_obs <- pmax(rnorm(N, mean = mu_ma7, sd = est_sd), 0)
      
    } else {
      stop("Unknown obs_type: ", obs_type)
    }
  }
  
  list(
    logw  = logw,
    X_obs = matrix(C_obs, nrow = N, ncol = 1L)
  )
}

# ============================================================
# Assemble SSM module
# ============================================================
EpiSSM_SEIHR <- list(
  InitState    = InitState,
  StateProcess = StateProcess,
  ObsProcess   = ObsProcess
)