
# ============================================================
# ParticleFilter.R
#
# ==========================================================

library(compiler)

#' Run a full bootstrap particle filter over t = 1 .. T (+ forecast horizon)
#'
#' @param SSM state-space model list (`InitState`, `StateProcess`, `ObsProcess`)
#' @param theta parameter vector
#' @param opts options list. Relevant fields:
#'   \describe{
#'     \item{N}{number of particles}
#'     \item{T}{number of observed time steps}
#'     \item{forecastingHorizon}{extra steps to simulate beyond T (default 0)}
#'     \item{L}{fixed-lag smoothing window (default 0)}
#'     \item{state_dim, obs_dim}{state / observation dimension}
#'     \item{PredictedObservation}{store predicted observations? (default FALSE)}
#'     \item{resample_method}{"systematic" (default), "multinomial", "stratified"}
#'   }
#' @param X_init optional initial particle matrix [N, state_dim]; otherwise
#'   sampled from `SSM$InitState()`
#' @param ... passed through to `SSM$InitState()`
#' @return list with `X` [N, T+h, state_dim], `Loglike` (scalar total),
#'   `loglik` (length-T vector of incremental log-likelihoods, t = 1..T),
#'   and optionally `X_obs`
BootstrapPF <- function(SSM, theta, opts, X_init = NULL, ...) {
  
  Y         <- opts$Data
  N         <- opts$N
  T         <- opts$T
  h         <- opts$forecastingHorizon %||% 0
  L         <- opts$L %||% 0
  state_dim <- opts$state_dim
  obs_dim   <- opts$obs_dim
  
  store_X_obs     <- opts$PredictedObservation %||% FALSE
  resample_method <- opts$resample_method %||% "systematic"
  
  X     <- array(0, dim = c(N, T + h, state_dim))
  X_obs <- if (store_X_obs) array(0, dim = c(N, T + h, obs_dim)) else NULL
  Loglike  <- 0
  loglik_t <- numeric(T)   # incremental log-likelihood, one entry per t = 1..T
  
  X[, 1, ] <- if (is.null(X_init)) SSM$InitState(N, theta, opts, Y, ...) else X_init
  
  for (t in seq_len(T)) {
    if (t > 1L) X[, t, ] <- SSM$StateProcess(X, t, theta, opts)
    
    obs  <- SSM$ObsProcess(X, Y, t, theta, opts)
    logw <- obs$logw
    if (store_X_obs) X_obs[, t, ] <- obs$X_obs
    
    m    <- max(logw)
    ll_t <- m + log(sum(exp(logw - m))) - log(N)
    loglik_t[t] <- ll_t
    Loglike <- Loglike + ll_t
    
    w_vec <- exp(logw - m)
    w_vec <- w_vec / sum(w_vec)
    indx  <- Resample(w_vec, method = resample_method)
    
    tmin <- max(1L, t - L)
    X[, tmin:t, ] <- X[indx, tmin:t, , drop = FALSE]
    if (store_X_obs) X_obs[, tmin:t, ] <- X_obs[indx, tmin:t, , drop = FALSE]
  }
  
  if (h > 0) {
    for (t in (T + 1L):(T + h)) {
      X[, t, ] <- SSM$StateProcess(X, t, theta, opts)
      if (store_X_obs) {
        obs <- SSM$ObsProcess(X, Y, t, theta, opts)
        X_obs[, t, ] <- obs$X_obs
      }
    }
  }
  
  res <- list(X = X, Loglike = Loglike, incloglik = loglik_t)
  if (store_X_obs) res$X_obs <- X_obs
  res
}
BootstrapPF <- cmpfun(BootstrapPF)




# ---------------------------
# Resample
# ---------------------------
# Inputs:
#   w: numeric vector of normalized or unnormalized weights (length N)
#   method: character, one of "systematic", "multinomial", "stratified"
# Output:
#   Integer vector of indices (length N) representing resampled particles
Resample <- function(w, method = "systematic") {
  methods <- c("multinomial", "systematic", "stratified")
  if (!method %in% methods)
    stop("Invalid method. Choose: ", paste(methods, collapse = ", "))
  
  N <- length(w)
  w[!is.finite(w) | w < 0] <- 0
  w <- w / sum(w)
  
  if (!is.finite(sum(w)) || sum(w) == 0)
    w <- rep(1 / N, N)
  
  w <- w / sum(w)
  
  if (method == "multinomial")
    return(sample.int(N, N, replace = TRUE, prob = w))
  
  cs <- cumsum(w)
  cs[N] <- 1
  
  positions <- if (method == "systematic")
    runif(1) / N + (0:(N - 1)) / N
  else
    (0:(N - 1) + runif(N)) / N
  
  indexes <- integer(N)
  i <- j <- 1
  
  while (i <= N) {
    if (positions[i] <= cs[j]) {
      indexes[i] <- j
      i <- i + 1
    } else {
      j <- min(j + 1, N)
    }
  }
  
  indexes
}