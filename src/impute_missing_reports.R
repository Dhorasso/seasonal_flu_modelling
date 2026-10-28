impute_small_counts <- function(season_data, col = "H_obs", lo = 0L, hi = 4L,
                                lambda = 1.26, seed = NULL) {
  
  miss <- is.na(season_data[[col]])
  n_miss <- sum(miss)
  
  if (n_miss == 0) return(season_data)
  
  # Set the seed temporarily, then restore the global RNG state afterwards
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  
  # Truncated Poisson: Poisson probabilities restricted to lo:hi and renormalised
  vals  <- lo:hi
  probs <- dpois(vals, lambda)
  probs <- probs / sum(probs)
  
  season_data[[col]][miss] <- sample(vals, n_miss, replace = TRUE, prob = probs)
  season_data
}