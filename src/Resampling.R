# ===========================
# Resampling.R
# ===========================

# Purpose: Resampling functions for particle filters

# ---------------------------
# Resample
# ---------------------------
# Inputs:
#   w: numeric vector of normalized or unnormalized weights (length N)
#   method: character, one of "systematic", "multinomial", "stratified"
# Output:
#   Integer vector of indices (length N) representing resampled particles

# ===========================
# Resampling.R
# ===========================

Resample <- function(w, method = "systematic") {
  
  valid_methods <- c("multinomial", "systematic", "stratified")
  if (!(method %in% valid_methods)) {
    stop(sprintf("Resampling method '%s' not implemented. Choose from: %s",
                 method, paste(valid_methods, collapse = ", ")))
  }
  
  N <- length(w)
  
  # --- Clean weights ---
  w[!is.finite(w) | w < 0] <- 0
  
  total <- sum(w)
  
  # --- If degenerate → use uniform weights ---
  if (total <= 0 || !is.finite(total)) {
    w <- rep(1 / N, N)
  } else {
    w <- w / total
  }
  
  # Guarantee exact normalization
  w <- w / sum(w)
  
  # ---------------------------
  # MULTINOMIAL
  # ---------------------------
  if (method == "multinomial") {
    return(sample.int(N, size = N, replace = TRUE, prob = w))
  }
  
  # Cumulative sum
  cs <- cumsum(w)
  cs[N] <- 1.0
  
  indexes <- integer(N)
  
  # ---------------------------
  # SYSTEMATIC
  # ---------------------------
  if (method == "systematic") {
    
    u0 <- runif(1) / N
    positions <- u0 + (0:(N-1)) / N
    
    i <- 1
    j <- 1
    
    while (i <= N) {
      if (positions[i] <= cs[j]) {
        indexes[i] <- j
        i <- i + 1
      } else {
        j <- j + 1
        if (j > N) j <- N  # clamp instead of stopping
      }
    }
    
    return(indexes)
  }
  
  # ---------------------------
  # STRATIFIED
  # ---------------------------
  if (method == "stratified") {
    
    positions <- ((0:(N-1)) + runif(N)) / N
    
    i <- 1
    j <- 1
    
    while (i <= N) {
      if (positions[i] <= cs[j]) {
        indexes[i] <- j
        i <- i + 1
      } else {
        j <- j + 1
        if (j > N) j <- N
      }
    }
    
    return(indexes)
  }
}
