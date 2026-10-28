source("src/ParticleFilter.R")


if (!requireNamespace("BayesianTools", quietly = TRUE)) {
  install.packages("BayesianTools")
}
library(BayesianTools)
library(compiler)

PMMH_inner <- function(SSM, opts) {
  
  # ---------------------------
  # Likelihood function
  # ---------------------------
  # Note: No caching because:
  # 1. Bootstrap PF is stochastic (different result each call)
  # 2. MCMC explores continuous space (exact matches are rare)
  # 3. Cache overhead > benefit in practice
  
  likelihood <- function(theta) {
    
    if (!is.null(opts$lower_bounds) && !is.null(opts$upper_bounds)) {
      if (any(theta < opts$lower_bounds) ||
          any(theta > opts$upper_bounds)) {
        return(-Inf)
      }
    }
    
    BootstrapPF(SSM, theta, opts)$Loglike
  }
  
  
  # ---------------------------
  # Prior density
  # ---------------------------
  logPriorDensity <- function(theta) {
    
    if (!is.null(opts$lower_bounds) && !is.null(opts$upper_bounds)) {
      if (any(theta < opts$lower_bounds) ||
          any(theta > opts$upper_bounds)) {
        return(-Inf)
      }
    }
    
    computeLogPrior(theta, opts$paramPriors)
  }
  
  
  # ---------------------------
  # Prior sampler
  # ---------------------------
  priorSampler <- function(n = 1) {
    nParams <- length(opts$paramPriors)
    samples <- matrix(NA, nrow = n, ncol = nParams)
    
    for (i in seq_len(nParams)) {
      priorSpec <- opts$paramPriors[[i]]
      
      if (is.list(priorSpec) && !is.null(priorSpec$dist)) {
        distFun <- get(paste0("r", priorSpec$dist))
        args <- priorSpec
        args$dist <- NULL
        args$n <- n
        samples[, i] <- do.call(distFun, args)
      } else if (is.function(priorSpec)) {
        samples[, i] <- priorSpec(n)
      } else {
        stop("Invalid prior specification")
      }
    }
    
    samples
  }
  
  # ---------------------------
  # Prior object
  # ---------------------------
  prior <- createPrior(
    density = logPriorDensity,
    sampler = priorSampler,
    lower   = opts$lower_bounds %||% NULL ,
    upper   = opts$upper_bounds  %||% NULL
  )
  
  # ---------------------------
  # Bayesian setup
  # ---------------------------
  bayesianSetup <- createBayesianSetup(
    likelihood = likelihood,
    prior = prior,
    names = opts$paramNames
    # parallel = TRUE
  )
  iterations <- opts$iterations
  settings <- list(
    iterations = iterations%||% 10000,
    nrChains   = opts$nChains %||% 2,
    # burnin     = opts$burnin %||% floor(0.3* iterations),
    # thin   = 3,
    message    = opts$message %||% TRUE
  )
  
  mcmcResults <- runMCMC(
    bayesianSetup = bayesianSetup,
    sampler = "DEzs",
    settings = settings
  )
}

PMMH <- cmpfun(PMMH_inner)




# --------------------------
# Compute log-prior including Jacobian and optional parameter limits
# --------------------------
# theta: parameter vector (original space)
# priors: list of prior specifications
#         - either a function returning log-density
#         - or a list with dist = "norm"/"gamma"/"unif"/"beta"/"lnorm" and parameters
# transforms: list per parameter: "none", "log", "logit", "bounded_logit"
# ParamLimits: optional matrix/data.frame with columns c(lower, upper) per parameter


computeLogPrior <- function(theta, priors) {
  nParams <- length(theta)
  logPriorSum <- 0
  
  for (i in seq_len(nParams)) {
    x <- theta[i]
    
    priorSpec <- priors[[i]]
    logp <- NA
    
    if (is.list(priorSpec) && !is.null(priorSpec$dist)) {
      distFun <- get(paste0("d", priorSpec$dist))
      args <- priorSpec
      args$x <- x
      args$log <- TRUE
      args$dist <- NULL
      logp <- do.call(distFun, args)
      
      if(is.na(logp)) logp <- -Inf
      
    } else if (is.function(priorSpec)) {
      logp <- priorSpec(x)
      if(is.na(logp)) logp <- -Inf
    } else {
      stop("Invalid prior specification")
    }
    
    logPriorSum <- logPriorSum + logp 
  }
  
  return(logPriorSum)
}
