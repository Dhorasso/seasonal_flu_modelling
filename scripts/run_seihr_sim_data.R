

# =============================================================================
# main_seihr_sim_data.R
#
# HOW TO RUN THIS SCRIPT
#   1. Pick (or write) a beta_t(t) scenario in section 4 -- a few examples
#      are included (sinusoid, sinusoid + intervention, constant,
#      decreasing logistic, piecewise). Set Tmax there too.
#   2. Source this file top to bottom. By default it propagates the model
#      using a wide default theta distribution (fast, no fitting) and
#      produces beta_t / reported-case / predictive-log-likelihood plots,
#      with the TRUE beta_t (used to simulate the data) overlaid on the
#      beta_t panel for comparison.
#   3. To fit theta to the data instead, uncomment the "Estimate theta with
#      PMMH" block below -- this replaces the default theta draws with
#      posterior draws from the fit, and everything downstream (plots +
#      Excel export) will then reflect the fitted model.
#   4. Optionally extend the run past the last simulated day by setting
#      opts$forecastingHorizon (see the "Optional: forecast" block).
# =============================================================================

#===================================
# IMPORTANT NOTE
#===================================

library(truncnorm)
library(ggplot2)
library(reshape2)
library(MASS)
library(tidyverse)
library(lubridate)
library(patchwork)
library(cowplot)
library(writexl)     # for exporting the final trajectories to .xlsx

source("src/EpiSSM_SEIHR.R")
source("src/ParticleFilter.R")
source("src/PMMH.R")
source("src/PosteriorMarginal.R")
# =============================================================================
# 2. USER CHOICES
# =============================================================================


# -----------------------------------------------------------------------------
# 3. Fixed epidemiological parameters (see Oza et al (2026) paper)
# -----------------------------------------------------------------------------
model_param <- list(
  Npop    = 5e6,
  epsilon = 0.5,
  delta   = 1,
  h       = 0.06,
  gamma0  = 0.25,
  gamma1  = 0.1,
  rf      = 0.425,
  n_sub   = 10     # Euler substeps per day for the deterministic ODE
)

# -----------------------------------------------------------------------------
# load simulated data
# -----------------------------------------------------------------------------
scenario  <- "sinusoid" # "sinusoid", "sinusoid_intervention", "piecewise"
sim <- read.csv(paste0("data/simulated/traj_", scenario, ".csv"))

sim_data <- data.frame(day = sim$day, H_obs = sim$reports, beta_t_true = sim$beta_t_true)


# Quick look at the simulated series
plot(sim_data$day, sim_data$H_obs,
     xlab = "day", ylab = "reported cases",
     main = paste("Simulated data -", scenario))
plot(sim_data$day, sim_data$beta_t_true,
     xlab = "day", ylab = "beta_t",
     main = paste("Simulated data -", scenario))




# -----------------------------------------------------------------------------
# 5. Options (these default options work well for influenza)
# -----------------------------------------------------------------------------
opts <- list(
  # --- Data ---
  Data = sim_data,
  T    = nrow(sim_data),
  obs_type = "raw", # "raw" for daily count and "ma7" for 7-day average (default "ma7")
  # --- Particle Filter Settings ---
  N                    = 1000, # Number of particles
  L                    = 40,   # Fixed-lag smoothing window
  state_dim            = 8,    # Number of states (see EpiSSM_SEIHR.R)
  obs_dim              = 1,    # Number of observation sources (here: hosp. data only)
  forecastingHorizon   = 0,    # Do not change  [ forecast past the last observed day (see below how to use) ]
  PredictedObservation = FALSE, # Do not change
  
  # --- Fixed epi parameters (see section 3) ---
  model_param = model_param,
  
  # --- Initial-state priors ---
  # Only I0 is customised here (the simulation truth starts at I0 = 1, so a
  # tight-ish prior around that is reasonable); S0, E0, H00, H10, R0 use the
  # model's defaults. To customise any of the others, add e.g.
  #   pE0 = list(function(N) runif(N, min = 0, max = 20))
  pI0     = list(function(N) runif(N, min = 0, max = 5)),
  pBeta0  = list(function(N) rlnorm(N, meanlog = log(0.35), sdlog = 0.2)),

  # --- PMMH priors / bounds ---
  paramNames   = c (
    expression(epsilon),
    # expression(delta),
    expression(h),
    expression(gamma[1]),
    bquote(sigma[beta])
  ),
  
  paramPriors  = list(
    list(dist = "lnorm", meanlog = log(0.6), sdlog = 0.2),  # epsilon
    list(dist = "lnorm", meanlog = log(0.07), sdlog = 0.5), # h
    list(dist = "lnorm", meanlog = log(0.12), sdlog = 0.5), # gamma_1
    list(dist = "unif", min = 0, max = 0.5)                 # sigma_beta                
  ),
  lower_bounds = c(1/4,0.001,1/30,0),
  upper_bounds = c(1,0.5,1/2,0.5),
  
  
  # ── MCMC settings ──────────────────────────────────────────────────────
  iterations = 15000,
  adapt      = TRUE,
  nChains    = 1,
  message    = TRUE
)

# =============================================================================
# STEP A: Quick default theta (used unless Step B is run below)
# =============================================================================
# Fast sanity check before spending time on a full PMMH fit: propagate the
# model using theta ~ Uniform(0.06, 0.12) for sigma (the beta_t random-walk
# volatility) and look at the resulting beta_t and reported-case fit.

opts$Nsub <- 100   # number of theta draws to propagate through the model
set.seed(123)

# thetaSamples <- matrix(runif(opts$Nsub, min = 0.08, max = 0.12), ncol = 1)
# colnames(thetaSamples) <- "sigma"

# -----------------------------------------------------------------------------
# STEP B: Estimate theta with PMMH (optional)
# Uncomment the block below to fit theta from the data instead of using the
# default draws above. thetaSamples is overwritten with posterior draws, so
# everything below this point automatically uses the fitted model.
# -----------------------------------------------------------------------------
mcmcResults <- PMMH(EpiSSM_SEIHR, opts)
marginalPlot(mcmcResults)

saveRDS(mcmcResults, file = paste0("pmcmc_sim_",scenario, ".rds"))
# mcmcResults <- readRDS(paste0("pmcmc_sim_",scenario, ".rds"))
summary(mcmcResults)

thetaSamples <- getSample(mcmcResults, parametersOnly = TRUE, thin = 1,
                          burnin = 1000, numSamples = opts$Nsub, coda = FALSE)

# -----------------------------------------------------------------------------
# Optional: forecast beyond the last simulated day
# -----------------------------------------------------------------------------
# Extend the run past the fit period with a forecast horizon:
# opts$forecastingHorizon <- 28
#
# By default (opts$beta_path left NULL) the forecast holds beta_t at its
# last fitted value. To explore a scenario instead, define beta as a
# vector over the forecast window (day T+1 .. day T+forecastingHorizon) and
# assign it to opts$beta_path, e.g.:

# beta_scenario <- rep(0, opts$forecastingHorizon)
# beta_scenario[1:10]  <- 0.13     # T+1  <= t <= T+10
# beta_scenario[11:20] <- 0.20     # T+11 <= t <= T+20
# beta_scenario[21:28] <- 0.10     # T+21 <= t <= T+28
#
# opts$beta_path <- beta_scenario

# =============================================================================
# STEP C: Propagate theta through the model, detect change points, and plot
# =============================================================================
marg <- PosteriorMarginal(thetaSamples, EpiSSM_SEIHR, opts)

# beta is at index 7 (see EpiSSM_SEIHR.R for the index of other states)
beta_mat   <- marg$X[, , 7]
# reported cases are the (only) observation stream, obs_dim = 1
report_mat <- marg$X_obs[, , 1]



# --- Plots
p_beta <- plot_marginal_ribbon(
  beta_mat, true_points = sim_data$beta_t_true,
  ylab  = bquote(bold("Transmission rate ") ~ beta[t]),
  color = "#CD6600"
) 
print(p_beta)

p_report <- plot_marginal_ribbon(
  report_mat, true_points = sim_data$H_obs,
  ylab  = "Inpatient cases",
  color = "blue"
) 
print(p_report)

fig_marginal <- (p_beta + p_report) +
  plot_annotation(title = paste0("Scenario - ", scenario),
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold")))
print(fig_marginal)

# =============================================================================
# STEP D: Save the median trajectories (and the true beta_t) to an Excel file
# =============================================================================
beta_median   <- apply(beta_mat, 2, median)
report_median <- apply(report_mat, 2, median)

# beta_median / report_median may be longer than sim_data (if a forecast
# horizon was used) -- extend the day index and pad the observed/true
# columns with NA for the forecast days rather than assume the lengths match.
n_total <- length(beta_median)
day_obs <- sim_data$day

day_full <- if (n_total > length(day_obs)) {
  c(day_obs, max(day_obs) + seq_len(n_total - length(day_obs)))
} else {
  day_obs
}

reported_full <- rep(NA_real_, n_total)
reported_full[seq_along(sim_data$H_obs)] <- sim_data$H_obs

beta_true_full <- rep(NA_real_, n_total)
beta_true_full[seq_along(sim$latent$beta)] <- sim_data$beta_t_true

results_df <- data.frame(
  day                 = day_full,
  reported_cases_obs  = reported_full,
  beta_t_true         = beta_true_full,
  beta_t_median       = beta_median,
  report_cases_median = report_median
)

