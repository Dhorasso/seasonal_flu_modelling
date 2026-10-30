
# =============================================================================
# main_seihr_real_data.R
#
# HOW TO RUN THIS SCRIPT
#   1. Set the season to fit or your own data -- it just needs a day column and a case-count
#      column (see the note in section 4).
#   2. Source this file top to bottom. By default it propagates the model
#      using a wide default theta distribution (fast, no fitting) and
#      produces beta_t / reported-case / predictive-log-likelihood plots
#      from that.
#   3. To fit theta to the data instead, uncomment the "Estimate theta with
#      PMMH" block below -- this replaces the default theta draws with
#      posterior draws from the fit (can take a while, ~1h40 in my HPSC PC), and everything downstream (plots + Excel export) will then
#      reflect the fitted model.
#   4. Optionally extend the run past the last observed day by setting
#      opts$forecastingHorizon (see the "Optional: forecast" block).
#   5. Sustained drops in predictive fit are flagged automatically (0-5
#      change points, see src/ChangePointDetection.R) and marked with
#      dashed vertical lines on every plot.
#   6. The final figure is saved as a PDF in Fig_hpsc_real_data/, and the
#      median beta_t / reported-case trajectories (plus the detected
#      change points) are saved to an Excel file in Results/.
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
library(zoo)
library(readxl)

source("src/EpiSSM_SEIHR.R")
source("src/ParticleFilter.R")
source("src/PMMH.R")
source("src/PosteriorMarginal.R")
source("src/impute_missing_reports.R")
# =============================================================================
# 2. USER CHOICES
# =============================================================================
season <- "2025/2026"     # <-- change this to fit a different season

# -----------------------------------------------------------------------------
# 3. Fixed epidemiological parameters (same as in the paper)
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
# 4. Load flu season data
#    season_data must end up as a 2-column data frame:
#      col 1 "day"   -- integer time index, 1, 2, 3, ...
#      col 2 "H_obs" -- observed case count for that day (NA on censored days
#                       -- ignored by ObsProcess when censored)
#    If you're using your own file instead of flu_seasons_extracted.csv,
#    just build season_data with those two columns directly and skip the
#    get_season() filtering step below.
# -----------------------------------------------------------------------------

# get_season() helper -- same one used to build the csv, kept here so this
# script is self-contained
get_season <- function(data, season_name) {
  data %>% filter(season == !!season_name)
}

flu_seasons <- read_csv("data/real/flu_seasons.csv", show_col_types = FALSE)

season_data <- get_season(flu_seasons, season)

if (nrow(season_data) == 0) {
  stop("No data found for season '", season, "' in flu_seasons.csv")
}

season_data <- data.frame(
  day   = season_data$day,
  H_obs = season_data$reports
)
## impute missing reports (skip this if you have the complete data)
season_data <- impute_small_counts(season_data, lambda = 1.26, seed = 123) # lamda is the poisson rate
# Quick look at the raw series
plot(season_data$day, season_data$H_obs,
     xlab = "day", ylab = "reported cases",
     main = paste("Season", season))

# -----------------------------------------------------------------------------
# 5. Options (these default options work well for influenza)
# -----------------------------------------------------------------------------
opts <- list(
  # --- Data ---
  Data = season_data,
  T    = nrow(season_data),
  obs_type = "ma7",  # "raw" for daily count and "ma7" for 7-day average (default "ma7")
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
  # Only I0 is customised here; S0, E0, H00, H10, R0 use the model's
  # defaults. To customise any of the others, add e.g.
  #   pE0 = list(function(N) runif(N, min = 0, max = 20))
  # Any distribution can be used here, e.g. runif(), rlnorm(), rgamma(), ...
  pI0     = list(function(N) runif(N, min = 0, max = 100)),
  pBeta0  = list(function(N) rlnorm(N, meanlog = log(0.35), sdlog = 0.2)),
  
  paramNames   = c (
    expression(epsilon),
    # expression(delta),
    expression(h),
    expression(gamma[1]),
    bquote(sigma[beta])
  ),
  # --- PMMH priors / bounds ---
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

# -----------------------------------------------------------------------------
# 6. Helper function (shared by every plot below)
# -----------------------------------------------------------------------------

# Ribbon plot of an Ntheta x T matrix of trajectories, with optional
# overlaid observed points. true_points may be shorter than state_mat's
# number of columns (e.g. when state_mat includes a forecast period) --
# only the first length(true_points) time steps get an overlay line.
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
      axis.title = element_text(face = "bold", size = 10),
      axis.text  = element_text(face = "bold", size = 10)
    )
  if (!is.null(true_points)) {
    p <- p + geom_line(data = data.frame(time = seq_len(length(true_points)),
                                         value = true_points),
                       aes(x = time, y = value), color = "black", linewidth = 1)
  }
  p
}

# =============================================================================
# STEP A: Quick default theta (used unless Step B is run below)
# =============================================================================
# Fast sanity check before spending time on a full PMMH fit: propagate the
# model using theta ~ Uniform(0.08, 0.12) for sigma (the beta_t random-walk
# volatility) and look at the resulting beta_t and reported-case fit.



# thetaSamples <- matrix(runif(opts$Nsub, min = 0.08, max = 0.12), ncol = 1)
# colnames(thetaSamples) <- "sigma"

# -----------------------------------------------------------------------------
# STEP B: Estimate theta with PMMH 
# -----------------------------------------------------------------------------
# 
mcmcResults <- PMMH(EpiSSM_SEIHR, opts)
marginalPlot(mcmcResults)

saveRDS(mcmcResults, file = paste0("pmcmc_real_", gsub("/", "-", season), ".rds"))
# mcmcResults <- readRDS(paste0("pmcmc_real_", gsub("/", "-", season), ".rds"))
summary(mcmcResults)

opts$Nsub <- 100   # number of theta draws to propagate through the model
set.seed(123)

thetaSamples <- getSample(mcmcResults, parametersOnly = TRUE, thin = 1,
                           burnin = 1000, numSamples = opts$Nsub, coda = FALSE)

# -----------------------------------------------------------------------------
# Optional: forecast beyond the last observed day
# -----------------------------------------------------------------------------
# Extend the run past the fit period with a forecast horizon:
#   opts$forecastingHorizon <- 14
#
# By default (opts$beta_path left NULL) the forecast holds beta_t at its
# last fitted value. To explore a scenario instead, define beta as a
# function of the *absolute* day number and evaluate it over the forecast
# window (day T+1 .. day T+forecastingHorizon), e.g.:
#
# beta_scenario <- function(day) {
#   dplyr::case_when(
#     day > 125 & day <= 135 ~ 0.12,
#     day > 135 & day <= 165 ~ 0.19,
#     day > 165              ~ 0.12,
#     TRUE                   ~ NA_real_
#   )
# }
# opts$forecastingHorizon <- 14
# opts$beta_path <- beta_scenario(seq(opts$T + 1, opts$T + opts$forecastingHorizon))


# =============================================================================
# Propagate theta through the model, detect change points, and plot
# =============================================================================
marg <- PosteriorMarginal(thetaSamples, EpiSSM_SEIHR, opts)

beta_mat   <- marg$X[, , 7]        # beta_t
report_mat <- marg$X_obs[, , 1]    # reported cases


# Medians from the PF
beta_median   <- apply(beta_mat, 2, median)
report_median <- apply(report_mat, 2, median)
n_total <- length(beta_median)

# Day index (extended if a forecast horizon was used)
day_obs  <- season_data$day
day_full <- if (n_total > length(day_obs)) {
  c(day_obs, max(day_obs) + seq_len(n_total - length(day_obs)))
} else day_obs



# =============================================================================
# plot
# =============================================================================

col_pf_beta   <- "#CD6600"
col_pf_report <- "blue3"


# ---- beta_t ---------------------------------------------------------
p_beta <- plot_marginal_ribbon(
  beta_mat,
  ylab  = bquote(bold("Transmission rate ") ~ beta[t]),
  color = col_pf_beta
) 

# ---- reported cases -----------------------------------

p_report <- plot_marginal_ribbon(
  report_mat, true_points = NULL,
  ylab  = "Inpatient cases",
  color = col_pf_report
) +
  geom_line(data = season_data, aes(x = day, y = H_obs, colour = "red"),
            linewidth = 0.7, linetype ="solid") 

fig_marginal <- (p_beta + p_report)+
  plot_annotation(title = paste0("Season - ", season),
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16,
                                                          face = "bold")))
print(fig_marginal)


# =============================================================================
# Save trajectories
# =============================================================================

reported_full <- rep(NA_real_, n_total)
reported_full[seq_along(season_data$H_obs)] <- season_data$H_obs

results_df <- data.frame(
  day            = day_full,
  cases_obs      = reported_full,
  beta_t_est     = beta_median,
  cases_est      = report_median
)
