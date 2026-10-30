# =============================================================================
# 
# Simulated-data results
# =============================================================================

library(truncnorm)
library(ggplot2)
library(reshape2)
library(MASS)
library(tidyverse)
library(lubridate)
library(patchwork)

source("src/EpiSSM_SEIHR.R")
source("src/ParticleFilter.R")
source("src/PMMH.R")
source("src/PosteriorMarginal.R")

scenarios <- c("sinusoid", "sinusoid_intervention", "piecewise")
panel_ids <- c("A", "B", "C")
scen_labels <- c("A: Sinusoid", "B: Sinusoid with intervention", "C: Piecewise")
traj_files <- paste0("data/simulated/traj_", scenarios, ".csv")      

burnin <- 1000       # iterations (per chain) discarded before sampling theta
Nsub   <- 100        # number of theta draws propagated through the model

stochastic <- TRUE   # must match what was used when the saved fits were produced

oza_name <- "Oza et al."
col_pf_beta   <- "#CD6600"
col_pf_report <- "blue3"
col_oza       <- "green4"
col_obs       <- "red"
col_true      <- "red"

beta_ylim   <- c(0, 0.6)
report_ylim <- NULL          # e.g. c(0, 160); NULL = automatic


col_report    <- "reports"      # simulated reported cases (data)
col_true_beta <- "beta_t_true"  # true beta_t
col_oza_beta  <- "beta"         # Oza et al. beta_t
col_oza_cases <- "C"            # Oza et al. fitted cases
col_day_cands <- c("day", "t", "time")   # day index; if absent, row number is used

# Fixed epidemiological parameters (same as in the paper)
model_param <- list(
  Npop = 5e6, epsilon = 0.5, delta = 1, h = 0.06, gamma0 = 0.25,
  gamma1 = 0.1, rf = 0.425, n_sub = 10
)

make_opts <- function(sim_data) {
  list(
    Data = sim_data, T = nrow(sim_data), obs_type = "raw",
    N = 1000, L = 40, state_dim = 8, obs_dim = 1,
    forecastingHorizon = 0, PredictedObservation = FALSE,
    model_param = model_param,
    stochastic = stochastic,
    pI0    = list(function(N) runif(N, min = 0, max = 5)),
    pBeta0 = list(function(N) rlnorm(N, meanlog = log(0.35), sdlog = 0.2)),
    paramNames = c(expression(epsilon), expression(h),
                   expression(gamma[1]), bquote(sigma[beta])),
    paramPriors  = list(
      list(dist = "lnorm", meanlog = log(0.6), sdlog = 0.2),  # epsilon
      list(dist = "lnorm", meanlog = log(0.07), sdlog = 0.5), # h
      list(dist = "lnorm", meanlog = log(0.12), sdlog = 0.5), # gamma_1
      list(dist = "unif", min = 0, max = 0.5)                 # sigma_beta                
    ),
    lower_bounds = c(1/4, 0.001, 0.01, 0),
    upper_bounds = c(1, 0.5, 1/2, 0.5),
    iterations = 15000, adapt = TRUE, nChains = 1, message = TRUE,
    Nsub = Nsub
  )
}


norm_nm <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

get_col <- function(df, cands, what) {
  hit <- match(norm_nm(cands), norm_nm(names(df)))
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0)
    stop("Column for '", what, "' not found. Available columns: ",
         paste(names(df), collapse = ", "))
  df[[hit[1]]]
}

probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)     # 95% and 50% intervals
quant_df <- function(mat, days) {
  q <- as.data.frame(t(apply(mat, 2, quantile, probs = probs, na.rm = TRUE)))
  names(q) <- c("lo", "q25", "q50", "q75", "hi")
  q$day <- days[seq_len(nrow(q))]
  q
}

rmse <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs)
  sqrt(mean((pred[ok] - obs[ok])^2))
}
mae <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs)
  mean(abs(pred[ok] - obs[ok]))
}

legend_theme <- function(pos = c(0.99, 0.99), just = c(1, 1)) {
  theme(
    legend.position      = pos,
    legend.justification = just,
    legend.background    = element_rect(fill = "white", colour = "grey30", linewidth = 0.5),
    legend.key           = element_rect(fill = "white", colour = NA),
    legend.margin        = margin(1, 2, 1, 2),
    legend.spacing       = unit(1, "pt"),
    legend.key.height    = unit(0.4, "cm"),
    legend.text          = element_text(face = "bold", size = 8)
  )
}

base_theme <- theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 11),
        axis.text  = element_text(face = "bold", size = 10),
        plot.tag   = element_text(face = "bold", size = 18))

lab_true <- "True"
lab_pf   <- "PF (median)"
lab_oza  <- oza_name
lab_obs  <- "Observed data"

one_scenario <- function(i, show_x) {
  scen <- scenarios[i]
  
  # ---- data and Oza estimates, read from the csv
  df <- read.csv(traj_files[i], check.names = FALSE)
  nr <- nrow(df)
  day_hit <- match(norm_nm(col_day_cands), norm_nm(names(df)))
  day_hit <- day_hit[!is.na(day_hit)]
  days <- if (length(day_hit)) df[[day_hit[1]]] else seq_len(nr)
  
  y_obs     <- get_col(df, col_report,    "reported cases")
  beta_true <- get_col(df, col_true_beta, "true beta")
  beta_oza  <- get_col(df, col_oza_beta,  "Oza beta")
  cases_oza <- get_col(df, col_oza_cases, "Oza cases (C)")
  rf <- 0.425
  cases_oza <- rf* cases_oza 
  sim_data <- data.frame(day = days, H_obs = y_obs)
  opts <- make_opts(sim_data)
  
  # ---- saved PMMH fit (never re-run) -> posterior marginals
  rds <- paste0("pmcmc_sim_", scen, ".rds")
  if (!file.exists(rds)) stop("Cannot find saved fit: ", rds)
  mcmcResults <- readRDS(rds)
  
  # NB: BayesianTools::getSample() uses `start` (not `burnin`) to drop burn-in
  thetaSamples <- getSample(mcmcResults, parametersOnly = TRUE, start = burnin,
                            thin = 1, numSamples = Nsub, coda = FALSE)
  set.seed(123)
  marg <- PosteriorMarginal(thetaSamples, EpiSSM_SEIHR, opts)
  
  qb <- quant_df(marg$X[, , 7],        days)    # beta_t
  qr <- quant_df(marg$X_obs[, , 1],    days)    # reported cases
  
  true_df <- data.frame(day = days, y = beta_true)
  oza_b   <- data.frame(day = days, y = beta_oza)
  oza_c   <- data.frame(day = days, y = cases_oza)
  obs_df  <- data.frame(day = days, y = y_obs)
  
  # ---- MAE / RMSE: beta vs TRUE beta; reported cases vs observed data
  #      (PF = posterior median, Oza = their point estimate)
  metrics <- data.frame(
    Scenario = scen_labels[i],
    Method   = c("PF (median)", oza_name),
    MAE_beta    = c(mae(qb$q50,  beta_true), mae(beta_oza,  beta_true)),
    RMSE_beta   = c(rmse(qb$q50, beta_true), rmse(beta_oza, beta_true)),
    MAE_report  = c(mae(qr$q50,  y_obs),     mae(cases_oza, y_obs)),
    RMSE_report = c(rmse(qr$q50, y_obs),     rmse(cases_oza, y_obs))
  )
  
  # ---------------- transmission rate ----------------
  p_beta <- ggplot() +
    geom_ribbon(data = qb, aes(day, ymin = lo,  ymax = hi),  fill = col_pf_beta, alpha = 0.15) +
    geom_ribbon(data = qb, aes(day, ymin = q25, ymax = q75), fill = col_pf_beta, alpha = 0.35) +
    geom_line(data = qb,      aes(day, q50, colour = lab_pf),   linewidth = 0.9) +
     geom_line(data = true_df, aes(day, y,   colour = lab_true), linewidth = 0.9) +
    geom_step(data = oza_b,   aes(day, y,   colour = lab_oza),  linewidth = 0.8) +
   
    scale_colour_manual(
      name   = NULL,
      values = setNames(c(col_true, col_pf_beta, col_oza), c(lab_true, lab_pf, lab_oza)),
      breaks = c(lab_true, lab_pf, lab_oza)) +
    coord_cartesian(xlim = range(days), ylim = beta_ylim) +
    labs(x = if (show_x) "Days" else NULL,
         y = bquote(bold("Transmission rate ") ~ beta[t]),
         tag = panel_ids[i]) +
    base_theme + legend_theme()
  
  # ---------------- reported cases ----------------
  p_rep <- ggplot() +
    geom_ribbon(data = qr, aes(day, ymin = lo,  ymax = hi),  fill = col_pf_report, alpha = 0.15) +
    geom_ribbon(data = qr, aes(day, ymin = q25, ymax = q75), fill = col_pf_report, alpha = 0.35) +
    geom_line(data = qr,     aes(day, q50, colour = lab_pf),  linewidth = 0.8) +
    geom_line(data = obs_df, aes(day, y,   colour = lab_obs), linewidth = 0.7) +
    geom_line(data = oza_c,  aes(day, y,   colour = lab_oza), linewidth = 0.8) +
    scale_colour_manual(
      name   = NULL,
      values = setNames(c(col_obs, col_pf_report, col_oza), c(lab_obs, lab_pf, lab_oza)),
      breaks = c(lab_obs, lab_pf, lab_oza)) +
    coord_cartesian(xlim = range(days), ylim = report_ylim) +
    labs(x = if (show_x) "Days" else NULL, y = "Inpatient cases") +
    base_theme + legend_theme()
  
  if (!show_x) {
    hide_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p_beta <- p_beta + hide_x
    p_rep  <- p_rep  + hide_x
  }
  list(beta = p_beta, rep = p_rep, metrics = metrics)
}

panels <- lapply(seq_along(scenarios), function(i)
  one_scenario(i, show_x = (i == length(scenarios))))

# =============================================================================
# plot
# =============================================================================
fig_sim <- wrap_plots(
  unlist(lapply(panels, function(p) list(p$beta, p$rep)), recursive = FALSE),
  nrow = 3, ncol = 2, byrow = TRUE
)
print(fig_sim)

dir.create("figures/fig_sim", showWarnings = FALSE)
ggsave("figures/fig_sim/fig_sim_3x2.pdf", fig_sim, width = 12, height = 10)

# =============================================================================
# Table: MAE and RMSE for beta_t (vs true beta) and for reported cases
#    (vs observed data), PF vs Oza et al., all scenarios
# =============================================================================
metrics_tab <- bind_rows(lapply(panels, `[[`, "metrics")) %>%
  mutate(MAE_beta = round(MAE_beta, 3), RMSE_beta = round(RMSE_beta, 3),
         MAE_report = round(MAE_report, 2), RMSE_report = round(RMSE_report, 2)) %>%
  rename(`MAE (beta)` = MAE_beta, `RMSE (beta)` = RMSE_beta,
         `MAE (cases)` = MAE_report, `RMSE (cases)` = RMSE_report)

print(as.data.frame(metrics_tab))
