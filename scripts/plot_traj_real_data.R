# =============================================================================
# Uses the SAVED PMMH fits (pmcmc_real_<season>.rds)
#   Fig 1  effective reproduction number
#   Fig 2  transmission rate beta_t                  
#   Fig 3  reported (inpatient) cases               
#   Fig 4  P(R_eff(t) > 1 | y_{1:T}), all seasons in a single graph
#   Table 1  MAE and RMSE of the reported cases (PF vs Oza)
# =============================================================================

library(truncnorm)
library(ggplot2)
library(patchwork)
library(reshape2)
library(MASS)
library(tidyverse)
library(lubridate)
library(readxl)
library(zoo)

source("src/EpiSSM_SEIHR.R")
source("src/ParticleFilter.R")
source("src/PMMH.R")
source("src/PosteriorMarginal.R")
source("src/impute_missing_reports.R")
# =============================================================================
# 1. USER CHOICES
# =============================================================================
seasons <- c("2022/2023", "2023/2024", "2024/2025", "2025/2026")

burnin  <- 1000          # iterations (per chain) discarded before sampling theta
Nsub    <- 100           # number of theta draws propagated through the model

# Must match what was used when the saved fits were produced
stochastic <- TRUE       # TRUE = chain-binomial model, FALSE = deterministic ODE


oza_name      <- "Oza et al."
cp_file       <- "data/real/traj_est_oza.xlsx"   # one sheet per season (e.g. "2025-2026")
cp_S_col      <- "S"                        # name of the S column in that file
cp_report_col <- "h_obs"                    # Oza reported cases column

col_pf_rt     <- "#CD6600"      # R_eff and beta 
col_pf_report <- "blue3"        # PF reported cases
col_oza       <- "green4"
col_obs       <- "red"

rt_ylim     <- c(0, 3)
beta_ylim   <- c(0, 0.6)
report_ylim <- c(0, 160)

# Fixed epidemiological parameters (same as in the paper)
model_param <- list(
  Npop = 5e6, epsilon = 0.5, delta = 1, h = 0.06, gamma0 = 0.25,
  gamma1 = 0.1, rf = 0.425, n_sub = 10
)

flu_seasons <- read_csv("data/real/flu_seasons.csv", show_col_types = FALSE)

make_opts <- function(season_data) {
  list(
    Data = season_data, T = nrow(season_data), obs_type = "ma7",
    N = 1000, L = 40, state_dim = 8, obs_dim = 1,
    forecastingHorizon = 0, PredictedObservation = FALSE,
    model_param = model_param,
    stochastic = stochastic,
    pI0    = list(function(N) runif(N, min = 0, max = 100)),
    pBeta0 = list(function(N) rlnorm(N, meanlog = log(0.35), sdlog = 0.2)),
    paramNames = c(expression(epsilon), expression(h),
                   expression(gamma[0]), bquote(sigma[beta])),
    paramPriors  = list(
      list(dist = "lnorm", meanlog = log(0.6), sdlog = 0.2),  # epsilon
      list(dist = "lnorm", meanlog = log(0.07), sdlog = 0.5), # h
      list(dist = "lnorm", meanlog = log(0.12), sdlog = 0.5), # gamma_1
      list(dist = "unif", min = 0, max = 0.5)                 # sigma_beta                
    ),
    lower_bounds = c(1/4, 0.01, 0.01, 0),
    upper_bounds = c(1, 0.3, 1/2, 0.5),
    iterations = 15000, adapt = TRUE, nChains = 1, message = TRUE,
    Nsub = Nsub
  )
}

rt_fun <- function(S_frac, beta, h, gamma0, delta) {
  S_frac * beta * (1 / delta + (1 - h) / gamma0)
}

rmse <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs)
  sqrt(mean((pred[ok] - obs[ok])^2))
}
mae <- function(pred, obs) {         
  ok <- is.finite(pred) & is.finite(obs)
  mean(abs(pred[ok] - obs[ok]))
}

probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)   
quant_df <- function(mat, season) {
  q <- as.data.frame(t(apply(mat, 2, quantile, probs = probs, na.rm = TRUE)))
  names(q) <- c("lo", "q25", "q50", "q75", "hi")
  q$day <- seq_len(nrow(q)); q$season <- season
  q
}


# "median (2.5%, 97.5%)" string
ci_str <- function(x, d) {
  q <- unname(quantile(x, c(0.5, 0.025, 0.975), na.rm = TRUE))
  sprintf(paste0("%.", d, "f (%.", d, "f, %.", d, "f)"), q[1], q[2], q[3])
}

# =============================================================================
# 2. Per season: load saved PMMH, propagate, build all quantities
# =============================================================================
one_season <- function(season) {
  tag <- gsub("/", "-", season)
  
  sd <- flu_seasons %>% filter(season == !!season)
  if (nrow(sd) == 0) stop("No data for season ", season)
  season_data <- data.frame(day = sd$day, H_obs = sd$reports) 

  season_data <- impute_small_counts(season_data, lambda = 1.26, seed = 123)
  opts <- make_opts(season_data)
  
  rds <- paste0("pmcmc_real_", tag, ".rds")
  if (!file.exists(rds)) stop("Cannot find saved fit: ", rds)
  mcmcResults <- readRDS(rds)
  
  theta_all <- as.matrix(getSample(mcmcResults, parametersOnly = TRUE,
                                   start = burnin, thin = 1, coda = FALSE))
  thetaSamples <- getSample(mcmcResults, parametersOnly = TRUE, start = burnin,
                            thin = 1, numSamples = Nsub, coda = FALSE)
  
  h_med  <- median(theta_all[, 2])
  
  set.seed(123)
  marg <- PosteriorMarginal(thetaSamples, EpiSSM_SEIHR, opts)
  
  S_mat      <- marg$X[, , 1]
  beta_mat   <- marg$X[, , 7]
  report_mat <- marg$X_obs[, , 1]
  rt_mat     <- rt_fun(S_mat / model_param$Npop, beta_mat, h_med, model_param$gamma0,
                       model_param$delta)
  
  pf_rt   <- quant_df(rt_mat,     season)
  pf_beta <- quant_df(beta_mat,   season)
  pf_obs  <- quant_df(report_mat, season)
  days    <- pf_rt$day
  
  # P(R_eff(t) > 1)
  pf_prob <- data.frame(season = season, day = days,
                        prob = colMeans(rt_mat > 1, na.rm = TRUE))
  
  # ---- Oza et al.
  cp_df <- read_excel(cp_file, sheet = tag)
  if (!cp_S_col %in% names(cp_df))
    stop("Column '", cp_S_col, "' not found in sheet ", tag, " of ", cp_file,
         ". Available columns: ", paste(names(cp_df), collapse = ", "))
  m <- match(days, cp_df$t)
  S_oza    <- cp_df[[cp_S_col]][m]
  beta_oza <- cp_df$beta[m]
  rep_oza  <- cp_df[[cp_report_col]][m]
  if (max(S_oza, na.rm = TRUE) > 1.5) S_oza <- S_oza / model_param$Npop
  
  rt_oza_v <- rt_fun(S_oza, beta_oza, model_param$h,
                     model_param$gamma0, model_param$delta)
  oza_rt   <- data.frame(season = season, day = days, y = rt_oza_v)
  oza_beta <- data.frame(season = season, day = days, y = beta_oza)
  oza_obs  <- data.frame(season = season, day = days, y = rep_oza)
  
  obs <- data.frame(season = season, day = season_data$day, y = season_data$H_obs)
  
  # ---- fit metrics on the median reported cases
  pf_med  <- pf_obs$q50[match(obs$day, pf_obs$day)]
  oza_med <- oza_obs$y[match(obs$day, oza_obs$day)]
  rm_tab <- data.frame(
    season   = season,
    mae_pf   = mae(pf_med,  obs$y), rmse_pf  = rmse(pf_med,  obs$y),
    mae_oza  = mae(oza_med, obs$y), rmse_oza = rmse(oza_med, obs$y)
  )
  

  
  list(pf_rt = pf_rt, pf_beta = pf_beta, pf_obs = pf_obs, pf_prob = pf_prob,
       oza_rt = oza_rt, oza_beta = oza_beta, oza_obs = oza_obs,
       obs = obs, rmse = rm_tab, h_med = h_med)
}

res <- lapply(seasons, one_season)
names(res) <- seasons

collect <- function(nm) {
  d <- bind_rows(lapply(res, `[[`, nm))
  d$season <- factor(d$season, levels = seasons)
  d
}
pf_rt <- collect("pf_rt");   pf_beta <- collect("pf_beta");   pf_obs <- collect("pf_obs")
oza_rt <- collect("oza_rt"); oza_beta <- collect("oza_beta"); oza_obs <- collect("oza_obs")
obs_df  <- collect("obs");   prob_df  <- collect("pf_prob")
rmse_tab <- bind_rows(lapply(res, `[[`, "rmse"))

print(data.frame(season = seasons,
                 h_median = sapply(res, `[[`, "h_med")))
print(rmse_tab)

# =============================================================================
# 3. Figures 1 and 2: R_eff(t) and beta_t 
# =============================================================================
lab_pf  <- "PF (median)"
lab_oza <- oza_name

facet_theme <- theme_bw() +
  theme(axis.title        = element_text(size = 14),
        axis.text         = element_text(size = 11),
        legend.position   = c(0.87, 0.90),
        legend.text       = element_text(size = 10),
        legend.key.width  = unit(1.8, "cm"),
        strip.text        = element_text(size = 12, face = "bold"),
        strip.background  = element_rect(fill = "grey95"))

make_facet_plot <- function(pf, oza, ylab, ylim, hline = NULL) {
  p <- ggplot() +
    geom_ribbon(data = pf, aes(day, ymin = lo,  ymax = hi),  fill = col_pf_rt, alpha = 0.15) +
    geom_ribbon(data = pf, aes(day, ymin = q25, ymax = q75), fill = col_pf_rt, alpha = 0.35)
  if (!is.null(hline))
    p <- p + geom_hline(yintercept = hline, linetype = "dashed",
                        colour = "black", linewidth = 0.6)
  p +
    geom_line(data = pf,  aes(day, q50, colour = lab_pf),  linewidth = 1) +
    geom_step(data = oza, aes(day, y,   colour = lab_oza), linewidth = 0.9) +
    scale_colour_manual(name = NULL,
                        values = setNames(c(col_pf_rt, col_oza), c(lab_pf, lab_oza)),
                        breaks = c(lab_pf, lab_oza)) +
    facet_wrap(~ season, nrow = 2, ncol = 2) +
    labs(x = "Days", y = ylab) +
    coord_cartesian(ylim = ylim) +
    guides(colour = guide_legend(override.aes = list(linewidth = 1.3))) +
    facet_theme
}

dir.create("figures/fig_real", showWarnings = FALSE)

p_rt <- make_facet_plot(pf_rt, oza_rt,
                        bquote(bold("Effective reproduction number ") ~ R[eff](t)),
                        rt_ylim, hline = 1)
print(p_rt)
ggsave("figures/fig_real/fig_Rt_real.pdf", p_rt, width = 10, height = 6)

p_beta <- make_facet_plot(pf_beta, oza_beta,
                          bquote(bold("Transmission rate ") ~ beta[t]),
                          beta_ylim)
print(p_beta)
ggsave("figures/fig_real/fig_beta_real.pdf", p_beta, width = 10, height = 6)

# =============================================================================
# 4. Figure 3: reported cases
# =============================================================================
lab_obs     <- "Observed data"
lab_pf_rep  <- "PF (median)"
lab_oza_rep <- oza_name

p_obs <- ggplot() +
  geom_ribbon(data = pf_obs, aes(day, ymin = lo,  ymax = hi),  fill = col_pf_report, alpha = 0.15) +
  geom_ribbon(data = pf_obs, aes(day, ymin = q25, ymax = q75), fill = col_pf_report, alpha = 0.35) +
  geom_line(data = pf_obs,  aes(day, q50, colour = lab_pf_rep),  linewidth = 0.8) +
  geom_line(data = obs_df,  aes(day, y,   colour = lab_obs),     linewidth = 0.7) +
  geom_line(data = oza_obs, aes(day, y,   colour = lab_oza_rep), linewidth = 0.8) +
  scale_colour_manual(
    name   = NULL,
    values = setNames(c(col_obs, col_pf_report, col_oza),
                      c(lab_obs, lab_pf_rep, lab_oza_rep)),
    breaks = c(lab_obs, lab_pf_rep, lab_oza_rep)
  ) +
  facet_wrap(~ season, nrow = 2, ncol = 2) +
  labs(x = "Days", y = "Inpatient cases") +
  coord_cartesian(ylim = report_ylim) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.3))) +
  facet_theme

print(p_obs)
ggsave("figures/fig_real/fig_obs_real.pdf", p_obs, width = 10, height = 6)

# =============================================================================
# 5. Figure 4: P(R_eff(t) > 1) by season
# =============================================================================
col_traj   <- "#CD950C"   
col_peak   <- "#8A2BE2"   
col_lo     <- "green4"    
col_hi     <- "red"      


peak_obs <- obs_df %>%
  dplyr::group_by(season) %>%
  dplyr::slice_max(y, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(season, peak_day = day)

peak_pts <- prob_df %>%
  dplyr::inner_join(peak_obs, by = "season") %>%
  dplyr::filter(day == peak_day)

print(peak_obs)

y_breaks <- c(0, 0.05, 0.5, 0.95, 1)

p_prob <- ggplot(prob_df, aes(day, prob)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = col_lo, linewidth = 0.6) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = col_hi, linewidth = 0.6) +
  geom_hline(yintercept = 0.5,  linetype = "dashed", colour = "black", linewidth = 0.6) +
  geom_vline(data = peak_obs, aes(xintercept = peak_day),
             colour = col_peak, linetype = "dashed", linewidth = 0.6) +
  geom_line(colour = col_traj, linewidth = 1.1) +
  geom_text(data = peak_obs, aes(x = peak_day, y = 0.12,
                                 label = paste0("peak: day ", peak_day)),
            colour = col_peak, size = 3.3, hjust = 1.08, vjust = 0) +
  scale_y_continuous(limits = c(0, 1), breaks = y_breaks,
                     labels = c("0", "0.05", "0.5", "0.95", "1")) +
  facet_wrap(~ season, nrow = 2, ncol = 2) +
  labs(x = "Days", y = bquote(P(R[eff](t) > 1))) +
  facet_theme +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 11,
                                   colour = c("black", col_lo, "black", col_hi, "black")))

print(p_prob)
ggsave("figures/fig_real/fig_prob_Reff_gt1.pdf", p_prob, width = 10, height = 6)
# =============================================================================
# 6. Table 1: MAE and RMSE of the median reported cases vs observed data
# =============================================================================
metrics_tab <- bind_rows(
  data.frame(Season = rmse_tab$season, Method = "PF (median)",
             MAE = rmse_tab$mae_pf,  RMSE = rmse_tab$rmse_pf),
  data.frame(Season = rmse_tab$season, Method = oza_name,
             MAE = rmse_tab$mae_oza, RMSE = rmse_tab$rmse_oza)
) %>%
  mutate(Season = factor(Season, levels = seasons)) %>%
  arrange(Season, desc(Method == "PF (median)")) %>%
  mutate(Season = as.character(Season),
         MAE = round(MAE, 2), RMSE = round(RMSE, 2))

print(metrics_tab)
