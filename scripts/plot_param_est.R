## ------------------------------------------------------------------
## Posterior densities of all calibrated parameters 
## Trace plot 
## ------------------------------------------------------------------

library(BayesianTools)
library(coda)
library(dplyr)
library(tidyr)
library(ggplot2)

## ---- 1. Settings ---------------------------------------------------

data_type <- "sim"   # "real" or "sim"   <-- the only switch you need

if (data_type == "real") {
  groups       <- c("2022/2023", "2023/2024", "2024/2025", "2025/2026")
  files        <- paste0("pmcmc_real_", gsub("/", "-", groups), ".rds")
  legend_title <- "Season"
  true_vals    <- NULL
  out_tag      <- "real"
} else {
  groups       <- c("sinusoid", "sinusoid_intervention", "piecewise")
  files        <- paste0("pmcmc_sim_", groups, ".rds")
  legend_title <- "Scenario"
  # True values used to simulate the data
  true_vals    <- c(epsilon = 0.5, h = 0.06, gamma1 = 0.1)
  out_tag      <- "sim"
}

# Parameter order MUST match opts$paramNames / opts$paramPriors
par_ids    <- c("epsilon", "h", "gamma1", "sigma_beta")
par_labels <- c("epsilon", "h", "gamma[1]", "sigma[beta]")  # parsed in facet strips

# Bounds of the prior / x-axis for each parameter.
# MUST match opts$lower_bounds and opts$upper_bounds (same order as par_ids).
lower_bounds <- c(1/4, 0.001, 0.001, 0)
upper_bounds <- c(1,   0.4,   0.5,    0.5)

bounds_df <- data.frame(
  param = factor(par_labels, levels = par_labels),
  lower = lower_bounds,
  upper = upper_bounds
)

burnin_frac <- 0.2   # fraction of EACH chain discarded as burn-in
thin        <- 5     # keep every `thin`-th iteration after burn-in

# Priors (copied from opts$paramPriors)
priors <- list(
  epsilon    = function(x) dlnorm(x, meanlog = log(0.6),  sdlog = 0.2),
  h          = function(x) dlnorm(x, meanlog = log(0.07), sdlog = 0.5),
  gamma1     = function(x) dlnorm(x, meanlog = log(0.12), sdlog = 0.5),
  sigma_beta = function(x) dunif(x, min = 0, max = 0.5)
)

## ---- 2. Load results (3 chains each) --------------------------------
## Returns, for one file: posterior samples (burn-in + thinning applied
## per chain), the full log-posterior trace of every chain, the parameter
## traces, and the convergence diagnostics (max R-hat, min ESS).

process_file <- function(file, group) {
  res <- readRDS(file)
  
  # coda = TRUE keeps chains separate (mcmc.list); parametersOnly = FALSE
  # also returns the log-posterior column
  smp <- getSample(res, parametersOnly = FALSE, coda = TRUE, start = 1, thin = 1)
  if (!inherits(smp, "mcmc.list")) smp <- list(smp)
  smp <- lapply(smp, as.matrix)
  
  npar <- length(par_ids)
  lp_col <- function(ch) {
    nm <- colnames(ch)
    if (!is.null(nm) && "Lposterior" %in% nm) "Lposterior" else npar + 1
  }
  
  # ---- posterior draws: burn-in + thin, chain by chain
  par_chains <- lapply(smp, function(ch) {
    n   <- nrow(ch)
    idx <- seq(from = max(1, floor(burnin_frac * n) + 1), to = n, by = thin)
    m   <- ch[idx, seq_len(npar), drop = FALSE]
    colnames(m) <- par_ids
    m
  })
  post <- data.frame(group = group, do.call(rbind, par_chains))
  
  # ---- log-posterior trace (all iterations, to show burn-in too)
  trace <- bind_rows(lapply(seq_along(smp), function(k) {
    ch <- smp[[k]]
    data.frame(group = group, chain = factor(k), iter = seq_len(nrow(ch)),
               logpost = ch[, lp_col(ch)],
               burn = floor(burnin_frac * nrow(ch)))
  }))
  
  # ---- parameter traces (all iterations, to show burn-in too)
  ptrace <- bind_rows(lapply(seq_along(smp), function(k) {
    m <- smp[[k]][, seq_len(npar), drop = FALSE]
    colnames(m) <- par_ids
    data.frame(group = group, chain = factor(k), iter = seq_len(nrow(m)), m)
  }))
  
  # ---- diagnostics on the retained (burned-in, thinned) draws
  mcl  <- mcmc.list(lapply(par_chains, mcmc))
  rhat <- tryCatch(max(gelman.diag(mcl, autoburnin = FALSE,
                                   multivariate = FALSE)$psrf[, 1]),
                   error = function(e) NA_real_)
  ess  <- min(effectiveSize(mcl))   # summed over chains, per parameter
  
  list(post = post, trace = trace, ptrace = ptrace,
       diag = data.frame(group = group, max_rhat = rhat, min_ess = ess))
}

out <- Map(process_file, files, groups)

post_df   <- bind_rows(lapply(out, `[[`, "post"))
trace_df  <- bind_rows(lapply(out, `[[`, "trace"))
diag_df   <- bind_rows(lapply(out, `[[`, "diag"))
ptrace_df <- bind_rows(lapply(out, `[[`, "ptrace")) %>%
  pivot_longer(all_of(par_ids), names_to = "param", values_to = "value") %>%
  mutate(param = factor(param, levels = par_ids, labels = par_labels))

post_df$group  <- factor(post_df$group,  levels = groups)
trace_df$group <- factor(trace_df$group, levels = groups)
diag_df$group  <- factor(diag_df$group,  levels = groups)

post_long <- post_df %>%
  pivot_longer(all_of(par_ids), names_to = "param", values_to = "value") %>%
  mutate(param = factor(param, levels = par_ids, labels = par_labels))

## ---- 3. Prior curves (evaluated over the parameter bounds) ---------

prior_df <- bind_rows(lapply(seq_along(par_ids), function(i) {
  xs <- seq(lower_bounds[i], upper_bounds[i], length.out = 400)
  data.frame(param   = factor(par_labels[i], levels = par_labels),
             value   = xs,
             density = priors[[par_ids[i]]](xs))
}))

# True values (simulated data only)
truth_df <- NULL
if (!is.null(true_vals)) {
  truth_df <- data.frame(
    param = factor(par_labels[match(names(true_vals), par_ids)], levels = par_labels),
    value = unname(true_vals)
  )
}

## ---- 4a. 2x2 posterior density plot (lines only) --------------------
## x-axis of every panel is fixed to the prior bounds, i.e. exactly the
## range over which the prior curve is drawn.

# Long data frame with the two bound values per parameter; used only to
# force each free-scale facet to span [lower, upper].
xlim_df <- bind_rows(
  transmute(bounds_df, param, value = lower),
  transmute(bounds_df, param, value = upper)
)

p <- ggplot(post_long, aes(x = value, colour = group)) +
  # stat_density + geom = "line" gives line (not square) legend keys
  stat_density(geom = "line", position = "identity",
               linewidth = 0.9, trim = TRUE) +
  geom_line(data = prior_df, aes(x = value, y = density),
            inherit.aes = FALSE, linetype = "dashed",
            colour = "black", linewidth = 0.7) +
  geom_blank(data = xlim_df, aes(x = value), inherit.aes = FALSE)

if (!is.null(truth_df)) {
  p <- p + geom_vline(data = truth_df, aes(xintercept = value),
                      colour = "black", linewidth = 0.8)
}

p <- p +
  scale_x_continuous(expand = c(0, 0)) +   # no padding: axis = bounds exactly
  facet_wrap(~ param, nrow = 2, ncol = 2, scales = "free",
             labeller = label_parsed) +
  labs(x = NULL, y = "Posterior density", colour = legend_title) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.3))) +
  theme(axis.title.y      = element_text(size = 14),
        legend.position   = "bottom",
        legend.title      = element_text(size = 13),
        legend.text       = element_text(size = 12),
        legend.key.width  = unit(1.8, "cm"),
        strip.text        = element_text(size = 15),   # parameter names
        strip.background  = element_rect(fill = "grey95"),
        panel.spacing.x   = unit(0.8, "lines"))        # room for edge tick labels

print(p)
dir.create(sprintf("figures/fig_%s", data_type), showWarnings = FALSE, recursive = TRUE)
ggsave(sprintf("figures/fig_%s/posteriors_%s_2x2.pdf",  data_type, out_tag),
       p, width = 8, height = 6.5)


## ---- 4b. Log-posterior trace plot (one panel per season/scenario) ---
## One line per chain (3 DEzs chains); dashed grey line = end of burn-in.

burn_df <- trace_df %>% distinct(group, burn)

p_lp <- ggplot(trace_df, aes(x = iter, y = logpost, colour = chain)) +
  geom_line(linewidth = 0.35, alpha = 0.85) +
  geom_vline(data = burn_df, aes(xintercept = burn),
             linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  facet_wrap(~ group, nrow = 2, ncol = 2, scales = "free_y") +
  scale_colour_brewer(palette = "Dark2", name = "Chain") +
  labs(x = "Iteration", y = "Log posterior") +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.3, alpha = 1))) +
  theme(axis.title        = element_text(size = 14),
        legend.position   = "bottom",
        legend.title      = element_text(size = 13),
        legend.text       = element_text(size = 12),
        legend.key.width  = unit(1.8, "cm"),
        strip.text        = element_text(size = 13),
        strip.background  = element_rect(fill = "grey95"))

print(p_lp)

ggsave(sprintf("figures/fig_%s/logposterior_%s_2x2.pdf",  data_type, out_tag),
       p_lp, width = 10, height = 6)

## ---- 4c. Parameter trace plots: one 2x2 figure per season/scenario --
## One panel per parameter, one line per chain; dashed grey line = end of
## burn-in; solid black horizontal line = true value (simulated data).

make_trace_plot <- function(g) {
  d  <- ptrace_df %>% filter(group == g)
  bn <- floor(burnin_frac * max(d$iter))
  
  pt <- ggplot(d, aes(x = iter, y = value, colour = chain)) +
    geom_line(linewidth = 0.3, alpha = 0.85) +
    geom_vline(xintercept = bn, linetype = "dashed",
               colour = "grey30", linewidth = 0.5)
  
  if (!is.null(truth_df)) {
    pt <- pt + geom_hline(data = truth_df, aes(yintercept = value),
                          colour = "black", linewidth = 0.7)
  }
  
  pt +
    facet_wrap(~ param, nrow = 2, ncol = 2, scales = "free_y",
               labeller = label_parsed) +
    scale_colour_brewer(palette = "Dark2", name = "Chain") +
    labs(x = "Iteration", y = NULL, title = paste(legend_title, "-", g)) +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(linewidth = 1.3, alpha = 1))) +
    theme(axis.title        = element_text(size = 14),
          plot.title        = element_text(hjust = 0.5, face = "bold", size = 15),
          legend.position   = "bottom",
          legend.title      = element_text(size = 13),
          legend.text       = element_text(size = 12),
          legend.key.width  = unit(1.8, "cm"),
          strip.text        = element_text(size = 15),
          strip.background  = element_rect(fill = "grey95"))
}

for (g in groups) {
  pt <- make_trace_plot(g)
  print(pt)
  ggsave(sprintf("figures/fig_%s/trace_params_%s_%s.pdf", data_type, out_tag,
                 gsub("/", "-", g)),
         pt, width = 10, height = 6)
}

## ---- 5. Summary table: median (95% CrI) per group and parameter,
##         plus max R-hat and min ESS ----------------------------------

summary_tab <- post_long %>%
  group_by(group, param) %>%
  summarise(med = median(value),
            lo  = quantile(value, 0.025),
            hi  = quantile(value, 0.975), .groups = "drop") %>%
  mutate(txt = sprintf("%.3f (%.3f, %.3f)", med, lo, hi)) %>%
  select(group, param, txt) %>%
  pivot_wider(names_from = param, values_from = txt) %>%
  left_join(diag_df, by = "group") %>%
  mutate(`Max R-hat` = sprintf("%.3f", max_rhat),
         `Min ESS`   = round(min_ess)) %>%
  select(-max_rhat, -min_ess) %>%
  rename(!!legend_title := group)

print(summary_tab)
write.csv(summary_tab, sprintf("posterior_summary_%s.csv", out_tag), row.names = FALSE)