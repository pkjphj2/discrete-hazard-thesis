# Thesis Section 3.3: simulation study (Figures 3.4 bias, 3.5 MSE).
# Run from the repository root: Rscript analysis/03_simulation.R   (about 5 minutes)
# Data generation and censoring are kept exactly as in the thesis code so set.seed(202507)
# reproduces the published figures draw for draw.
suppressMessages({library(dplyr); library(tidyr); library(ggplot2)})
dir.create("results", showWarnings = FALSE)
set.seed(202507)

## 1) Baseline S0 and derived quantities
S0_geom <- function(k, p)        (1 - p)^(k + 1)
S0_dwi  <- function(k, theta, p) p^((k + 1)^theta)
S0_dmw  <- function(k, theta, p, c) p^(((k + 1)^theta) * c^(k + 1))

S_from_S0 <- function(k, x, beta, S0fun, ...) S0fun(k, ...) / exp(k * x * beta)
haz_from_S <- function(k, x, beta, S0fun, ...) {
  Skm1 <- S_from_S0(k - 1, x, beta, S0fun, ...)
  Sk   <- S_from_S0(k,     x, beta, S0fun, ...)
  1 - ifelse(Skm1 <= 0, 0, Sk / Skm1)
}

## 2) Event times: first success of step-wise Bernoulli draws with hazard r(k|x)
simulate_by <- function(n, beta, max_time, px, S0fun, ...) {
  x <- rbinom(n, 1, px); time <- integer(n); event <- integer(n)
  for (i in seq_len(n)) {
    xi <- x[i]; k <- 0L; hit <- FALSE
    while (k <= max_time) {
      r_k <- pmin(pmax(haz_from_S(k, xi, beta, S0fun, ...), 0), 1)
      if (runif(1) < r_k) { time[i] <- k; event[i] <- 1L; hit <- TRUE; break }
      k <- k + 1L
    }
    if (!hit) { time[i] <- max_time; event[i] <- 0L }
  }
  data.frame(time = time, event = event, covar = x)
}

## 3) Censoring to the target proportion: censor the latest events, time ~ floor(U * t)
apply_censor_exact <- function(dat, target_censor = 0.4, discrete = TRUE) {
  n <- nrow(dat)
  n_cens_cur <- sum(dat$event == 0L)
  n_cens_target <- floor(n * target_censor)
  if (n_cens_cur < n_cens_target) {
    idx_event <- which(dat$event == 1L)
    n_extra <- min(n_cens_target - n_cens_cur, length(idx_event))
    idx_extra <- idx_event[order(dat$time[idx_event], decreasing = TRUE)[seq_len(n_extra)]]
    dat$event[idx_extra] <- 0L
    if (discrete) dat$time[idx_extra] <- pmax(0L, floor(runif(length(idx_extra)) * dat$time[idx_extra]))
  }
  attr(dat, "censor_real") <- mean(dat$event == 0L)
  dat
}

## 4) Maximum likelihood fits (log guarded against 0)
safe_ll <- function(f, s, data) -sum(data$event * log(pmax(f, .Machine$double.eps)) +
                                       (1 - data$event) * log(pmax(s, .Machine$double.eps)))
neglog_geom <- function(par, data) {
  p <- par[1]; beta <- par[2]
  s0 <- function(t) (1 - p)^(t + 1); s <- function(t, x) s0(t) / exp(t * x * beta)
  f <- function(t, x) (s0(t - 1) * exp(x * beta) - s0(t)) / exp(t * x * beta)
  safe_ll(f(data$time, data$covar), s(data$time, data$covar), data)
}
neglog_dwi <- function(par, data) {
  theta <- par[1]; p <- par[2]; beta <- par[3]
  s0 <- function(t) p^((t + 1)^theta); s <- function(t, x) s0(t) / exp(t * x * beta)
  safe_ll(s(data$time - 1, data$covar) - s(data$time, data$covar), s(data$time, data$covar), data)
}
neglog_dmw <- function(par, data) {
  theta <- par[1]; p <- par[2]; cc <- par[3]; beta <- par[4]
  s0 <- function(t) p^(((t + 1)^theta) * cc^(t + 1)); s <- function(t, x) s0(t) / exp(t * x * beta)
  safe_ll(s(data$time - 1, data$covar) - s(data$time, data$covar), s(data$time, data$covar), data)
}
fit_beta <- function(fn, init, dat) tryCatch(optim(init, fn = fn, data = dat)$par[length(init)],
                                             error = function(e) NA_real_)

## 5) Scenarios
max_time <- 20
n_set    <- c(50, 200)
beta_set <- c(0.04, 0.08, 0.16)
cens_set <- c(0.2, 0.4)
gens     <- c("Geom", "DMW")
nsim     <- 100
BL <- list(Geom = list(p = 0.04), DMW = list(theta = 1.0, p = 0.99, c = 0.99))

rows <- list()
for (gen in gens) for (n_ in n_set) for (cens in cens_set) for (beta_true in beta_set) {
  for (i in 1:nsim) {
    dat <- switch(gen,
      Geom = simulate_by(n_, beta_true, max_time, px = 0.5, S0fun = S0_geom, p = BL$Geom$p),
      DMW  = simulate_by(n_, beta_true, max_time, px = 0.5, S0fun = S0_dmw,
                         theta = BL$DMW$theta, p = BL$DMW$p, c = BL$DMW$c))
    dat <- apply_censor_exact(dat, target_censor = cens)
    rows[[length(rows) + 1]] <- data.frame(
      gen = gen, n = n_, censor_target = cens, censor_real = attr(dat, "censor_real"),
      beta_true = beta_true, sim = i,
      Geom = fit_beta(neglog_geom, c(0.04, 0.08), dat),
      DWI  = fit_beta(neglog_dwi,  c(0.8, 0.99, 0.08), dat),
      DMW  = fit_beta(neglog_dmw,  c(1.0, 0.99, 0.99, 0.08), dat))
  }
}
sim <- bind_rows(rows)
write.csv(sim, "results/simulation_estimates.csv", row.names = FALSE)

## 6) Bias / MSE summaries and figures
long <- sim %>%
  pivot_longer(c(Geom, DWI, DMW), names_to = "Model", values_to = "beta_hat") %>%
  mutate(bias = beta_hat - beta_true, mse = (beta_hat - beta_true)^2,
         obs = 1 - censor_target,
         scenario = sprintf("%s_beta=%.2f_n=%d_obs=%.1f", gen, beta_true, n, obs)) %>%
  arrange(gen, beta_true, desc(n), obs) %>%
  mutate(scenario = factor(scenario, levels = rev(unique(scenario))))

summary_tbl <- long %>% group_by(gen, beta_true, n, censor_target, Model) %>%
  summarise(mean_bias = mean(bias, na.rm = TRUE), median_bias = median(bias, na.rm = TRUE),
            mse = mean(mse, na.rm = TRUE), .groups = "drop")
write.csv(summary_tbl, "results/simulation_summary.csv", row.names = FALSE)

box <- function(var, xlab, xlim) ggplot(long, aes(.data[[var]], scenario)) +
  geom_boxplot(fill = "grey85", outlier.shape = NA, width = 0.6) +
  geom_point(size = 0.6, alpha = 0.5) + facet_grid(. ~ Model) +
  coord_cartesian(xlim = xlim) + labs(x = xlab, y = NULL) + theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 8),
        strip.background = element_rect(fill = "white"), strip.text = element_text(face = "bold"))
ggsave("results/simulation_bias.png", box("bias", expression("Bias of estimated " * beta), c(-0.3, 0.3)),
       width = 10, height = 8, dpi = 150)
ggsave("results/simulation_mse.png", box("mse", expression("MSE of estimated " * beta), c(0, 0.05)),
       width = 10, height = 8, dpi = 150)
cat("simulation done:", nrow(sim), "datasets\n")
