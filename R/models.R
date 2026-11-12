# Additive alternative discrete hazard model
#   S(k|x) = S0(k) * exp(-k * x * beta),  f(k|x) = S(k-1|x) - S(k|x)
# Parameter order follows the thesis tables: theta, p, c, beta.

baseline <- list(
  geom   = function(t, q) (1 - q[["p"]])^(t + 1),
  # Thesis Table 4.6 (CML) fitted the geometric model with S0(k) = (1-p)^k; kept to reproduce it.
  geom_k = function(t, q) (1 - q[["p"]])^t,
  dw     = function(t, q) q[["p"]]^((t + 1)^q[["theta"]]),
  dmw    = function(t, q) q[["p"]]^(((t + 1)^q[["theta"]]) * q[["c"]]^(t + 1))
)

par_names <- list(geom = c("p", "beta"), geom_k = c("p", "beta"),
                  dw = c("theta", "p", "beta"), dmw = c("theta", "p", "c", "beta"))

surv <- function(model, t, x, q) baseline[[model]](t, q) / exp(t * x * q[["beta"]])

neg_loglik <- function(par, model, data) {
  q <- setNames(par, par_names[[model]])
  s <- surv(model, data$time, data$covar, q)
  f <- surv(model, data$time - 1, data$covar, q) - s
  -sum(data$event * log(f) + (1 - data$event) * log(s))
}

# Nelder-Mead from the given start values, as in the thesis code.
fit_model <- function(model, data, init) {
  o <- suppressWarnings(optim(init, neg_loglik, model = model, data = data))
  k <- length(init); n <- nrow(data); ll <- -o$value
  est <- setNames(o$par, par_names[[model]])
  data.frame(model = model,
             theta = unname(est["theta"]), p = unname(est["p"]),
             c = unname(est["c"]), beta = unname(est["beta"]),
             LogLik = ll, AIC = 2 * k - 2 * ll, BIC = log(n) * k - 2 * ll,
             convergence = o$convergence)
}

fit_all <- function(data, inits) {
  do.call(rbind, lapply(names(inits), function(m) fit_model(m, data, inits[[m]])))
}
