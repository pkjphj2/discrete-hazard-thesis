suppressMessages({library(ggplot2); library(patchwork)})

# Estimated survival by covariate group and relative survival exp(-k*beta), from fitted parameters.
plot_application <- function(fit_row, k_max = 50) {
  m <- fit_row$model
  q <- unlist(fit_row[c("theta", "p", "c", "beta")])
  k <- 1:k_max
  surv_df <- rbind(data.frame(k = k, x = "0", S = surv(m, k, 0, q)),
                   data.frame(k = k, x = "1", S = surv(m, k, 1, q)))
  rel_df <- data.frame(k = 0:k_max, value = exp(-(0:k_max) * q[["beta"]]))
  base <- theme_minimal(base_size = 10.5) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          panel.grid.minor = element_blank())
  p1 <- ggplot(surv_df, aes(k, S, color = x)) + geom_line(linewidth = 1) +
    scale_color_manual(values = c("#E41A1C", "#377EB8")) +
    labs(title = sprintf("Estimated S(k|x), %s model", toupper(sub("_k", "", m))),
         x = "Time (k)", y = "Survival probability S(k|x)", color = "x") +
    base + theme(legend.position = c(0.85, 0.85))
  p2 <- ggplot(rel_df, aes(k, value)) + geom_line(linewidth = 1, color = "#4DAF4A") +
    geom_point(color = "#4DAF4A", size = 1.7) +
    labs(title = sprintf("Relative survival exp(-k x %.4f)", q[["beta"]]),
         x = "Time (k)", y = "Relative survival") + base
  p1 | p2
}

# PMF and hazard of a baseline survival curve S(0..kmax) (thesis Figures 3.1-3.3).
pmf_hazard <- function(S) {
  data.frame(pmf = pmax(c(1 - S[1], head(S, -1) - S[-1]), 0),
             hazard = 1 - S / c(1, head(S, -1)))
}

plot_pmf_hazard <- function(specs, model, title, kmax_pmf = 14, kmax_haz = 20) {
  build <- function(kmax) do.call(rbind, lapply(specs, function(s) {
    k <- 0:kmax
    cbind(data.frame(spec = s$label, k = k), pmf_hazard(baseline[[model]](k, s$q)))
  }))
  pmf <- build(kmax_pmf); haz <- build(kmax_haz)
  lv <- vapply(specs, `[[`, "", "label")
  pmf$spec <- factor(pmf$spec, lv); haz$spec <- factor(haz$spec, lv)
  stem <- function(d, y, ylab, ttl) ggplot(d, aes(k, .data[[y]], color = spec, shape = spec)) +
    geom_segment(aes(xend = k, y = 0, yend = .data[[y]]), linewidth = 0.5, alpha = 0.85, show.legend = FALSE) +
    geom_point(size = 2.6) + labs(x = "k", y = ylab, title = ttl, color = NULL, shape = NULL) +
    theme_bw() + theme(plot.title = element_text(face = "bold", hjust = 0.5),
                       panel.grid.minor = element_blank())
  (stem(pmf, "pmf", "Probability mass functions", paste(title, "- PMF")) |
     stem(haz, "hazard", "Hazard functions", paste(title, "- HF"))) +
    plot_layout(guides = "collect") & theme(legend.position = "top")
}
