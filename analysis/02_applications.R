# Thesis Chapter 4: fit Geometric / DW / DMW to the three datasets and compare with the thesis tables.
# Run from the repository root: Rscript analysis/02_applications.R
suppressMessages({library(dplyr); library(readxl)})
source("R/models.R"); source("R/plots.R")
dir.create("results", showWarnings = FALSE)

datasets <- list(
  HALLUCA = list(
    file = "data/halluca.csv",
    load = function(f) read.csv(f) %>% filter(TNMmis == 0) %>%
      mutate(time = floor(t), covar = TNMIV, event = d),
    inits = list(geom = c(0.05, 0.05), dw = c(0.9, 0.9, 0.05), dmw = c(1, 0.9, 1, 0.05)),
    best = "dmw"),
  CLRC = list(
    file = "data/private/clrc.csv",
    load = function(f) read.csv(f) %>% mutate(time = floor(t / 365), covar = TNMIV, event = d),
    inits = list(geom = c(0.01, 0.03), dw = c(0.7, 0.99, 0.03), dmw = c(0.7, 0.99, 1.06, 0.03)),
    best = "dmw"),
  CML = list(
    file = "data/private/cml.xlsx",
    load = function(f) read_excel(f) %>% filter(medicine %in% c("imatinib", "nilotinib")) %>%
      mutate(time = floor(survival * 30 / 365), covar = as.numeric(medicine == "imatinib"), event = death),
    # geom_k reproduces Table 4.6; geom is the same model with S0(k) = (1-p)^(k+1) like the other datasets.
    inits = list(geom_k = c(0.05, 0.05), geom = c(0.05, 0.05), dw = c(0.8, 0.9, 0.03),
                 dmw = c(1.3, 0.99, 0.93, 0.008)),
    best = "geom_k")
)

paper <- read.csv("analysis/paper_values.csv")
checks <- list()
for (nm in names(datasets)) {
  ds <- datasets[[nm]]
  if (!file.exists(ds$file)) { message(nm, ": ", ds$file, " not found, skipped (see data/README.md)"); next }
  d <- ds$load(ds$file)
  fits <- fit_all(d, ds$inits)
  write.csv(fits, sprintf("results/table_%s.csv", tolower(nm)), row.names = FALSE)
  cat(sprintf("\n== %s (n = %d, events = %d)\n", nm, nrow(d), sum(d$event)))
  print(fits[, c("model", "theta", "p", "c", "beta", "LogLik", "AIC", "BIC")], digits = 6, row.names = FALSE)
  ggsave(sprintf("results/%s_application.png", tolower(nm)),
         plot_application(fits[fits$model == ds$best, ]), width = 10, height = 5.5, dpi = 300)
  ref <- paper[paper$dataset == nm, ]
  m <- merge(ref, fits, by = "model", suffixes = c("_paper", "_repro"))
  checks[[nm]] <- data.frame(dataset = nm, model = m$model,
    LogLik_paper = m$LogLik_paper, LogLik_repro = round(m$LogLik_repro, 3),
    match = abs(m$LogLik_paper - m$LogLik_repro) < 0.0015 &
      abs(m$beta_paper - round(m$beta_repro, 4)) < 1e-4 & abs(m$p_paper - round(m$p_repro, 4)) < 1e-4)
}
if (length(checks)) {
  chk <- do.call(rbind, checks)
  write.csv(chk, "results/reproduction_check.csv", row.names = FALSE)
  cat("\n== Reproduction check vs thesis tables\n"); print(chk, row.names = FALSE)
}
