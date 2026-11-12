# Thesis Figures 3.1-3.3: PMF and hazard of the Geometric, DW and DMW baselines.
# Run from the repository root: Rscript analysis/01_distributions.R
source("R/models.R"); source("R/plots.R")
dir.create("results", showWarnings = FALSE)

spec <- function(label, ...) list(label = label, q = c(...))
figs <- list(
  geometric = plot_pmf_hazard(list(spec("Geo(p=0.2)", p = 0.2), spec("Geo(p=0.4)", p = 0.4),
                                   spec("Geo(p=0.6)", p = 0.6)), "geom", "Geometric"),
  dw  = plot_pmf_hazard(list(spec("DW(q=0.9, θ=0.6)", p = 0.9, theta = 0.6),
                             spec("DW(q=0.9, θ=1)",   p = 0.9, theta = 1.0),
                             spec("DW(q=0.9, θ=1.6)", p = 0.9, theta = 1.6)), "dw", "DW"),
  dmw = plot_pmf_hazard(list(spec("DMW(q=0.9, θ=0.6, c=0.5)",   p = 0.9,  theta = 0.6, c = 0.5),
                             spec("DMW(q=0.95, θ=0.4, c=1.3)",  p = 0.95, theta = 0.4, c = 1.3),
                             spec("DMW(q=0.5, θ=0.5, c=1.05)",  p = 0.5,  theta = 0.5, c = 1.05)), "dmw", "DMW"))
for (nm in names(figs)) ggsave(sprintf("results/%s_pmf_hazard.png", nm), figs[[nm]], width = 10, height = 4, dpi = 300)
cat("saved", length(figs), "figures\n")
