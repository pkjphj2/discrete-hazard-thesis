# A novel discrete parametric survival model based on the alternative discrete hazard rate

Code for my M.S. thesis (Department of Statistics, Hankuk University of Foreign Studies, February 2026).
Full text on [RISS](https://www.riss.kr/search/detail/DetailView.do?p_mat_type=be54d9b8bc7cdb09&control_no=6454e230431be04cffe0bdc3ef48d419).

The thesis proposes an **additive alternative discrete hazard model** for discrete-time survival data:

$$S(k \mid x) = S_0(k)\,\exp(-k\,x^\top\beta), \qquad r(k \mid x) = r_0(k) + x^\top\beta$$

The covariate shifts the alternative hazard additively, so its effect reads as **relative survival**
$\exp(-k\beta)$ that accumulates over time. Three baselines are derived:

| Baseline | $S_0(k)$ | Hazard shape |
|---|---|---|
| Geometric | $(1-p)^{k+1}$ | constant |
| Discrete Weibull (DW) | $p^{(k+1)^\theta}$ | monotone |
| Discrete Modified Weibull (DMW) | $p^{(k+1)^\theta c^{k+1}}$ | monotone or bathtub |

## Repository

```
R/models.R                 baselines, likelihood, fitting (Nelder-Mead, as in the thesis)
R/plots.R                  survival / relative survival and PMF / hazard plots
analysis/01_distributions.R   Figures 3.1-3.3
analysis/02_applications.R    Tables 4.2, 4.4, 4.6 and Figures 4.1-4.3
analysis/03_simulation.R      Figures 3.4-3.5 (seed 202507, ~5 min)
analysis/paper_values.csv     numbers printed in the thesis, used for the reproduction check
data/                         HALLUCA included; see data/README.md for the other two
```

```r
# from the repository root
Rscript analysis/01_distributions.R
Rscript analysis/02_applications.R
Rscript analysis/03_simulation.R
```

Packages: `dplyr`, `tidyr`, `readxl`, `ggplot2`, `patchwork`. Tested with R 4.5.2.

## Reproducibility

Every table and figure in the thesis is reproduced. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the comparison and the three places where the thesis differs from a straight re-run.

## Results at a glance

| Dataset | n | Events | Best model (AIC) | $\hat\beta$ |
|---|---:|---:|---|---:|
| HALLUCA (lung cancer) | 1,360 | 1,111 | DMW | 0.0835 |
| CLRC (colorectal cancer) | 54,311 | 1,696 | DMW | 0.0334 |
| HIRA CML | 2,655 | 134 | Geometric | 0.0066 |

## Citation

Park, H. (2026). *A novel discrete parametric survival model based on the alternative discrete hazard rate*. M.S. thesis, Hankuk University of Foreign Studies.

## License
Code: All rights reserved. 열람만 허용하며, 허락 없이 복제·수정·배포·사용할 수 없습니다. 자세한 내용은 [LICENSE](LICENSE).
Data: `data/halluca.csv` keeps its original license (CC BY 4.0, see [data/README.md](data/README.md)).
