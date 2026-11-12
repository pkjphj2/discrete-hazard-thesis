# Reproducibility report

Re-run of the thesis analyses with this repository (R 4.5.2, Windows 11), compared with the printed thesis.
`analysis/02_applications.R` writes the comparison to `results/reproduction_check.csv`.

## Summary

| Thesis item | Status |
|---|---|
| Figures 3.1-3.3 (PMF / hazard of the baselines) | Identical values |
| Table 4.2 HALLUCA, all three models | Identical (LogLik, AIC, BIC, parameters to 4 decimals) |
| Table 4.4 CLRC, Geometric and DW | Identical LogLik / AIC / BIC |
| Table 4.4 CLRC, DMW | Re-run finds a slightly higher maximum (LogLik -9089.476 vs -9089.943); same conclusion |
| Table 4.6 CML, all three models | Identical, with the Geometric note below |
| Figures 3.4-3.5 simulation (2,400 datasets, seed 202507) | Identical estimates (max difference 4.5e-13) |

## Differences found and their causes

### 1. CML Geometric row uses S0(k) = (1-p)^k

Section 3.2.1 defines the geometric baseline as S0(k) = (1-p)^(k+1), and the HALLUCA and CLRC rows use that form.
Table 4.6 (CML) was fitted with S0(k) = (1-p)^k, the form written in the S(k|x) line under it and in equation (4.7).
`geom_k` in `R/models.R` reproduces Table 4.6 exactly (p = 0.0068, beta = 0.0066, LogLik = -748.864).

With the (k+1) form used for the other datasets the CML Geometric fit is p = 0.0051, beta = 0.0083, LogLik = -764.406,
AIC = 1532.81, BIC = 1544.58. Geometric remains the best model by AIC and BIC (DW 1534.69 / 1552.34, DMW 1534.57 / 1558.11),
so the conclusion of Section 4.3 does not change. `02_applications.R` reports both rows.

### 2. CLRC DMW depends on the optimizer start value

The thesis parameters give LogLik -9089.946 on the same data (printed -9089.943), so data and likelihood agree.
The DMW likelihood is very flat for CLRC: from 80 start values Nelder-Mead stopped anywhere between
-9089.476 and -9091.507. Several different starts converge to the same best point,

| | theta | p | c | beta | LogLik | AIC | BIC |
|---|---:|---:|---:|---:|---:|---:|---:|
| Thesis | 0.6458 | 0.9886 | 1.067 | 0.0334 | -9089.943 | 18187.88 | 18223.49 |
| Re-run (best of 80 starts) | 0.6315 | 0.9883 | 1.0643 | 0.0325 | -9089.476 | 18186.95 | 18222.56 |

so the thesis value is an earlier stopping point of the optimizer; its exact start value could not be recovered.
DMW stays the best model by AIC and BIC and beta changes only in the fourth decimal, so Section 4.2 is unaffected.
This repository uses the start value from the final thesis script, which reaches the best point.

### 3. Transcription differences with identical likelihoods

| Table | Parameter | Thesis | Re-run | LogLik (both) |
|---|---|---:|---:|---:|
| 4.4 CLRC DW | theta | 0.8104 | 0.8109 | -9097.681 |
| 4.6 CML DMW | theta, c | 1.3035, 0.9324 | 1.3025, 0.9323 | -763.285 |

Section 4.2 reports 1,996 deaths in CLRC; the data contain 1,696 (the HALLUCA sample size), most likely a typo.

## Other notes

- The original script did not run as written: `gens <- c("Geom","DMW").` had a trailing period. Fixed here.
- In the original `SF&RS by appropriate distributions.Rmd`, the CLRC curve swapped theta and c
  (`(k+1)^1.067 * 0.6458^(k+1)`). The thesis Figure 4.2 was drawn from the correct form in the `.R` script;
  here every curve is computed from the fitted parameters, so this cannot recur.
- In Figure 3.3 the DMW(q = 0.9, theta = 0.6, c = 0.5) curve has slightly negative hazards for k >= 1
  (S(k) increases), so that parameter combination is not a proper survival distribution. Reproduced as published.
- Table 4.5's 20 sample survival times all match `survival * 30 / 365` of the CML data, which confirms the time conversion.
