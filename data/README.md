# Data

| Dataset | File | Included | Source |
|---|---|:---:|---|
| HALLUCA (Halle lung carcinoma) | `halluca.csv` | Yes | `halluca_short.csv` from Voeltz D, Hoyer A (2022), *Parametric additive hazard model for time-to-event analysis*, Zenodo, [doi:10.5281/zenodo.7124989](https://doi.org/10.5281/zenodo.7124989), CC BY 4.0. Used in Voeltz et al. (2024), *BMC Med Res Methodol* 24:48. 1,696 patients; 1,360 after excluding missing TNM stage. |
| CLRC (colorectal cancer, National Cancer Center Korea) | `private/clrc.csv` | No | Downloaded from the National Cancer Center public data portal. Download it there and save it at this path. |
| HIRA CML | `private/cml.xlsx` | No | Derived from Korean nationwide claims data (HIRA). Not redistributable. |

`private/` is git-ignored. Scripts skip a dataset whose file is missing.

Expected columns:

- `halluca.csv`: `t` (years), `TNMII`, `TNMIIIa`, `TNMIIIb`, `TNMIV`, `TNMmis`, `id_count`, `d` (1 = death)
- `clrc.csv`: `t` (days), `TNMIV`, `d`
- `cml.xlsx`: `medicine` (imatinib / nilotinib / dasatinib), `death`, `survival` (months)
