# macodgent

**MA Codgent Toolkit for Routine Medical Data Science Analysis**

A curated collection of reusable analysis functions distilled from the second
part ("Routine Analysis Techniques") of the column *Exploring Medical Data
Science with R*. Every exported function uses the `mc_` prefix and returns a tidy
`data.frame`, so results drop straight into reports and downstream pipelines.

## Design principles

1. Every exported function is prefixed with `mc_` to avoid name clashes with
   base R and other packages.
2. Analysis functions return clean `data.frame` objects (not just `print`
   output), ready for reporting and piping.
3. Heavy / specialised methods (GAM, Lasso, competing-risk, multistate, …) live
   in `Suggests` and degrade gracefully through namespace guards; the core
   depends only on `tidyr`, `survival`, `ggplot2`, `broom`, `stats`, and
   `utils`.

## Modules and functions

| Chapter | Functions |
|---|---|
| 1 Meet the data | `mc_mode()`, `mc_inspect()`, `mc_theme()`, `mc_palette()` |
| 2 Preprocessing | `mc_impute_knn()`, `mc_impute_mice()`, `mc_dedup()`, `mc_scale()`, wide/long transforms |
| 3 Quantitative description | `mc_desc_quant()` (20 metrics at once), `mc_norm_test()`, `mc_var_test()`, `mc_ttest()`, `mc_anova()`, `mc_nonparam()` |
| 4 Qualitative description | `mc_contingency()` (with OR), `mc_corr()` (with partial correlation), `mc_tableone()` |
| 5 Conventional regression | `mc_logistic()`, `mc_logistic_screen()` (univariate screen → multivariable), `mc_stepwise()`, `mc_ordinal_logistic()` |
| 6 Survival analysis | `mc_survfit()` (KM + log-rank), `mc_cox()` (PH test + forest plot) |
| 7 Advanced regression | `mc_gam()`, `mc_lasso()`, `mc_ridge()`, `mc_elastic_net()`, `mc_rcs()` (self-built spline basis), `mc_quantile_reg()` |
| 8 Advanced survival | `mc_competing_risk()` (Fine-Gray), `mc_multistate()` |

## Installation

```r
# From CRAN (once published)
install.packages("macodgent")

# Development version from source
# R CMD INSTALL macodgent
# or: devtools::install("macodgent")
```

## Quick examples

```r
library(macodgent)
library(MASS)

# One-shot quantitative description (20 metrics)
mc_desc_quant(birthwt, vars = c("age", "lwt", "bwt"))

# Univariate screen -> multivariable logistic
bw <- birthwt
mc_logistic_screen(bw, y = "low",
  x_vars = c("age", "lwt", "race", "smoke", "ht", "ui"))

# Restricted cubic spline (self-built basis, no rms needed)
rcs <- mc_rcs(birthwt, x = "age", y = "bwt", nk = 3)
rcs$plot
```

## License

MIT © MA Codgent OPC
