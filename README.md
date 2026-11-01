# macodgent

**MA Codgent Toolkit for Routine Medical Data Science Analysis**

A curated collection of reusable analysis functions implementing the Routine
Analysis Techniques for medical-data-science analysis. Every exported function uses the `mc_` prefix and returns a tidy
`data.frame`, so results drop straight into reports and downstream pipelines.

---

## Design principles

1. Every exported function is prefixed with `mc_` to avoid name clashes.
2. Analysis functions return clean `data.frame` objects, ready for reporting
   and piping into downstream workflows.
3. Heavy methods live in `Suggests` with graceful namespace guards; the core
   depends only on `tidyr`, `survival`, `ggplot2`, `broom`, `stats`, `utils`.

---

## Installation

```r
# Install from source
# R CMD INSTALL macodgent
# or with devtools:
# devtools::install("macodgent")

# Load
library(macodgent)
```

---

## Modules & function index

33 exported functions in 8 modules.

### Meet the data

| Function       | Purpose                                              |
| -------------- | ---------------------------------------------------- |
| `mc_mode()`    | Mode, supports multiple modes                        |
| `mc_inspect()` | Per-column data overview (type / missing / unique)   |
| `mc_theme()`   | MA Codgent light ggplot2 theme                       |
| `mc_palette()` | Brand colour palette                                 |

### Preprocessing

| Function         | Purpose                                          |
| ---------------- | ------------------------------------------------ |
| `mc_impute_knn()`  | KNN imputation (self-implemented, no extra deps) |
| `mc_impute_mice()` | Multiple imputation (MICE wrapper)               |
| `mc_dedup()`       | Remove duplicates (by key columns or all)        |
| `mc_scale()`       | Feature scaling (z / min-max / robust)           |
| `mc_wide_long()`   | Wide to long (tidyr wrapper)                     |
| `mc_long_wide()`   | Long to wide (tidyr wrapper)                     |

### Quantitative description

| Function         | Purpose                                                  |
| ---------------- | -------------------------------------------------------- |
| `mc_desc_quant()` | One-shot ~20 metrics, supports grouping                  |
| `mc_norm_test()`  | Normality test (Shapiro-Wilk / KS)                       |
| `mc_var_test()`   | Variance homogeneity (Bartlett / Fligner / Levene)       |
| `mc_ttest()`      | Unified t-test (one-sample / independent / paired)       |
| `mc_anova()`      | One-way ANOVA                                            |
| `mc_nonparam()`   | Non-parametric group tests (Wilcoxon / Kruskal-Wallis)   |

### Qualitative description

| Function          | Purpose                               |
| ----------------- | ------------------------------------- |
| `mc_contingency()` | Contingency table + chi-square (2x2 OR) |
| `mc_corr()`        | Correlation (with partial)            |
| `mc_tableone()`    | Baseline Table One (supports strata)  |

### Conventional regression

| Function               | Purpose                                    |
| --------------------- | ------------------------------------------ |
| `mc_logistic()`       | Binary logistic regression (OR + CI)       |
| `mc_logistic_screen()`| Univariate -> multivariable (screen first) |
| `mc_stepwise()`       | Stepwise selection (AIC based)             |
| `mc_ordinal_logistic()`| Ordinal logistic (MASS::polr)            |

### Survival analysis

| Function      | Purpose                                  |
| ------------- | ---------------------------------------- |
| `mc_survfit()`| Kaplan-Meier + log-rank + plot          |
| `mc_cox()`    | Cox PH model + forest plot              |

### Advanced regression

| Function           | Purpose                               |
| ----------------- | ------------------------------------- |
| `mc_gam()`        | Generalized additive model (mgcv)     |
| `mc_lasso()`      | Lasso regression (glmnet, alpha=1)    |
| `mc_ridge()`      | Ridge regression (glmnet, alpha=0)    |
| `mc_elastic_net()`| Elastic net (glmnet)                  |
| `mc_rcs()`        | Restricted cubic spline (self-impl.)  |
| `mc_quantile_reg()`| Quantile regression (quantreg)       |

### Advanced survival

| Function             | Purpose                              |
| ------------------- | ------------------------------------ |
| `mc_competing_risk()`| Competing-risk (Fine-Gray, cmprsk)  |
| `mc_multistate()`    | Multi-state model (mstate)          |

---

## Quick start

```r
library(macodgent)
library(MASS)   # ships the birthwt example dataset

# 1) Data overview
mc_inspect(birthwt)

# 2) One-shot quantitative description
mc_desc_quant(birthwt, vars = c("age", "lwt", "bwt"))
mc_desc_quant(birthwt, vars = c("age", "lwt"), group = "smoke")

# 3) Baseline Table One, stratified by smoke
mc_tableone(birthwt, vars = c("age", "lwt", "smoke", "low"), group = "smoke")

# 4) Binary logistic regression (OR + 95% CI)
bw <- birthwt
bw$low <- factor(bw$low)
mc_logistic(bw, low ~ age + smoke)

# 5) Univariate screen -> multivariable
mc_logistic_screen(bw, y = "low",
  x_vars = c("age", "lwt", "race", "smoke", "ht", "ui"))

# 6) Restricted cubic spline (self-implemented basis, no rms needed)
rcs <- mc_rcs(birthwt, x = "age", y = "bwt", nk = 3)
rcs$plot

# 7) Cox model + forest plot
fit <- mc_cox(survival::lung, survival::Surv(time, status) ~ age + sex + ph.ecog)
fit$forest
```

---

## Reading the help

Every function ships with a detailed help page (description, arguments,
value, runnable examples, see-also, references).

```r
?mc_desc_quant          # help for one function
help(package = "macodgent")   # package index
```

Press F1 on a function name in RStudio.

---

## License

MIT (c) MA Codgent OPC
