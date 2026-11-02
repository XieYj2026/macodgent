# macodgent

**MA Codgent Toolkit for Routine Medical Data Science Analysis**

A curated collection of reusable analysis functions implementing the Routine
Analysis Techniques for medical-data-science analysis. Every exported function
uses the `mc_` prefix, and every result is shaped so it can go straight into a
report or into the next step of a pipeline.

Current release: **1.0.2** — 39 exported functions in 8 modules.

---

## Design principles

1. Every exported function is prefixed with `mc_` to avoid name clashes.
2. Results are shaped for reporting. Descriptive and testing functions (e.g.
   `mc_desc_quant()`, `mc_ttest()`, `mc_contingency()`, `mc_corr()`,
   `mc_tableone()`, `mc_model_compare()`) return tidy `data.frame` objects.
   Modelling functions return a `list` that holds the fitted model, a tidy
   coefficient table and — where it makes sense — a `ggplot` object, so the same
   call serves both reporting and inspection.
3. Heavy methods live in `Suggests` with graceful namespace guards; the core
   depends only on `survival`, `ggplot2`, `broom`, `stats` and `utils`.
4. Where it lowers the barrier without adding a dependency, the algorithm is
   implemented in-house rather than wrapped: the restricted cubic spline basis in
   `mc_rcs()` (so no `rms::datadist()` step is needed), principal components
   regression in `mc_pcr()`, Levene's test and skewness/kurtosis in the
   descriptive functions, VIF and bootstrap confidence intervals in `mc_linear()`,
   and wide/long reshaping in `mc_wide_long()` / `mc_long_wide()`.

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

39 exported functions in 8 modules.

### Meet the data

| Function       | Purpose                                              |
| -------------- | ---------------------------------------------------- |
| `mc_mode()`    | Mode, supports multiple modes                        |
| `mc_inspect()` | Per-column data overview (type / missing / unique)   |
| `mc_theme()`   | MA Codgent light ggplot2 theme                       |
| `mc_palette()` | Brand colour palette (with interpolated ramps)       |

### Preprocessing

| Function         | Purpose                                                    |
| ---------------- | ---------------------------------------------------------- |
| `mc_impute_knn()`| KNN imputation (self-implemented, no extra deps)           |
| `mc_impute_mice()`| Multiple imputation (MICE wrapper)                        |
| `mc_dedup()`     | Remove duplicates (by key columns or all)                  |
| `mc_scale()`     | Feature scaling (z / min-max / robust)                     |
| `mc_wide_long()` | Wide to long (self-implemented, no pivoting package)       |
| `mc_long_wide()` | Long to wide (self-implemented, requires `id_cols`)        |

### Quantitative description

| Function         | Purpose                                                    |
| ---------------- | ---------------------------------------------------------- |
| `mc_desc_quant()`| One-shot ~20 metrics, supports grouping                    |
| `mc_norm_test()` | Normality test (Shapiro-Wilk / KS)                         |
| `mc_var_test()`  | Variance homogeneity (Bartlett / Fligner / Levene)         |
| `mc_ttest()`     | Unified t-test (one-sample / independent / paired)         |
| `mc_anova()`     | One-way ANOVA                                              |
| `mc_nonparam()`  | Non-parametric group tests (Wilcoxon / Kruskal-Wallis)     |

### Qualitative description

| Function          | Purpose                                                  |
| ----------------- | -------------------------------------------------------- |
| `mc_contingency()`| Contingency table + chi-square / Fisher / CMH (2x2 OR)   |
| `mc_corr()`       | Correlation, Pearson / Spearman, incl. partial           |
| `mc_tableone()`   | Baseline Table One (grouped; descriptive only, no p)     |

### Conventional regression

| Function                | Purpose                                          |
| ----------------------- | ------------------------------------------------ |
| `mc_linear()`           | Linear regression + bootstrap CI, VIF, 4 diag plots |
| `mc_logistic()`         | Binary logistic regression (OR + CI)             |
| `mc_logistic_screen()`  | Univariate -> multivariable (screen first)       |
| `mc_stepwise()`         | Stepwise selection (AIC based)                   |
| `mc_ordinal_logistic()` | Ordinal logistic (MASS::polr)                    |
| `mc_clogistic()`        | Conditional logistic (matched case-control)      |
| `mc_model_compare()`    | AIC / BIC / AICc / Akaike weights across models  |

### Survival analysis

| Function        | Purpose                                          |
| --------------- | ------------------------------------------------ |
| `mc_survfit()`  | Kaplan-Meier + log-rank + plot                   |
| `mc_cox()`      | Cox PH model + forest plot                       |
| `mc_ph_test()`  | Proportional-hazards assumption check            |

### Advanced regression

| Function           | Purpose                                        |
| ------------------ | ---------------------------------------------- |
| `mc_gam()`         | Generalized additive model (mgcv)              |
| `mc_lasso()`       | Lasso regression (glmnet, alpha=1)             |
| `mc_ridge()`       | Ridge regression (glmnet, alpha=0)             |
| `mc_elastic_net()` | Elastic net (glmnet)                           |
| `mc_rcs()`         | Restricted cubic spline (self-impl. basis)     |
| `mc_quantile_reg()`| Quantile regression (quantreg)                 |
| `mc_pcr()`         | Principal components regression (self-impl.)   |
| `mc_nnet_reg()`    | Neural-network regression (nnet)               |

### Advanced survival

| Function               | Purpose                              |
| ---------------------- | ------------------------------------ |
| `mc_competing_risk()`  | Competing-risk (Fine-Gray, cmprsk)   |
| `mc_multistate()`      | Multi-state model (mstate)           |

---

## Quick start

```r
library(macodgent)
library(MASS)      # ships the birthwt example dataset
library(survival)  # ships the lung example dataset

# 1) Data overview: type / n / missing / unique per column
mc_inspect(birthwt)

# 2) One-shot quantitative description (~20 metrics), grouped by smoking
mc_desc_quant(birthwt, vars = c("bwt", "age", "lwt"), group = "smoke")

# 3) Baseline Table 1.
#    mc_tableone() is descriptive only; run mc_ttest() / mc_contingency()
#    separately for group-comparison p-values.
bw <- birthwt
bw$race  <- factor(bw$race,  levels = 1:3, labels = c("white", "black", "other"))
bw$smoke <- factor(bw$smoke, levels = 0:1, labels = c("no", "yes"))
bw$low   <- factor(bw$low,   levels = 0:1, labels = c("normal", "low"))
mc_tableone(bw, vars = c("age", "lwt", "race", "smoke", "bwt"), group = "low")

# 4) Linear regression with self-implemented diagnostics
#    (bootstrap CI + VIF + four diagnostic plots)
lin <- mc_linear(birthwt, bwt ~ age + lwt + smoke + factor(race), bootstrap = 200)
lin$coef        # coefficient table
lin$vif         # multicollinearity
lin$plot        # residuals vs fitted

# 5) Binary logistic regression (OR + 95% CI)
mc_logistic(bw, low ~ age + smoke)

# 6) Univariate screen -> multivariable logistic
screen <- mc_logistic_screen(bw, y = "low",
                             x_vars = c("age", "lwt", "race", "smoke", "ht", "ui"))
screen$selected      # variables passing the univariate screen
screen$multivariable # multivariable coefficients

# 7) Restricted cubic spline (self-implemented basis, no rms needed)
rcs <- mc_rcs(birthwt, x = "age", y = "bwt", nk = 4)
rcs$plot

# 8) Cox model + forest plot
cox <- mc_cox(survival::lung, Surv(time, status) ~ age + sex + ph.ecog)
cox$forest
```

---

## Reading the help

Every function ships with a detailed help page (description, arguments, value,
runnable examples, see-also, references).

```r
?mc_desc_quant          # help for one function
help(package = "macodgent")   # package index
```

Press F1 on a function name in RStudio. See `NEWS.md` for the changelog.

---

## License

MIT (c) MA Codgent OPC
