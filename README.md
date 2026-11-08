# macodgent

**MA Codgent Toolkit for Routine Medical Data Science Analysis**

A curated collection of reusable analysis functions implementing the Routine
Analysis Techniques for medical-data-science analysis. Every exported function
uses the `mc_` prefix. Tabular results are returned as tidy `data.frame`
objects; model-fitting functions return a `list` holding the fitted model, a
coefficient table and, where applicable, a plot.

The analysis layer reports **effect sizes, confidence intervals and
multiplicity-adjusted p values** alongside every test, and model functions
ship with explicit methodological checks (discrimination, calibration,
overfitting, collinearity, influence, proportional hazards).

**One call, many results:** `mc_auto_analyze()` runs the whole first-pass
workflow (audit, descriptives, Table One, group comparison, correlation,
univariable screening, multivariable model, diagnostics, written findings and
a decision log). `mc_auto_model()` routes any of 12 model families through one
uniform interface; `mc_auto_report()` and `mc_auto_export()` turn that bundle
into a standalone HTML report or a folder of CSV tables and a Word-readable
`report.doc`.

---

## Design principles

1. Every exported function is prefixed with `mc_` to avoid name clashes.
2. Descriptive and test functions return clean `data.frame` objects, ready for
   reporting. Model-fitting functions return a `list` (fitted model +
   coefficient table + diagnostics + optional plot), so the underlying model
   remains available for further analysis.
3. A p value is never returned alone. Tests also give an effect size with a
   confidence interval; families of tests give adjusted p values.
4. Heavy methods live in `Suggests` with graceful namespace guards; the core
   imports `survival`, `ggplot2`, `broom`, `stats` and `utils`. Reshaping,
   effect sizes, ROC/AUC and all diagnostics are implemented in base R.
5. Names state what the function does. Ambiguous abbreviations were replaced
   in 1.0.2; the old names remain available as exact aliases.
6. Methodological checks are centralised in one shared layer, so linear,
   logistic and Cox models return the same `PASS` / `REVIEW` verdict shape
   instead of three incompatible ad-hoc summaries.
7. Validation metrics are computed out of sample where scaling or tuning is
   involved (holdout / CV split, train-only scaling), so reported performance
   is not optimistically biased.

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

47 current functions in 11 modules, plus 15 backward-compatible aliases
(62 exported names).

### Automated one-call analysis

| Function | Purpose |
| --- | --- |
| `mc_auto_analyze()` | Infer model family, audit columns, describe, compare groups, correlate, screen, fit, diagnose and write findings + a decision log in one call |
| `mc_auto_model()` | Uniform dispatcher for 12 model families, always returning `coefficients`, `performance`, `checks`, `diagnostics`, `fit` and `formula` |
| `mc_auto_report()` | Self-contained light-theme HTML report with a forest plot, diagnostic badges, sortable tables and print CSS |
| `mc_auto_export()` | Write every result table as CSV plus `analysis.rds`; `word = TRUE` adds a Word / WPS readable `report.doc` |

```r
res <- mc_auto_analyze(MASS::birthwt, outcome = "low", group = "low",
                       report = "birthwt-report.html")
res$findings          # plain-language conclusions
res$model$checks      # PASS / REVIEW per methodological check
res$decisions         # every automatic choice, with its reason

# Force a family instead of inferring it, or call the dispatcher directly
pcr <- mc_auto_analyze(mtcars, outcome = "mpg", family = "pca")
fit <- mc_auto_model(mtcars, outcome = "mpg", family = "elastic_net")
fit$performance
print(fit)            # compact console summary

mc_auto_export(res, dir = "birthwt-tables", word = TRUE)
```

The family is inferred from the outcome (linear / logistic, or Cox with
`time` + `status`). `family =` additionally accepts `ordinal`, `quantile`,
`gam`, `rcs`, `pca`, `elastic_net`, `nnet`, `competing_risk` and `multistate`.
Outcome-derived predictors, collinear numeric predictors and model-size
violations are detected and excluded, with each decision logged. Variable
pre-selection is explicitly flagged as exploratory.

### Meet the data

| Function | Purpose |
| --- | --- |
| `mc_stat_mode()` | Statistical mode, supports multiple modes |
| `mc_inspect()` | Per-column overview; `detailed = TRUE` adds outlier count, constant flag, distribution and a suggested analysis role |
| `mc_theme()` | MA Codgent light ggplot2 theme |
| `mc_palette()` | Brand colour palette |

### Preprocessing

| Function | Purpose |
| --- | --- |
| `mc_impute_knn()` | KNN imputation (self-implemented, no extra deps) |
| `mc_impute_mice()` | Multiple imputation (MICE wrapper) + chain-mean convergence diagnostics |
| `mc_dedup()` | Remove duplicates (by key columns or all) |
| `mc_scale_numeric()` | Feature scaling (z / min-max / robust); `train =` reuses training parameters to prevent leakage |
| `mc_pivot_longer()` | Wide to long (base R), duplicate-key check, multiple measure columns |
| `mc_pivot_wider()` | Long to wide (base R), `values_fill` and `values_fn` aggregation |

### Quantitative description

| Function | Purpose |
| --- | --- |
| `mc_describe_numeric()` | ~20 metrics + distribution verdict, supports grouping |
| `mc_test_normality()` | Shapiro-Wilk / KS + skewness, kurtosis, verdict |
| `mc_test_equal_variance()` | Bartlett / Fligner / Levene + group counts |
| `mc_ttest()` | One-sample / independent / paired, with Hedges' g and CI |
| `mc_anova_oneway()` | One-way ANOVA + eta-squared / omega-squared |
| `mc_test_rank()` | Wilcoxon / Kruskal-Wallis + rank-biserial / epsilon-squared |

### Qualitative description

| Function | Purpose |
| --- | --- |
| `mc_crosstab_test()` | Chi-square / Fisher / CMH + Cramer's V + sparse-cell diagnostics |
| `mc_correlate()` | Pearson / Spearman / partial, with CI and adjusted p |
| `mc_baseline_table()` | Table One with p values, SMD and effect sizes |

### Inference & reporting

| Function | Purpose |
| --- | --- |
| `mc_compare_groups()` | One-call group comparison: automatic test selection, effect size + CI, Holm-adjusted p |
| `mc_adjust_p()` | Bonferroni / Holm / BH / BY multiplicity correction with significance flags |
| `mc_posthoc_tukey()` | Tukey HSD pairwise contrasts after ANOVA |
| `mc_logistic_report()` | Full logistic suite: AUC, calibration, VIF, influence, EPV, explicit verdicts |

### Conventional regression

| Function | Purpose |
| --- | --- |
| `mc_linear()` | Linear regression with bootstrap CI, VIF, four diagnostics and shared PASS / REVIEW checks |
| `mc_logistic()` | Binary logistic regression (OR + CI); `diagnostics = TRUE` reuses the full `mc_logistic_report()` suite |
| `mc_logistic_screen()` | Univariate screening then multivariable, with adjusted p and bootstrap selection stability |
| `mc_stepwise()` | AIC stepwise selection for any GLM family, with bootstrap selection stability |
| `mc_ordinal_logistic()` | Ordinal logistic (MASS::polr) with OR 95% CI and an approximate proportional-odds check |
| `mc_clogistic()` | Conditional logistic (matched / stratified) with OR CI and stratum / event checks |
| `mc_model_compare()` | Rank models by AICc weights, with LRT and a comparability audit across families |

### Survival analysis

| Function | Purpose |
| --- | --- |
| `mc_kaplan_meier()` | KM + log-rank + median survival + landmark rates + risk table |
| `mc_cox()` | Cox PH model + forest plot + concordance + PH verdict |
| `mc_test_proportional_hazards()` | Schoenfeld residual test + plot, global p, adjusted per-term p and an explicit verdict |

### Advanced regression

| Function | Purpose |
| --- | --- |
| `mc_gam()` | GAM (mgcv) with smooth-term significance table, parametric terms and performance |
| `mc_lasso()` | Lasso regression (glmnet, alpha=1) with CV table and coefficient stability |
| `mc_ridge()` | Ridge regression (glmnet, alpha=0) with CV table |
| `mc_elastic_net()` | Elastic net (glmnet) with lambda.1se, CV table and stability |
| `mc_rcs()` | Restricted cubic spline with a response-scale curve (correct link-scale CI) and a nonlinearity test |
| `mc_quantile_reg()` | Quantile regression (quantreg) with CI, inference and a cross-quantile trend test |
| `mc_pca_regression()` | PCR with original-scale CI and an untouched external validation set |
| `mc_nnet_reg()` | Neural-network regression (nnet) reporting out-of-sample RMSE / R2 |

### Advanced survival

| Function | Purpose |
| --- | --- |
| `mc_competing_risk()` | Competing-risk (Fine-Gray, cmprsk) with HR 95% CI, CIF and cause counts |
| `mc_multistate()` | Multi-state model (mstate) with state-occupation probabilities and checks |

---

## Example dataset — SIMULATED

**`sim_htn_trial` is simulated data, not real patient data.** It contains
3000 synthetic participants and 47 columns, generated by
`data-raw/sim_htn_trial.R` with `set.seed(20261105)`. It has no real patients,
no real trial results and no personally identifiable information. Every
help-page example, this README and the vignette use this one dataset. All
estimates, p values and survival probabilities shown are simulation artefacts
and must not be used for clinical, regulatory or epidemiological inference.
`attr(sim_htn_trial, "simulated")` is `TRUE`; `attr(..., "provenance")` and
`comment(...)` repeat the warning. See `?sim_htn_trial` for the 47-column
dictionary.

| Variable role | Columns | Used by |
| --- | --- | --- |
| Continuous endpoints | `sbp_week12`, `dbp_week12`, `sbp_baseline` | linear, GAM, RCS, quantile, PCR, penalised, neural-network regression |
| Binary / ordinal endpoints | `bp_control`, `adverse_event`, `bp_grade`, `compliance` | logistic, ordinal logistic, screening, stepwise, crosstabs |
| Time-to-event | `followup_days`, `cv_event`, `comp_event` | Kaplan-Meier, Cox, PH test, competing risk, multistate |
| Grouping / strata | `treatment`, `site`, `race`, `dose_level`, `match_id`, `case_status` | Table One, group comparison, ANOVA, CMH, conditional logistic |
| Repeated measures | `sbp_baseline`, `sbp_week4`, `sbp_week8`, `sbp_week12` | wide-to-long / long-to-wide pivots |
| Missing values | `potassium`, `quality_score`, `adherence_pct`, `urine_albumin` | KNN / MICE imputation, data audit |

## Quick start

```r
library(macodgent)

# SIMULATED DATA -- 3000 synthetic hypertension-trial subjects, 47 columns.
# Not real patients; see ?sim_htn_trial before interpreting any result.
data(sim_htn_trial)
d <- sim_htn_trial
dim(d)
attr(d, "simulated")          # TRUE

# 0) One call: audit, descriptives, comparison, correlation, screening,
#    multivariable model, diagnostics, written findings and a decision log
auto <- mc_auto_analyze(d, outcome = "bp_control", group = "treatment",
  vars = c("age", "bmi", "sbp_baseline", "ldl", "glucose", "smoke", "diabetes"),
  auc_boot = 0L)
auto$findings
auto$model$checks
auto$decisions
mc_auto_export(auto, dir = tempfile(), word = TRUE)

# 1) Data overview and distribution verdict
mc_inspect(d[, c("age", "bmi", "crp", "smoke")], detailed = TRUE)
mc_describe_numeric(d, vars = c("age", "bmi", "sbp_baseline"),
                    group = "treatment")

# 2) Group comparison with effect sizes and multiplicity adjustment
mc_compare_groups(d, group = "treatment",
                  vars = c("age", "bmi", "crp", "smoke", "sex"))

# 3) Table One with p values, adjusted p and SMD
mc_baseline_table(d, vars = c("age", "bmi", "sex", "smoke", "sbp_baseline"),
                  group = "treatment")

# 4) Full logistic report: AUC, calibration, VIF, influence, EPV, verdicts
rep <- mc_logistic_report(d, bp_control ~ age + bmi + treatment + diabetes,
                          auc_boot = 0L)
rep$checks
rep$performance
rep$coefficients

# 5) ANOVA with effect size, then the post-hoc contrasts it requires
mc_anova_oneway(sbp_week12 ~ race, data = d)
mc_posthoc_tukey(data = d, formula = sbp_week12 ~ race)

# 6) Correlation matrix with CIs and Holm-adjusted p values
mc_correlate(d, vars = c("age", "bmi", "sbp_baseline", "ldl", "glucose"))

# 7) Survival: KM with numbers at risk, Cox with concordance and PH verdict
km <- mc_kaplan_meier(d, time = "followup_days", status = "cv_event",
                      group = "treatment")
km$median_survival
km$risk_table
cox <- mc_cox(d, survival::Surv(followup_days, cv_event) ~
                age + bmi + treatment + diabetes)
cox$performance
cox$checks

# 8) One uniform dispatcher for 12 model families
m <- mc_auto_model(d, outcome = "sbp_week12", family = "linear")
print(m)
```

---

## Renamed functions

15 ambiguous names were replaced by self-explanatory ones in 1.0.2. The old
names are **still exported as exact aliases**; existing scripts keep working.

| Superseded | Use instead | Reason |
| --- | --- | --- |
| `mc_mode()` | `mc_stat_mode()` | collided with base `mode()` (storage type) |
| `mc_pcr()` | `mc_pca_regression()` | PCR reads as polymerase chain reaction in medicine |
| `mc_desc_quant()` | `mc_describe_numeric()` | "quant" reads as quantile |
| `mc_norm_test()` | `mc_test_normality()` | "norm" reads as normalisation |
| `mc_var_test()` | `mc_test_equal_variance()` | "var" reads as variable |
| `mc_nonparam()` | `mc_test_rank()` | did not name the test family |
| `mc_anova()` | `mc_anova_oneway()` | did not state the design |
| `mc_survfit()` | `mc_kaplan_meier()` | duplicated `survival::survfit` for a different object |
| `mc_wide_long()` | `mc_pivot_longer()` | direction was unclear |
| `mc_long_wide()` | `mc_pivot_wider()` | direction was unclear |
| `mc_scale()` | `mc_scale_numeric()` | did not state that only numeric columns change |
| `mc_tableone()` | `mc_baseline_table()` | spelled-out clinical term |
| `mc_ph_test()` | `mc_test_proportional_hazards()` | abbreviation |
| `mc_corr()` | `mc_correlate()` | verb-form consistency |
| `mc_contingency()` | `mc_crosstab_test()` | did not state that a test is run |

See `?mc_legacy_aliases` for the complete mapping.

---

## Reading the help

Every function ships with a detailed help page (description, arguments,
value, runnable examples, see-also).

```r
?mc_logistic_report          # help for one function
help(package = "macodgent")   # package index
```

Press F1 on a function name in RStudio.

---

## License

MIT (c) MA Codgent OPC
