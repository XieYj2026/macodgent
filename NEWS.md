# macodgent NEWS

## macodgent 1.1.0 - bundled synthetic dataset & GitHub release (2026-11-08)

- **Bundled a built-in synthetic dataset `sim_htn_trial`** (3000 subjects x 47
  variables) so every function example and vignette runs out of the box with
  `data(sim_htn_trial)`. The data are fully simulated (seed 20261105), carry
  `attr(., "simulated")` / `provenance` markers, and are documented in
  `?sim_htn_trial`; they are for demonstration only and must not be used for
  clinical, regulatory, or epidemiological inference.
- All 62 exported functions (47 canonical names + 15 legacy aliases) now use
  `sim_htn_trial` as their uniform example dataset.
- Documentation (`man/`) and `README.md` / `vignettes/macodgent-intro.Rmd`
  updated to the synthetic dataset with an explicit "SIMULATED DATA" notice.
- This release is the source-of-truth snapshot published to
  `github.com/XieYj2026/macodgent`.

## macodgent 1.0.2 - full optimisation pass (2026-11-05)

### Correctness fixes (P0)

- `mc_rcs()` now builds GLM confidence intervals on the link scale and
  back-transforms through `family(fit)$linkinv`; adding +/- `se.fit` to a
  `type = "response"` prediction was methodologically wrong.
- `mc_nnet_reg()` reports out-of-sample (holdout / CV) RMSE and R2 instead of
  training-set error; `mc_pca_regression()` evaluates on an untouched
  validation set and back-transforms coefficients to the original scale.
- Fixed a dimension error in PCR (`rot %*% beta[-1]`) that made
  `family = "pca"` fail with "non-conformable arguments".
- `mc_logistic(diagnostics = TRUE)` reuses the full `mc_logistic_report()`
  suite rather than returning a bare coefficient table; `mc_linear()` now
  exposes the same shared PASS / REVIEW checks and has direct unit tests.
- `mc_clogistic()` returns scalar, correctly shaped stratum and event checks;
  `mc_ordinal_logistic()` tags intercept rows so they can be filtered.

### New uniform model dispatcher

- `mc_auto_model()` routes 12 families (`linear`, `logistic`, `cox`, `ordinal`,
  `quantile`, `gam`, `rcs`, `pca`, `elastic_net`, `nnet`, `competing_risk`,
  `multistate`) and always returns the same shape: `family`, `formula`,
  `predictors`, `n`, `events`, `coefficients`, `performance`, `checks`,
  `diagnostics`, `fit`, `raw`. A `print.mc_auto_model()` method is registered.
- `mc_auto_analyze(family = ...)` accepts any of those families; the classic
  linear / logistic / Cox inference remains the default when `family = "auto"`.
- The 1228-line automation file was split into `mc-auto-core.R`,
  `mc-auto.R`, `mc-auto-report.R` and `mc-auto-model.R`.

### Shared methodological check layer

- `R/mc-diagnostics.R` centralises PASS / REVIEW verdicts for linear
  (normality, homoscedasticity, VIF, Cook's D, global F), logistic (AUC, EPV,
  VIF, LRT) and Cox (proportional hazards, concordance) models.

### Enhanced functions

- Regression: `mc_logistic_screen()` and `mc_stepwise()` gain adjusted p values,
  family selection and bootstrap stability; `mc_model_compare()` gains LRT and
  a comparability audit; `mc_ordinal_logistic()` gains OR CIs and a
  proportional-odds check.
- Advanced: `mc_gam()` gains a smooth-term table, parametric terms and
  performance; `mc_quantile_reg()` gains CIs, inference and a trend test;
  `mc_rcs()` gains a response-scale curve and nonlinearity p; lasso / ridge /
  elastic net gain CV tables, `lambda.1se` and coefficient selection frequency.
- Survival: `mc_competing_risk()` gains HR CIs, CIF and cause counts;
  `mc_multistate()` gains state-occupation probabilities; PH testing gains
  global and adjusted p values with an explicit verdict.
- Preprocessing: MICE convergence diagnostics, train-only scaling via
  `mc_scale_numeric(train = ...)`, duplicate-key checks and multi-measure
  pivoting, `values_fill` / `values_fn` aggregation; `mc_inspect(detailed =
  TRUE)` adds outlier, constant, distribution and role flags.

### Reporting and export

- `mc_auto_report()` HTML now has sortable table headers and `@media print`
  rules. `mc_auto_export(word = TRUE)` writes `report.doc`, an HTML table that
  Word / WPS open directly, without `officer` or `rmarkdown`.

### Verification

- `R CMD build` + `R CMD check --no-manual`: **Status: OK** (0 errors,
  0 warnings, 0 notes). `testthat`: 203 assertions in 57 test blocks, 0
  failures, 0 warnings, 0 skips. All 62 exported names: examples run with
  0 errors and every name has a help page; all 47 current names also have a
  direct unit test. Number of exports rises from 61 to 62 (47 current names +
  15 backward-compatible aliases); `man/mc_auto_model.Rd` added.

## macodgent 1.0.2 - one-call automated analysis (2026-11-04)

### New automated workflow

- `mc_auto_analyze()` - one call returns 13 components: data overview,
  per-column audit with automatic role assignment, descriptive statistics
  with distribution verdicts, a stratified baseline table, group comparisons,
  correlation analysis, univariable screening, a multivariable model, model
  diagnostics, plain-language findings, a decision log, caveats and the
  cleaned data. The model family is inferred from the outcome (linear /
  logistic, or Cox from `time` + `status`).
- Automatic safeguards: constant / mostly-missing / identifier-like /
  high-cardinality columns are excluded; numeric 0/1 columns become factors;
  survival status is recoded to 0/1; outcome-derived predictors (association
  >= 0.98) and collinear numeric predictors are dropped; model size is capped
  by events-per-variable; rank-deficient fits are retried with fewer
  predictors. Every choice appears in `decisions`.
- Diagnostics: linear (residual normality, Breusch-Pagan, VIF, Cook's D,
  global F), logistic (via `mc_logistic_report()`: AUC, calibration, EPV, VIF,
  influence, LRT), Cox (global and per-covariate proportional hazards,
  concordance, EPV, LRT). Each check carries a PASS / REVIEW verdict.
- `mc_auto_report()` - dependency-free, self-contained HTML report with KPI
  cards, result tables, a CSS forest plot and diagnostic badges. Also
  available directly via `mc_auto_analyze(report = "...")`.
- `mc_auto_export()` - writes every table as CSV plus `analysis.rds`.
- `print.mc_auto_analysis()` - compact console summary.

## macodgent 1.0.2 - analysis layer and naming (2026-11-04)

### Naming

15 ambiguous function names were replaced by self-explanatory names. The old
names remain exported as exact aliases, so existing code is unaffected (see
`?mc_legacy_aliases`). The most consequential changes are `mc_pcr()` to
`mc_pca_regression()` (PCR reads as polymerase chain reaction in a medical
context), `mc_mode()` to `mc_stat_mode()` (base `mode()` means storage type),
`mc_survfit()` to `mc_kaplan_meier()` (the old name duplicated
`survival::survfit` for a different object), and `mc_wide_long()` /
`mc_long_wide()` to `mc_pivot_longer()` / `mc_pivot_wider()`.

### New functions

- `mc_compare_groups()` - one-call group comparison that selects an appropriate
  test per variable from its type and distribution, and returns an effect size
  with a confidence interval, the raw p value and a Holm-adjusted p value for
  the whole family of tests.
- `mc_adjust_p()` - Bonferroni / Holm / BH / BY correction with significance
  flags.
- `mc_posthoc_tukey()` - Tukey HSD pairwise contrasts after a significant
  one-way ANOVA.
- `mc_logistic_report()` - full diagnostic suite for a binary logistic model:
  AUC with a bootstrap CI, ROC curve, Hosmer-Lemeshow calibration with decile
  plot, approximate VIF, influential observations, events per variable,
  likelihood-ratio test, Nagelkerke R-squared, Brier score, and an explicit
  pass/review verdict for each methodological check. Implemented in base R,
  without pROC, car or rms.

### Enhanced functions

- `mc_ttest()` now returns Hedges' g (or Cohen's d / dz), its approximate CI,
  group sizes, means and standard deviations alongside the original columns.
- `mc_anova_oneway()` adds eta-squared and omega-squared.
- `mc_test_rank()` adds rank-biserial correlation (with Fisher-z CI) or
  epsilon-squared, plus sample sizes and degrees of freedom.
- `mc_test_normality()` adds skewness, kurtosis, n and a plain-language verdict.
- `mc_test_equal_variance()` adds degrees of freedom and group counts; the
  Levene branch now returns the same columns as the other two methods.
- `mc_describe_numeric()` adds a distribution verdict.
- `mc_correlate()` adds confidence intervals, per-pair n and Holm-adjusted p
  values, and handles missing values by complete cases.
- `mc_crosstab_test()` adds Cramer's V / common OR as an effect-size table and
  expected-count (sparse cell) diagnostics with a recommendation.
- `mc_baseline_table()` adds group-comparison p values, Holm-adjusted p values,
  effect sizes and the standardized mean difference.
- `mc_kaplan_meier()` adds median survival with CI, landmark survival
  probabilities and a numbers-at-risk table; the log-rank result now carries df.
- `mc_cox()` adds model performance (n, events, concordance, likelihood-ratio
  test, AIC, R-squared) and explicit PH/global-fit verdicts.
- Replaced a deprecated ggplot2 `size` aesthetic on lines with `linewidth`.

## macodgent 1.0.2 (2026-11-01)

- Fixed all help-page examples so they run cleanly under `example(mc_xxx)`
  (verified: 39/39 exported functions, 0 errors).
- Restored R constants that a prior batch edit had wrongly lower-cased inside
  examples (`FALSE`/`TRUE`/`NA`/`Inf`/`OR`/`HR`), which caused
  `object 'fALSE' not found` type errors when copying examples from the help page.
- `mc_clogistic` now evaluates `survival::clogit` from within the survival
  namespace; survival 3.8.x's `clogit` calls `coxph()` via
  `eval(call, parent.frame())`, so it failed with "could not find function
  coxph" unless survival was attached. The workaround keeps the package
  namespace clean (no `library(survival)` side effect).
- `mc_impute_knn` example corrected to use the real `iris` column names
  (`Sepal.Length` / `Sepal.Width`) instead of the non-existent lowercase
  `sepal.length` / `sepal.width`, which previously threw
  "replacement has 127 rows, data has 150".
- Documentation audit (2026-11-04): fixed 16 malformed roxygen tags
  (`@param...` x5 and `@param.data` x11) that had left the `.data` and `...`
  arguments undocumented in the help pages; corrected README / DESCRIPTION /
  package-level docs which incorrectly claimed every function returns a
  `data.frame` and that the core depends on `tidyr`; added the 6 missing
  functions to the README index (now 39/39, matching NAMESPACE).
- Packaging: declared the `nnet` dependency used by `mc_nnet_reg`, added
  `importFrom(stats, coef, model.response)` and NSE global declarations,
  removed a broken `[0, 1]` markdown link, and excluded `.github/` plus the
  non-portable Chinese-named HTML manual from the source tarball.
  `R CMD build` + `R CMD check --no-manual` now reports **Status: OK**
  (0 errors / 0 warnings / 0 notes).

## macodgent 1.0.1 (2026-11-01)

- Patch release aligning all help-page examples (`\examples`) with the teaching
  datasets ("Routine Analysis Techniques"): `birthwt` (MASS), `ovarian` (survival),
  `infert` (datasets), `mtcars`, and `iris`. Examples now load data via the
  `pkg::dataset` idiom (e.g. `MASS::birthwt`) so they run cleanly under
  `example(mc_xxx)` with no spurious "data set not found" warning.
- Fixed a real bug in `mc_ph_test`: `cox.zph`'s residual matrix has events as
  rows and covariates as columns, so the function now selects the covariate
  **column** (was selecting by row, causing `arguments imply differing number of
  rows`) and plots the standardized Schoenfeld residuals against (transformed)
  event time correctly.
- Cleaned up example output: `broom` tibbles are coerced to `data.frame` before
  printing (no `row.names = FALSE` warning), and the partial-correlation
  example no longer double-controls a variable.
- `mc_pcr` example simplified to numeric predictors only (PCR requires an
  all-numeric model matrix).

## macodgent 0.1.1 (2026-10-26)

- First release, resubmitted after CRAN incoming pre-test feedback. The earlier
  `0.1.0` tarball was uploaded once but never published on CRAN, so this is the
  version users can actually install.
- A curated toolkit of reusable analysis functions implementing the Routine Analysis Techniques.
- Eight modules: meet the data, preprocessing, quantitative / qualitative
  descriptive statistics, conventional and advanced regression, survival
  analysis, and advanced survival analysis.
- All exported functions use the `mc_` prefix and return tidy `data.frame`
  objects for direct use in reports and downstream workflows.
- 33 exported functions, backed by a `testthat` unit-test suite.
- Layered dependencies: core depends only on `tidyr`, `survival`, `ggplot2`,
  `broom`, `stats`, and `utils`; heavier methods (GAM, Lasso, competing-risk,
  multistate, ...) sit in `Suggests` with `requireNamespace()` guards for graceful
  degradation.
- Fixes addressing the CRAN incoming pre-test NOTEs on `0.1.0`:
  - `cran-comments.md` is now listed in `.Rbuildignore`, so it is no longer a
    non-standard file at the top level of the built tarball. (Its content moved
    to the submission comment field.)
  - `DESCRIPTION` wording tightened so nothing is flagged by the spell check.
