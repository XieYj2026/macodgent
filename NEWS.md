# macodgent NEWS

## macodgent 1.0.1 (2026-11-01)

- Patch release aligning all help-page examples (`\examples`) with the teaching
  datasets from the column *Exploring Medical Data Science with R*
  ("Routine Analysis Techniques"): `birthwt` (MASS), `ovarian` (survival),
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
- A curated toolkit of reusable analysis functions distilled from the second
  part ("Routine Analysis Techniques") of the column *Exploring Medical Data
  Science with R*.
- Eight modules: meet the data, preprocessing, quantitative / qualitative
  descriptive statistics, conventional and advanced regression, survival
  analysis, and advanced survival analysis.
- All exported functions use the `mc_` prefix and return tidy `data.frame`
  objects for direct use in reports and downstream workflows.
- 33 exported functions, backed by a `testthat` unit-test suite.
- Layered dependencies: core depends only on `tidyr`, `survival`, `ggplot2`,
  `broom`, `stats`, and `utils`; heavier methods (GAM, Lasso, competing-risk,
  multistate, …) sit in `Suggests` with `requireNamespace()` guards for graceful
  degradation.
- Fixes addressing the CRAN incoming pre-test NOTEs on `0.1.0`:
  - `cran-comments.md` is now listed in `.Rbuildignore`, so it is no longer a
    non-standard file at the top level of the built tarball. (Its content moved
    to the submission comment field.)
  - `DESCRIPTION` wording tightened so nothing is flagged by the spell check.
