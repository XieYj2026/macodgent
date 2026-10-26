# macodgent NEWS

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
