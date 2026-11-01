# =============================================================================
# Software name : macodgent medical data-science research toolkit
# Copyright holder : Dr.Xie
# Version : 1.0.0
# Language : R
# Completion date : 2026-10-30
# Copyright notice : This program is independently developed by the copyright
#                   holder and protected under the Software Protection
#                   Regulations; without written permission it may not be
#                   copied, distributed, rented, or used for commercial purposes.
# -----------------------------------------------------------------------------
# Module : qualitative description & association analysis module (desc-qual)
# This file contains : mc_contingency / mc_corr / mc_tableone
#              mc_contingency covers three independence tests: chi-square,
#              Fisher's exact, and CMH stratified; mc_corr's partial correlation
#              is computed by regressing on controls and taking residuals, an
#              original implementation.
# =============================================================================

#' Contingency table & independence test
#'
#' Corresponds to Chapter 4, Section 4.3 (chi-square, Fisher's exact, and the
#' Cochran-Mantel-Haenszel test). Builds a contingency table and runs one of
#' three independence tests selected by `method`.
#'
#' @details
#' Purpose: build contingency tables for categorical variables and run
#' independence/association tests; the first step in analyzing categorical
#' associations, covering the three methods required by Section 4.3.
#' Methods:
#' \itemize{
#'   \item \code{"chisq"} (default): chi-square independence test for 2x2 and
#'   larger tables; \code{correct = TRUE} applies Yates' continuity correction to
#'   2x2 tables to avoid overstating significance in small samples. For 2x2
#'   tables, Fisher's exact method additionally gives the OR and its 95% CI.
#'   \item \code{"fisher"}: \code{stats::fisher.test} exact test for any r x c
#'   table (suitable for small samples or sparse expected counts).
#'   \item \code{"cmh"}: requires a \code{strata} variable; builds a 3-way
#'   contingency table and runs \code{stats::mantelhaen.test} (Cochran-Mantel-
#'   Haenszel) to judge whether the row-column association persists after
#'   controlling for the stratifying factor, returning the common OR.
#' }
#' Arguments: `.data` is a data frame; `row`/`col` are the row/column variable
#' names; `correct` controls chi-square continuity correction; `method` selects
#' the test type; `strata` is the layer variable name for CMH.
#' Result: returns a list with table (contingency table, stratified when 3-way),
#' test (a data frame: test/statistic/df/p.value), and method-specific OR /
#' commonOR.
#' Note: chi-square requires sufficiently large expected counts (generally a
#' sufficient proportion of cells with expected > 5), otherwise use
#' \code{method = "fisher"}. CMH is suitable for estimating "adjusted"
#' association when there is a confounding stratum (e.g. multi-center studies).
#' OR > 1 suggests a positive association; a CI crossing 1 means the association
#' is not significant.
#' @param .data A `data.frame`.
#' @param row Row variable name, categorical.
#' @param col Column variable name, categorical.
#' @param correct Continuity correction for chi-square; default `TRUE` (only
#'   applies to `method = "chisq"`).
#' @param method Test type: `"chisq"` (default), `"fisher"`, or `"cmh"`.
#' @param strata Stratifying variable name for the CMH test.
#' @return A list with `table` (contingency table), `test` (a data frame:
#'   test/statistic/df/p.value); for `method = "chisq"` with a 2x2 table it also
#'   returns `OR` (OR/CI.low/CI.high/fisher.p), for `method = "fisher"` it
#'   returns `OR` (with 95% CI and p), and for `method = "cmh"` it returns
#'   `commonOR` (common odds ratio and its CI).
#' @examples
#' # Chapter 4, Section 4.3: demonstrate three independence tests with birthwt
#' birthwt <- MASS::birthwt
#'
#' ## ---- Chi-square independence test (2x2 with Fisher OR) ----
#' res <- mc_contingency(birthwt, row = "smoke", col = "low")
#' cat("---- Contingency table / table ----\n"); print(res$table)
#' cat("\n---- Chi-square test / chi-square ----\n"); print(res$test, row.names = FALSE)
#' cat("\n---- Odds ratio OR / odds ratio ----\n"); print(res$OR, row.names = FALSE)
#'
#' ## ---- Fisher's exact test (small samples / sparse expected counts) ----
#' res2 <- mc_contingency(birthwt, row = "smoke", col = "low", method = "fisher")
#' print(res2$test, row.names = FALSE); print(res2$OR, row.names = FALSE)
#'
#' ## ---- CMH stratified test: common OR of smoking & low birth weight after controlling for hypertension (ht) ----
#' res3 <- mc_contingency(birthwt, row = "smoke", col = "low",
#'                        method = "cmh", strata = "ht")
#' cat("\n---- CMH test / Cochran-Mantel-Haenszel ----\n")
#' print(res3$test, row.names = FALSE)
#' print(res3$commonOR, row.names = FALSE)
#' @export
mc_contingency <- function(.data, row, col, correct = TRUE,
                           method = c("chisq", "fisher", "cmh"),
                           strata = NULL) {
  method <- match.arg(method)
  if (method == "cmh") {
    if (is.null(strata)) stop("`method = 'cmh'` requires a `strata` variable.")
    if (!requireNamespace("survival", quietly = TRUE))
      stop("`method = 'cmh'` requires the survival package.")
    xtab <- stats::xtabs(~ .data[[row]] + .data[[col]] + .data[[strata]])
    mh <- stats::mantelhaen.test(xtab)
    test <- data.frame(test = "Cochran-Mantel-Haenszel",
                       statistic = as.numeric(mh$statistic),
                       df = as.numeric(mh$parameter),
                       p.value = as.numeric(mh$p.value),
                       stringsAsFactors = FALSE)
    ## Common odds ratio (MH point estimate); its 95% CI is obtained with an
    ## exact conditional logistic model (survival::coxph, method = "exact",
    ## equivalent to clogit), more robust than the MH asymptotic variance;
    ## survival is already in Imports.
    cor <- as.numeric(mh$estimate)
    dn <- dimnames(xtab)
    rows <- list()
    for (kk in seq_len(dim(xtab)[3])) {
      t2 <- xtab[, , kk]
      a <- t2[1, 1]; b <- t2[1, 2]; c <- t2[2, 1]; d <- t2[2, 2]
      rows[[length(rows) + 1]] <- data.frame(stratum = dn[[3]][kk], event = 1, exposure = 1, n = a)
      rows[[length(rows) + 1]] <- data.frame(stratum = dn[[3]][kk], event = 0, exposure = 1, n = b)
      rows[[length(rows) + 1]] <- data.frame(stratum = dn[[3]][kk], event = 1, exposure = 0, n = c)
      rows[[length(rows) + 1]] <- data.frame(stratum = dn[[3]][kk], event = 0, exposure = 0, n = d)
    }
    long <- do.call(rbind, rows)
    long <- long[rep(seq_len(nrow(long)), long$n), c("stratum", "event", "exposure")]
    long$stratum <- factor(long$stratum)
    ci <- tryCatch({
      cl <- survival::coxph(
        survival::Surv(rep(1, nrow(long)), long$event) ~ long$exposure +
          survival::strata(long$stratum), method = "exact")
      as.numeric(stats::confint(cl))
    }, error = function(e) c(NA_real_, NA_real_))
    res <- list(
      table = xtab,
      test = test,
      commonOR = data.frame(commonOR = cor,
                            CI.low = if (length(ci) == 2) exp(ci[1]) else NA_real_,
                            CI.high = if (length(ci) == 2) exp(ci[2]) else NA_real_,
                            stringsAsFactors = FALSE))
    return(res)
  }
  tab <- table(.data[[row]], .data[[col]])
  if (method == "fisher") {
    ft <- stats::fisher.test(tab)
    test <- data.frame(test = "Fisher's exact",
                       statistic = NA_real_,
                       df = NA_real_,
                       p.value = as.numeric(ft$p.value),
                       stringsAsFactors = FALSE)
    a <- tab[1, 1]; b <- tab[1, 2]; c <- tab[2, 1]; d <- tab[2, 2]
    or_point <- if (nrow(tab) == 2 && ncol(tab) == 2) (a * d) / (b * c) else NA_real_
    res <- list(
      table = tab,
      test = test,
      OR = data.frame(OR = or_point,
                      CI.low = as.numeric(ft$conf.int[1]),
                      CI.high = as.numeric(ft$conf.int[2]),
                      fisher.p = as.numeric(ft$p.value),
                      stringsAsFactors = FALSE))
    return(res)
  }
  ## default: chi-square
  ct <- stats::chisq.test(tab, correct = correct)
  res <- list(
    table = tab,
    test = data.frame(test = "Chi-square",
                      statistic = as.numeric(ct$statistic),
                      df = as.numeric(ct$parameter),
                      p.value = as.numeric(ct$p.value),
                      stringsAsFactors = FALSE)
  )
  if (nrow(tab) == 2 && ncol(tab) == 2) {
    a <- tab[1, 1]; b <- tab[1, 2]; c <- tab[2, 1]; d <- tab[2, 2]
    or <- (a * d) / (b * c)
    ft <- stats::fisher.test(tab)
    res$OR <- data.frame(OR = or,
                         CI.low = as.numeric(ft$conf.int[1]),
                         CI.high = as.numeric(ft$conf.int[2]),
                         fisher.p = as.numeric(ft$p.value),
                         stringsAsFactors = FALSE)
  }
  res
}

#' Correlation analysis (with partial correlation)
#'
#' Corresponds to Chapter 4: Pearson / Spearman correlation; optionally compute
#' partial correlation by residualising all variables on the control variable(s).
#'
#' @details
#' Purpose: compute correlation coefficients between variables, supporting
#' Pearson (linear) and Spearman (rank/monotonic), plus partial correlation
#' (after controlling for confounders).
#' Methods: pearson measures linear association, assumes approximate normality;
#' spearman is rank-based, robust to non-linear monotonic relations and outliers;
#' control specifies one or more control variables to compute partial correlation
#' (the net correlation after removing their effect).
#' Arguments: `.data` is a data frame; `vars` is a vector of variable names to
#' correlate; `method` selects the correlation type; `control` is an optional
#' control variable name.
#' Result: returns a correlation matrix (coefficients); for partial correlation
#' returns the controlled matrix.
#' Note: the correlation coefficient only measures association strength and
#' direction, not causation; near 0 does not mean no relation (could be
#' non-linear). Spearman suits ordinal data or clearly skewed data. Partial
#' correlation helps identify "spurious correlation" caused by a common confounder.
#' @param .data A `data.frame`.
#' @param vars Variable names; `NULL` = all numeric columns.
#' @param method `"pearson"` (default) or `"spearman"`.
#' @param control Optional control variable names for partial correlation.
#' @return A long `data.frame` with columns `var1`, `var2`, `r`, `p.value`.
#' @examples
#' # Chapter 4, Section 4.5: correlation of birth weight & pregnancy weight in birthwt (Spearman)
#' birthwt <- MASS::birthwt
#' cat("==== Spearman correlation / Spearman ====\n")
#' # Note: birthwt has tied observations; cor.test switches to the asymptotic
#' # method for p-values and warns, which is normal; we suppress the warning here
#' # to simplify the example output
#' suppressWarnings(print(mc_corr(birthwt, vars = c("bwt","lwt"), method = "spearman"),
#'                          row.names = FALSE))
#' cat("\n==== Partial correlation after controlling for mother's age / partial on age ====\n")
#' print(mc_corr(birthwt, vars = c("bwt","lwt"), control = "age"),
#'       row.names = FALSE)
#' @export
mc_corr <- function(.data, vars = NULL, method = c("pearson", "spearman"),
                    control = NULL) {
  method <- match.arg(method)
  if (is.null(vars)) vars <- names(.data)[sapply(.data, is.numeric)]
  mat <- .data[vars]
  if (!is.null(control)) {
    if (method != "pearson") warning("Partial correlation is based on residuals; method forced to pearson.")
    method <- "pearson"
    for (v in vars) {
      f <- stats::as.formula(paste(v, "~", paste(control, collapse = "+")))
      fit <- stats::lm(f, data = .data)
      mat[[v]] <- stats::residuals(fit)
    }
  }
  n <- length(vars)
  pairs <- utils::combn(n, 2)
  rows <- list()
  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1, k]; j <- pairs[2, k]
    t <- stats::cor.test(mat[[i]], mat[[j]], method = method)
    rows[[length(rows) + 1]] <- data.frame(
      var1 = vars[i], var2 = vars[j],
      r = as.numeric(t$estimate), p.value = as.numeric(t$p.value),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

#' Baseline characteristics table
#'
#' Corresponds to Chapter 4 ("Table One"). Numeric variables are summarised as
#' mean +/- SD; categorical variables as counts and percentages; can be
#' stratified by a grouping variable.
#'
#' @details
#' Purpose: generate the common "Table 1" baseline characteristics table, one row
#' per variable, shown by exposure/group column, ready for the manuscript.
#' Methods: numeric variables summarised as "mean +/- SD" (when approximately
#' normal); categorical variables as counts and percentages.
#' Arguments: `.data` is a data frame; `vars` are the variable names to include;
#' `group` is an optional grouping column (e.g. treatment), summarising per group;
#' `digits` controls decimal places, default 2.
#' Result: returns a tidy data frame, directly usable with knitr::kable() or
#' print() to make a table.
#' Note: p-values for group comparisons are usually not in this table (compute
#' separately with mc_ttest/mc_nonparam/mc_contingency). Too many levels for a
#' categorical variable widens the table; consider merging rare categories first.
#' This is a descriptive table, not a substitute for statistical inference.
#' @param .data A `data.frame`.
#' @param vars Variable names to include.
#' @param group Optional grouping variable name.
#' @param digits Number of decimals; default `2`.
#' @return A `data.frame` with columns `variable`, `level`, `stratum`, `value`.
#' @examples
#' # Chapter 4: birthwt baseline characteristics table (grouped by low birth weight)
#' birthwt <- MASS::birthwt
#' bw <- birthwt
#' bw$race  <- factor(bw$race, levels = 1:3, labels = c("White","Black","Other"))
#' bw$smoke <- factor(bw$smoke, levels = 0:1, labels = c("No","Yes"))
#' bw$low   <- factor(bw$low, levels = 0:1, labels = c("Normal","Low"))
#' mc_tableone(bw, vars = c("age","lwt","race","smoke","bwt"), group = "low")
#' @export
mc_tableone <- function(.data, vars, group = NULL, digits = 2) {
  rows <- list()
  for (v in vars) {
    x <- .data[[v]]
    if (is.numeric(x)) {
      rows[[length(rows) + 1]] <- data.frame(
        variable = v, level = "", stratum = "Overall",
        value = sprintf("%.*f +/- %.*f", digits, mean(x, na.rm = TRUE),
                        digits, stats::sd(x, na.rm = TRUE)),
        stringsAsFactors = FALSE)
      if (!is.null(group)) for (lv in levels(factor(.data[[group]]))) {
        xx <- x[.data[[group]] == lv]
        rows[[length(rows) + 1]] <- data.frame(
          variable = v, level = "", stratum = as.character(lv),
          value = sprintf("%.*f +/- %.*f", digits, mean(xx, na.rm = TRUE),
                          digits, stats::sd(xx, na.rm = TRUE)),
          stringsAsFactors = FALSE)
      }
    } else {
      x <- factor(x)
      tb <- table(x)
      for (lvl in names(tb)) {
        pct <- 100 * tb[lvl] / sum(tb)
        rows[[length(rows) + 1]] <- data.frame(
          variable = v, level = lvl, stratum = "Overall",
          value = sprintf("%d (%.1f%%)", tb[lvl], pct),
          stringsAsFactors = FALSE)
      }
      if (!is.null(group)) for (g in levels(factor(.data[[group]]))) {
        xx <- factor(x[.data[[group]] == g])
        tb <- table(xx)
        for (lvl in names(tb)) {
          pct <- 100 * tb[lvl] / sum(tb)
          rows[[length(rows) + 1]] <- data.frame(
            variable = v, level = lvl, stratum = as.character(g),
            value = sprintf("%d (%.1f%%)", tb[lvl], pct),
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  do.call(rbind, rows)
}
