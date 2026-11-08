# =============================================================================
# Software name : macodgent medical data-science research toolkit
# Copyright holder : Dr.Xie
# Version : 1.0.0
# Language : R
# Completion date : 2026-10-30
# Copyright notice : This program is independently developed by the copyright
# holder and protected under the Software Protection
# Regulations; without written permission it may not be
# copied, distributed, rented, or used for commercial purposes.
# -----------------------------------------------------------------------------
# Module : qualitative description & association analysis module (desc-qual)
# This file contains : mc_crosstab_test / mc_correlate / mc_baseline_table
# mc_crosstab_test covers three independence tests: chi-square,
# Fisher's exact, and CMH stratified; mc_correlate's partial correlation
# is computed by regressing on controls and taking residuals, an
# original implementation.
# =============================================================================

#' Contingency table & independence test
#'
#' Corresponds to (chi-square, Fisher's exact, and the
#' Cochran-Mantel-Haenszel test). Builds a contingency table and runs one of
#' Three independence tests selected by `method`.
#'
#' @details
#' Purpose: build contingency tables for categorical variables and run
#' Independence/association tests; the first step in analyzing categorical
#' Associations, covering the three methods required by Section 4.3.
#' Methods:
#' \itemize{
#' \item \code{"chisq"} (default): chi-square independence test for 2x2 and
#' Larger tables; \code{correct = TRUE} applies Yates' continuity correction to
#' 2x2 tables to avoid overstating significance in small samples. For 2x2
#' Tables, Fisher's exact method additionally gives the OR and its 95% CI.
#' \item \code{"fisher"}: \code{stats::fisher.test} exact test for any r x c
#' Table (suitable for small samples or sparse expected counts).
#' \item \code{"cmh"}: requires a \code{strata} variable; builds a 3-way
#' Contingency table and runs \code{stats::mantelhaen.test} (Cochran-Mantel-
#' Haenszel) to judge whether the row-column association persists after
#' Controlling for the stratifying factor, returning the common OR.
#' }
#' Arguments: `.data` is a data frame; `row`/`col` are the row/column variable
#' Names; `correct` controls chi-square continuity correction; `method` selects
#' The test type; `strata` is the layer variable name for CMH.
#' Result: returns a list with table (contingency table, stratified when 3-way),
#' Test (a data frame: test/statistic/df/p.value), and method-specific OR /
#' CommonOR.
#' Note: chi-square requires sufficiently large expected counts (generally a
#' Sufficient proportion of cells with expected > 5), otherwise use
#' \code{method = "fisher"}. CMH is suitable for estimating "adjusted"
#' Association when there is a confounding stratum (e.g. multi-center studies).
#' OR > 1 suggests a positive association; a CI crossing 1 means the association
#' Is not significant.
#' @param .data A `data.frame`.
#' @param row Row variable name, categorical.
#' @param col Column variable name, categorical.
#' @param correct Continuity correction for chi-square; default `TRUE` (only
#' Applies to `method = "chisq"`).
#' @param method Test type: `"chisq"` (default), `"fisher"`, or `"cmh"`.
#' @param strata Stratifying variable name for the CMH test.
#' @return A list with `table` (contingency table), `test` (a data frame:
#'   test/statistic/df/p.value), `effect` (Cramer's V or the common OR, with a
#'   CI where available), `diagnostics` (sample size, minimum expected count,
#'   number of sparse cells and a recommendation), plus the method-specific
#'   `OR` or `commonOR` tables.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' mc_crosstab_test(d, row = "treatment", col = "bp_control")
#' mc_crosstab_test(d, row = "treatment", col = "bp_control", method = "fisher")
#' mc_crosstab_test(d, row = "treatment", col = "bp_control",
#'                  method = "cmh", strata = "site")
#' @export
mc_crosstab_test <- function(.data, row, col, correct = TRUE,
 method = c("chisq", "fisher", "cmh"),
 strata = NULL) {
 method <- match.arg(method)
 if (method == "cmh") {
 if (is.null(strata)) stop("`method = 'cmh'` requires a `strata` variable.")
 if (!requireNamespace("survival", quietly = TRUE))
 stop("`method = 'cmh'` requires the survival package.")
 xtab <- stats::xtabs(~.data[[row]] +.data[[col]] +.data[[strata]])
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
 stringsAsFactors = FALSE),
 effect = data.frame(effect = "common odds ratio", estimate = cor,
 CI.low = if (length(ci) == 2) exp(ci[1]) else NA_real_,
 CI.high = if (length(ci) == 2) exp(ci[2]) else NA_real_,
 stringsAsFactors = FALSE),
 diagnostics = data.frame(n = sum(xtab),
 note = "stratified table; inspect stratum-specific counts",
 stringsAsFactors = FALSE))
 return(res)
 }
 tab <- table(.data[[row]],.data[[col]])
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
 stringsAsFactors = FALSE),
 effect = data.frame(effect = if (nrow(tab) == 2 && ncol(tab) == 2)
 "phi / Cramer's V" else "Cramer's V",
 estimate = NA_real_,
 note = "exact test; effect size estimated from the chi-square approximation",
 stringsAsFactors = FALSE),
 diagnostics = .mc_table_diagnostics(tab))
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
 v <-.mc_cramers_v(as.numeric(ct$statistic), sum(tab), nrow(tab), ncol(tab))
 res$effect <- data.frame(
 effect = if (nrow(tab) == 2 && ncol(tab) == 2) "phi / Cramer's V" else "Cramer's V",
 estimate = v, CI.low = NA_real_, CI.high = NA_real_,
 note = "point estimate; no analytic CI without a resampling step",
 stringsAsFactors = FALSE)
 res$diagnostics <-.mc_table_diagnostics(tab, ct$expected)
 res
}

#' Correlation analysis (with partial correlation)
#'
#' Pearson / Spearman correlation; optionally compute
#' Partial correlation by residualising all variables on the control variable(s).
#'
#' @details
#' Purpose: compute correlation coefficients between variables, supporting
#' Pearson (linear) and Spearman (rank/monotonic), plus partial correlation
#' (after controlling for confounders).
#' Methods: pearson measures linear association, assumes approximate normality;
#' Spearman is rank-based, robust to non-linear monotonic relations and outliers;
#' Control specifies one or more control variables to compute partial correlation
#' (the net correlation after removing their effect).
#' Arguments: `.data` is a data frame; `vars` is a vector of variable names to
#' Correlate; `method` selects the correlation type; `control` is an optional
#' Control variable name.
#' Result: returns a correlation matrix (coefficients); for partial correlation
#' Returns the controlled matrix.
#' Note: the correlation coefficient only measures association strength and
#' Direction, not causation; near 0 does not mean no relation (could be
#' Non-linear). Spearman suits ordinal data or clearly skewed data. Partial
#' Correlation helps identify "spurious correlation" caused by a common confounder.
#' @param .data A `data.frame`.
#' @param vars Variable names; `NULL` = all numeric columns.
#' @param method `"pearson"` (default) or `"spearman"`.
#' @param control Optional control variable names for partial correlation.
#' @param conf Confidence level for the correlation interval; default `0.95`.
#' @param p_adjust Multiplicity method for the family of pairwise tests;
#'   default `"holm"`.
#' @return A long `data.frame` with columns `var1`, `var2`, `r`, `conf.low`,
#'   `conf.high`, `n`, `method`, `p.value`, `p.adjusted`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' mc_correlate(d, vars = c("age", "bmi", "sbp_baseline", "ldl", "glucose"))
#' mc_correlate(d, vars = c("sbp_week12", "dbp_week12"),
#'              method = "spearman", control = "age")
#' @export
mc_correlate <- function(.data, vars = NULL, method = c("pearson", "spearman"),
 control = NULL, conf = 0.95,
 p_adjust = c("holm", "bonferroni", "BH", "BY", "none")) {
 method <- match.arg(method)
 p_adjust <- match.arg(p_adjust)
 if (is.null(vars)) vars <- names(.data)[sapply(.data, is.numeric)]
 mat <-.data[vars]
 if (!is.null(control)) {
 if (method != "pearson") warning("Partial correlation is based on residuals; method forced to pearson.")
 method <- "pearson"
 for (v in vars) {
 f <- stats::as.formula(paste(v, "~", paste(control, collapse = "+")))
 fit <- stats::lm(f, data =.data)
 mat[[v]] <- stats::residuals(fit)
 }
 }
 n <- length(vars)
 if (n < 2L) stop("`vars` must contain at least two numeric columns.")
 pairs <- utils::combn(n, 2)
 rows <- list()
 for (k in seq_len(ncol(pairs))) {
 i <- pairs[1, k]; j <- pairs[2, k]
 a <- mat[[i]]; b <- mat[[j]]
 ok <- !is.na(a) & !is.na(b)
 t <- suppressWarnings(stats::cor.test(a[ok], b[ok], method = method,
 conf.level = conf))
 ci <- tryCatch(as.numeric(t$conf.int), error = function(e) c(NA_real_, NA_real_))
 rows[[length(rows) + 1]] <- data.frame(
 var1 = vars[i], var2 = vars[j],
 r = as.numeric(t$estimate),
 conf.low = ci[1], conf.high = ci[2],
 n = sum(ok), method = method,
 p.value = as.numeric(t$p.value), stringsAsFactors = FALSE)
 }
 out <- do.call(rbind, rows)
 ## A correlation matrix implies many simultaneous tests; report adjusted p.
 out$p.adjusted <- as.numeric(stats::p.adjust(out$p.value, method = p_adjust))
 rownames(out) <- NULL
 out
}

#' Baseline characteristics table
#'
#' Numeric variables are summarised as
#' Mean +/- SD; categorical variables as counts and percentages; can be
#' Stratified by a grouping variable.
#'
#' @details
#' Purpose: generate the common "Table 1" baseline characteristics table, one row
#' Per variable, shown by exposure/group column, ready for the manuscript.
#' Methods: numeric variables summarised as "mean +/- SD" (when approximately
#' Normal); categorical variables as counts and percentages.
#' Arguments: `.data` is a data frame; `vars` are the variable names to include;
#' `group` is an optional grouping column (e.g. treatment), summarising per group;
#' `digits` controls decimal places, default 2.
#' Result: returns a tidy data frame, directly usable with knitr::kable() or
#' print() to make a table.
#' Note: p-values for group comparisons are usually not in this table (compute
#' Separately with mc_ttest/mc_test_rank/mc_crosstab_test). Too many levels for a
#' Categorical variable widens the table; consider merging rare categories first.
#' This is a descriptive table, not a substitute for statistical inference.
#' @param .data A `data.frame`.
#' @param vars Variable names to include.
#' @param group Optional grouping variable name.
#' @param digits Number of decimals; default `2`.
#' @param p_value Add group-comparison p values (Holm-adjusted across `vars`);
#'   default `TRUE`.
#' @param smd Add the standardized mean difference for numeric variables;
#'   default `TRUE`.
#' @param p_adjust Multiplicity method for the p-value column; default `"holm"`.
#' @return A `data.frame` with columns `variable`, `level`, `stratum`, `value`,
#'   and, when `group` is given, `test`, `p.value`, `p.adjusted`, `effect`,
#'   `effect_size` and `smd`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' mc_baseline_table(d,
#'   vars = c("age", "bmi", "sex", "smoke", "sbp_baseline"),
#'   group = "treatment")
#' @export
mc_baseline_table <- function(.data, vars, group = NULL, digits = 2,
 p_value = TRUE, smd = TRUE,
 p_adjust = c("holm", "bonferroni", "BH", "BY", "none")) {
 p_adjust <- match.arg(p_adjust)
 cmp <- NULL
 if (!is.null(group) && (p_value || smd)) {
 cmp <- tryCatch(mc_compare_groups(.data, group = group, vars = vars,
 p_adjust = p_adjust),
 error = function(e) NULL)
 }
 smd_of <- function(x, g) {
 lv <- levels(droplevels(factor(g)))
 if (length(lv) != 2L || !is.numeric(x)) return(NA_real_)
 a <- x[g == lv[1] & !is.na(x)]; b <- x[g == lv[2] & !is.na(x)]
 den <- sqrt((stats::var(a) + stats::var(b)) / 2)
 if (!is.finite(den) || den == 0) return(NA_real_)
 (mean(a) - mean(b)) / den
 }
 rows <- list()
 for (v in vars) {
 x <-.data[[v]]
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
 out <- do.call(rbind, rows)
 if (!is.null(cmp)) {
 out$test <- cmp$test[match(out$variable, cmp$variable)]
 out$p.value <- cmp$p.value[match(out$variable, cmp$variable)]
 out$p.adjusted <- cmp$p.adjusted[match(out$variable, cmp$variable)]
 out$effect <- cmp$effect[match(out$variable, cmp$variable)]
 out$effect_size <- cmp$estimate[match(out$variable, cmp$variable)]
 }
 if (smd && !is.null(group)) {
 out$smd <- vapply(seq_len(nrow(out)), function(i)
 smd_of(.data[[as.character(out$variable[i])]], .data[[group]]),
 numeric(1))
 }
 out
}
