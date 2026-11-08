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
# Module : quantitative description & hypothesis testing module (desc-quant)
# This file contains :.desc_one (internal) / mc_describe_numeric / mc_test_normality /
# mc_test_equal_variance / mc_ttest / mc_anova_oneway / mc_test_rank
# where.desc_one self-implements skewness, kurtosis and
# Shapiro-Wilk normality statistics, and mc_test_equal_variance's Levene
# method is an original lm-based implementation.
# =============================================================================

# Internal: compute a set of descriptive statistics for one variable
# (Fisher-Pearson skewness / kurtosis)
.desc_one <- function(x, conf = 0.95) {
 n_missing <- sum(is.na(x))
 x <- x[!is.na(x)]
 n <- length(x)
 if (n < 3) {
 return(setNames(rep(NA_real_, 18),
 c("n", "n_missing", "mean", "median", "sd", "var", "min", "max",
 "range", "IQR", "SE.mean", "CI.mean.low", "CI.mean.high", "CV",
 "skewness", "kurtosis", "shapiro.W", "shapiro.p")))
 }
 m <- mean(x); s <- sd(x)
 se <- s / sqrt(n)
 qt <- stats::qt(1 - (1 - conf) / 2, df = n - 1)
 z <- (x - m) / s
 g1 <- (n / ((n - 1) * (n - 2))) * sum(z^3)
 g2 <- ((n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3))) * sum(z^4) -
 3 * (n - 1)^2 / ((n - 2) * (n - 3))
 sh <- tryCatch(stats::shapiro.test(x), error = function(e) NULL)
 c(
 n = n, n_missing = n_missing, mean = m, median = stats::median(x),
 sd = s, var = stats::var(x), min = min(x), max = max(x),
 range = max(x) - min(x), IQR = stats::IQR(x),
 SE.mean = se, CI.mean.low = m - qt * se, CI.mean.high = m + qt * se,
 CV = if (m != 0) s / abs(m) else NA_real_,
 skewness = g1, kurtosis = g2,
 shapiro.W = if (is.null(sh)) NA_real_ else as.numeric(sh$statistic),
 shapiro.p = if (is.null(sh)) NA_real_ else as.numeric(sh$p.value)
 )
}

#' One-shot quantitative description
#'
#' Produces about 20 descriptive statistics
#' In one call (n, n missing, mean, median, sd, variance, min, max, range, IQR,
#' SE, 95% CI, CV, skewness, kurtosis, Shapiro-Wilk normality). Supports grouped
#' Computation.
#'
#' @details
#' Purpose: one-call descriptive statistics for quantitative variables, covering
#' About 20 common metrics; the main function for baseline tables and pre-modeling
#' Exploration.
#' Arguments: `.data` is a data frame; `vars` selects the numeric columns to
#' Describe (all numeric columns by default); `group` optionally specifies a
#' Grouping column for per-group output; `conf` is the confidence level for the
#' Mean, default 0.95.
#' Result: returns a data frame with sample size, mean, sd, median, IQR, min, max,
#' Skewness, kurtosis, normality-test statistics, and the mean confidence interval.
#' Note: skewness/kurtosis help judge distribution shape (skewness far from 0 or
#' Kurtosis much greater than 3 suggests non-normality); if planning t-tests/ANOVA,
#' Check normality and variance homogeneity first. When grouping, group sizes may
#' Differ, so mind the fairness of comparisons.
#' @param .data A `data.frame`.
#' @param vars Numeric variable names; `NULL` = all numeric columns.
#' @param group Optional grouping variable name.
#' @param conf Confidence level; default `0.95`.
#' @return A wide `data.frame`, one row per variable (one row per group level
#'   when grouped), with columns `variable`, `group` (when grouped), the 18
#'   statistic columns and a plain-language `distribution` verdict.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' mc_describe_numeric(d, vars = c("age", "bmi", "sbp_baseline"))
#' mc_describe_numeric(d, vars = c("crp", "triglycerides"), group = "treatment")
#' @export
mc_describe_numeric <- function(.data, vars = NULL, group = NULL, conf = 0.95) {
 if (!is.data.frame(.data)) stop("`.data` must be a data.frame.")
 if (is.null(vars)) vars <- names(.data)[sapply(.data, is.numeric)]
 stopifnot(all(vars %in% names(.data)))

 rows <- list()
 if (is.null(group)) {
 for (v in vars) rows[[length(rows) + 1]] <- c(variable = v,.desc_one(.data[[v]], conf))
 } else {
 gvec <-.data[[group]]
 for (lv in unique(gvec)) {
 sub <-.data[gvec == lv, vars, drop = FALSE]
 for (v in vars) {
 one <-.desc_one(sub[[v]], conf)
 row <- c(setNames(as.character(lv), group), variable = v, one)
 rows[[length(rows) + 1]] <- row
 }
 }
 }
 mat <- do.call(rbind, rows)
 df <- as.data.frame(mat, stringsAsFactors = FALSE)
 for (cn in names(df)) if (!cn %in% c("variable", group)) {
 df[[cn]] <- as.numeric(df[[cn]])
 }
 if (!is.null(group)) df <- df[, c(group, setdiff(names(df), group))]
 if ("shapiro.p" %in% names(df)) {
 skewed <- !is.na(df$skewness) & abs(df$skewness) > 1
 heavy <- !is.na(df$kurtosis) & df$kurtosis > 3
 sh <- !is.na(df$shapiro.p) & df$shapiro.p < (1 - conf) / 2 * 2
 df$distribution <- ifelse(is.na(df$shapiro.p), "not tested",
 ifelse(skewed | heavy | sh,
 "non-normal: inspect a Q-Q plot, consider a rank test",
 "approximately normal"))
 }
 df
}

#' Normality test
#'
#' Shapiro-Wilk (default, robust for n <= 5000) or
#' Kolmogorov-Smirnov. Returns the test name, statistic and p-value.
#'
#' @details
#' Purpose: test whether a quantitative variable follows a normal distribution;
#' A prerequisite for choosing parametric vs non-parametric methods.
#' Methods: `shapiro` (default) is the Shapiro-Wilk test, sensitive and powerful
#' For small-to-moderate samples (typically n < 5000); `ks` is the
#' Kolmogorov-Smirnov test, by default compared against the standard normal,
#' Sensitive to location/scale.
#' Arguments: `x` is the numeric vector to test; `method` selects one of the two.
#' Result: returns the test statistic, p-value, etc.
#' Note: with large samples normality tests easily become "significant" (tiny p)
#' Even when deviation is practically harmless, so do not rely on p alone - combine
#' With histograms, Q-Q plots and skewness/kurtosis. The Shapiro test may be
#' Unreliable or very slow for n > 5000.
#' @param x Numeric vector.
#' @param method `"shapiro"` (default) or `"ks"` (test method).
#' @param alpha Significance level for the verdict; default `0.05`.
#' @param ... Extra arguments passed to `stats::ks.test()`.
#' @return A `data.frame` with columns `test`, `statistic`, `p.value`, `n`,
#'   `skewness`, `kurtosis` and `verdict`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' mc_test_normality(d$crp)             # simulated skewed marker
#' mc_test_normality(d$age)
#' @export
mc_test_normality <- function(x, method = c("shapiro", "ks"), alpha = 0.05,...) {
 method <- match.arg(method)
 x <- x[!is.na(x)]
 if (method == "shapiro") {
 t <- stats::shapiro.test(x)
 test <- "Shapiro-Wilk"
 stat <- as.numeric(t$statistic)
 pval <- as.numeric(t$p.value)
 } else {
 t <- stats::ks.test(x, "pnorm", mean(x), stats::sd(x),...)
 test <- "Kolmogorov-Smirnov"
 stat <- as.numeric(t$statistic)
 pval <- as.numeric(t$p.value)
 }
 d <-.desc_one(x)
 sk <- unname(d["skewness"]); ku <- unname(d["kurtosis"])
 verdict <- if (!is.finite(pval)) "not tested" else
 if (pval >= alpha) "no evidence against normality" else
 "departure from normality; consider a rank test or transformation"
 data.frame(test = test, statistic = stat, p.value = pval,
 n = length(x), skewness = sk, kurtosis = ku,
 verdict = verdict, stringsAsFactors = FALSE)
}

#' Variance homogeneity test
#'
#'
#' Bartlett, Fligner-Killeen, and Levene (mean-deviation method). Returns the
#' Test name, statistic and p-value.
#'
#' @details
#' Purpose: test whether several groups have equal population variances (variance
#' Homogeneity), a prerequisite before choosing ANOVA / t-tests (which usually
#' Assume equal variance).
#' Methods: `bartlett` (default) is most sensitive for normal data; `fligner`
#' (Fligner-Killeen) is robust to non-normality; `levene` (Levene) is based on
#' Deviations and robust to heavy-tailed distributions, the most commonly used.
#' Arguments: `x` is a numeric vector; `g` is the grouping factor/vector (same
#' Length as `x`); `method` selects one of the three.
#' Result: returns the test statistic, degrees of freedom and p-value.
#' Note: if variances are unequal (small p), use Welch correction (e.g. mc_ttest
#' Defaults to Welch for two samples) or a non-parametric method instead of the
#' Standard t-test/ANOVA. Bartlett relies strongly on normality; for non-normal
#' Data prefer Levene or Fligner.
#' @param x Numeric vector.
#' @param g A grouping factor.
#' @param method `"bartlett"` (default) / `"fligner"` / `"levene"` (test method).
#' @return A `data.frame` with columns `test`, `statistic`, `df`, `p.value`,
#'   `n`, `n_groups`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' mc_test_equal_variance(d$sbp_week12, d$treatment)
#' mc_test_equal_variance(d$crp, d$site, method = "levene")
#' @export
mc_test_equal_variance <- function(x, g, method = c("bartlett", "fligner", "levene")) {
 method <- match.arg(method)
 g <- as.factor(g)
 if (method == "bartlett") {
 t <- stats::bartlett.test(x, g)
 nm <- "Bartlett"
 } else if (method == "fligner") {
 t <- stats::fligner.test(x, g)
 nm <- "Fligner-Killeen"
 } else {
 means <- tapply(x, g, mean, na.rm = TRUE)
 d <- abs(x - means[as.character(g)])
 a <- stats::anova(stats::lm(d ~ g))
 return(data.frame(test = "Levene",
 statistic = as.numeric(a["g", "F value"]),
 df = as.numeric(a["g", "Df"]),
 p.value = as.numeric(a["g", "Pr(>F)"]),
 n = sum(!is.na(x)), n_groups = nlevels(g),
 stringsAsFactors = FALSE))
 }
 data.frame(test = nm, statistic = as.numeric(t$statistic),
 df = as.numeric(t$parameter), p.value = as.numeric(t$p.value),
 n = sum(!is.na(x)), n_groups = nlevels(g),
 stringsAsFactors = FALSE)
}

#' T-test (unified interface)
#'
#' One-sample, independent two-sample, and paired
#' T-tests through a single entry point.
#'
#' @details
#' Purpose: a unified entry point for t-tests covering three common cases:
#' One-sample (vs a known constant), independent two-sample, and paired.
#' Methods: two independent samples use Welch's corrected t-test by default (no
#' Equal-variance assumption), closer to practice; can switch to Student's t.
#' Arguments: `x` is the first group or one-sample vector; `y` is the second
#' Vector (for independent/paired; left NULL for one-sample); `paired = TRUE`
#' Runs a paired t-test; `mu` is the reference value for the one-sample test,
#' Default 0; `conf`/`alternative` control confidence level and the alternative
#' Direction.
#' Result: returns the mean difference (or sample mean), standard error, t
#' Statistic, degrees of freedom, p-value and confidence interval.
#' Note: paired tests require x, y to be paired (e.g. same person measured
#' Before/after); using an independent test when groups are not independent
#' Understates error. For clearly skewed or small samples, consider mc_test_rank()'s
#' Wilcoxon test.
#' @param x Numeric vector.
#' @param y A second vector (for independent/paired); `NULL` for one-sample.
#' @param paired Logical; paired? default `FALSE`.
#' @param mu Hypothesised mean for one-sample; default `0`.
#' @param conf Confidence level; default `0.95`.
#' @param alternative `"two.sided"` (default) / `"less"` / `"greater"`.
#' @return A `data.frame` with the original columns `method`, `estimate`,
#'   `conf.low`, `conf.high`, `statistic`, `p.value`, plus `n1`, `n2`,
#'   `mean1`, `mean2`, `sd1`, `sd2`, `effect`, `cohens_d`, `hedges_g`,
#'   `effect_ci.low` and `effect_ci.high`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' a <- d$sbp_week12[d$treatment == "Active"]
#' p <- d$sbp_week12[d$treatment == "Placebo"]
#' mc_ttest(a, p)                       # Hedges' g and its CI
#' mc_ttest(a, mu = 135)                # one-sample test
#' @export
mc_ttest <- function(x, y = NULL, paired = FALSE, mu = 0, conf = 0.95,
 alternative = "two.sided") {
 if (is.null(y) && !paired) {
 tt <- stats::t.test(x, mu = mu, conf.level = conf, alternative = alternative)
 } else {
 tt <- stats::t.test(x, y, paired = paired, conf.level = conf, alternative = alternative)
 }
 est <- tt$estimate
 if (length(est) == 2) est <- est[1] - est[2]
 zc <- stats::qnorm(1 - (1 - conf) / 2)
 ## Effect size with an approximate CI: a bare p value is not reportable.
 if (is.null(y) && !paired) {
 xv <- x[!is.na(x)]; n1 <- length(xv); n2 <- NA_integer_
 d <- (mean(xv) - mu) / stats::sd(xv)
 g <- NA_real_; eff <- "Cohen d (one-sample)"
 se <- sqrt(1 / n1 + d^2 / (2 * (n1 - 1)))
 eci <- d + c(-1, 1) * zc * se
 m1 <- mean(xv); m2 <- NA_real_; s1 <- stats::sd(xv); s2 <- NA_real_
 } else if (paired) {
 dv <- (x - y)[!is.na(x - y)]; n1 <- length(dv); n2 <- n1
 d <- mean(dv) / stats::sd(dv)
 g <- NA_real_; eff <- "Cohen dz (paired)"
 se <- sqrt(1 / n1 + d^2 / (2 * (n1 - 1)))
 eci <- d + c(-1, 1) * zc * se
 m1 <- mean(x, na.rm = TRUE); m2 <- mean(y, na.rm = TRUE)
 s1 <- stats::sd(x, na.rm = TRUE); s2 <- stats::sd(y, na.rm = TRUE)
 } else {
 a <- x[!is.na(x)]; b <- y[!is.na(y)]
 n1 <- length(a); n2 <- length(b)
 sp <- sqrt(((n1 - 1) * stats::var(a) + (n2 - 1) * stats::var(b)) /
 (n1 + n2 - 2))
 d <- (mean(a) - mean(b)) / sp
 j <-.mc_hedges_j(n1, n2)
 g <- d * j; eff <- "Hedges g (two-sample)"
 eci <-.mc_d_ci(d, n1, n2, conf) * j
 m1 <- mean(a); m2 <- mean(b)
 s1 <- stats::sd(a); s2 <- stats::sd(b)
 }
 data.frame(method = tt$method, estimate = as.numeric(est),
 conf.low = as.numeric(tt$conf.int[1]),
 conf.high = as.numeric(tt$conf.int[2]),
 statistic = as.numeric(tt$statistic),
 p.value = as.numeric(tt$p.value),
 n1 = n1, n2 = n2, mean1 = m1, mean2 = m2, sd1 = s1, sd2 = s2,
 effect = eff, cohens_d = d, hedges_g = g,
 effect_ci.low = eci[1], effect_ci.high = eci[2],
 stringsAsFactors = FALSE)
}

#' One-way ANOVA
#'
#' Built on `stats::aov`, returns the
#' Degrees of freedom, sum of squares, mean square, F statistic and p-value of
#' Each term.
#'
#' @details
#' Purpose: one-way ANOVA, comparing whether the means of three or more groups
#' Are equal overall.
#' Methods: decomposes variation via a linear model, comparing between-group and
#' Within-group mean squares to obtain the F statistic.
#' Arguments: `formula` like y ~ g (response ~ grouping factor); `data` is the
#' Data frame containing both columns.
#' Result: returns the degrees of freedom, sum of squares, mean square, F value
#' And p-value of each term.
#' Note: ANOVA only tells you "at least two groups differ"; which two require
#' Post-hoc multiple comparisons (e.g. Tukey HSD). Its assumptions are normality
#' And equal variance per group - verify first with mc_test_normality() and mc_test_equal_variance().
#' If variances are unequal or distributions skewed, consider Kruskal-Wallis
#' (see mc_test_rank).
#' @param formula A formula like `y ~ group`.
#' @param data A `data.frame` containing the variables.
#' @return A `data.frame` with columns `term`, `df`, `sumsq`, `meansq`,
#'   `statistic`, `p.value`, `eta_sq`, `omega_sq`. Use [mc_posthoc_tukey()]
#'   for the pairwise contrasts that a significant omnibus F requires.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' mc_anova_oneway(sbp_week12 ~ race, data = d)
#' mc_anova_oneway(crp ~ site, data = d)
#' @export
mc_anova_oneway <- function(formula, data) {
 a <- stats::aov(formula, data = data)
 s <- summary(a)[[1]]
 eta <- rep(NA_real_, nrow(s)); omega <- rep(NA_real_, nrow(s))
 ssb <- s[1, "Sum Sq"]; ssw <- s[nrow(s), "Sum Sq"]
 dfb <- s[1, "Df"]; msw <- s[nrow(s), "Mean Sq"]
 ## eta-squared and omega-squared quantify how much variance the factor
 ## explains; omega-squared is the less biased of the two.
 if (is.finite(ssb) && is.finite(ssw) && ssb + ssw > 0) {
 eta[1] <- ssb / (ssb + ssw)
 omega[1] <- (ssb - dfb * msw) / (ssb + ssw + msw)
 }
 out <- data.frame(term = rownames(s), df = s[, "Df"], sumsq = s[, "Sum Sq"],
 meansq = s[, "Mean Sq"], statistic = s[, "F value"],
 p.value = s[, "Pr(>F)"], eta_sq = eta, omega_sq = omega,
 stringsAsFactors = FALSE)
 out
}

#' Non-parametric group difference test
#'
#' Wilcoxon
#' Rank-sum (two groups) or Kruskal-Wallis (multiple groups). Returns the test
#' Name, statistic and p-value.
#'
#' @details
#' Purpose: when data violate normality/equal-variance assumptions, or are merely
#' Ordinal, use non-parametric methods to compare groups.
#' Methods: `wilcoxon` is the Wilcoxon rank-sum test (Mann-Whitney U) for two
#' Independent samples; `kruskal` is the Kruskal-Wallis test for three or more
#' Groups. They compare distribution location (ranks), not means.
#' Arguments: `x` is a numeric vector; `g` is the grouping vector (wilcoxon for
#' Two groups, kruskal for multiple); `method` selects the test.
#' Result: returns the test statistic and p-value (Kruskal also includes df).
#' Note: non-parametric tests usually have lower power than their parametric
#' Counterparts (when data are truly normal) but are more robust; they test
#' "whether distributions are identical / locations shifted", so conclusions
#' Should be phrased as "median/distribution difference" not "mean difference".
#' For paired two groups use the Wilcoxon signed-rank test (not covered here;
#' Use stats::wilcox.test(paired = TRUE)).
#' @param x Numeric vector.
#' @param g Grouping vector/factor.
#' @param method `"wilcoxon"` (default) or `"kruskal"` (test method).
#' @param conf Confidence level for the rank-biserial interval; default `0.95`.
#' @return A `data.frame` with columns `test`, `statistic`, `df`, `p.value`,
#'   `n`, `effect`, `effect_size`, `effect_ci.low`, `effect_ci.high`, `note`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' mc_test_rank(d$crp, d$treatment)                       # two groups, Wilcoxon
#' mc_test_rank(d$crp, d$site, method = "kruskal")       # several groups
#' @export
mc_test_rank <- function(x, g, method = c("wilcoxon", "kruskal"),
 conf = 0.95) {
 method <- match.arg(method)
 g <- droplevels(as.factor(g))
 keep <- !is.na(x) & !is.na(g)
 x <- x[keep]; g <- g[keep]
 if (method == "wilcoxon") {
 t <- stats::wilcox.test(x ~ g, exact = FALSE, conf.level = conf)
 nm <- "Wilcoxon rank-sum"
 lv <- levels(g)
 a <- x[g == lv[1]]; b <- x[g == lv[2]]
 n1 <- length(a); n2 <- length(b)
 u <- sum(rank(c(a, b))[seq_len(n1)]) - n1 * (n1 + 1) / 2
 r <- 2 * u / (n1 * n2) - 1
 ci <-.mc_r_ci(r, n1 + n2, conf)
 eff <- "rank-biserial r"
 val <- r; lo <- ci[1]; hi <- ci[2]
 dfv <- NA_real_
 note <- paste0("n=", n1, "/", n2, "; Fisher-z CI")
 } else {
 t <- stats::kruskal.test(x ~ g)
 nm <- "Kruskal-Wallis"
 k <- nlevels(g); n <- length(x)
 h <- as.numeric(t$statistic)
 val <- (h - k + 1) / (n - k)
 eff <- "epsilon-squared"
 lo <- NA_real_; hi <- NA_real_
 dfv <- if (is.null(t$parameter)) NA_real_ else as.numeric(t$parameter)
 note <- paste0("k=", k, " groups; n=", n,
 "; pairwise contrasts need multiplicity adjustment")
 }
 data.frame(test = nm, statistic = as.numeric(t$statistic),
 df = dfv, p.value = as.numeric(t$p.value),
 n = length(x), effect = eff, effect_size = val,
 effect_ci.low = lo, effect_ci.high = hi, note = note,
 stringsAsFactors = FALSE)
}
