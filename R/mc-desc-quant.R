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
# This file contains :.desc_one (internal) / mc_desc_quant / mc_norm_test /
# mc_var_test / mc_ttest / mc_anova / mc_nonparam
# where.desc_one self-implements skewness, kurtosis and
# Shapiro-Wilk normality statistics, and mc_var_test's Levene
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
#' @param.data A `data.frame`.
#' @param vars Numeric variable names; `NULL` = all numeric columns.
#' @param group Optional grouping variable name.
#' @param conf Confidence level; default `0.95`.
#' @return A wide `data.frame`, one row per variable (one row per group level
#' When grouped), with columns `variable`, `group` (when grouped), and the 18
#' Statistic columns.
#' @examples
#' # low birth weight study (birthwt, MASS package)
#' birthwt <- MASS::birthwt
#' cat("==== overall: birth weight / age / pregnancy weight ====\n")
#' print(mc_desc_quant(birthwt, vars = c("bwt","age","lwt")), row.names = FALSE)
#' cat("\n==== grouped by mother's smoking status ====\n")
#' print(mc_desc_quant(birthwt, vars = c("bwt","age"), group = "smoke"), row.names = FALSE)
#' @export
mc_desc_quant <- function(.data, vars = NULL, group = NULL, conf = 0.95) {
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
#' @param... Extra arguments passed to `stats::ks.test()`.
#' @return A `data.frame` with columns `test`, `statistic`, `p.value`.
#' @examples
#' # test whether quantitative variables in birthwt are normal
#' birthwt <- MASS::birthwt
#' mc_norm_test(birthwt$bwt) # birth weight
#' mc_norm_test(birthwt$age) # mother's age
#' @export
mc_norm_test <- function(x, method = c("shapiro", "ks"),...) {
 method <- match.arg(method)
 x <- x[!is.na(x)]
 if (method == "shapiro") {
 t <- stats::shapiro.test(x)
 data.frame(test = "Shapiro-Wilk", statistic = as.numeric(t$statistic),
 p.value = as.numeric(t$p.value), stringsAsFactors = FALSE)
 } else {
 t <- stats::ks.test(x, "pnorm", mean(x), stats::sd(x),...)
 data.frame(test = "Kolmogorov-Smirnov", statistic = as.numeric(t$statistic),
 p.value = as.numeric(t$p.value), stringsAsFactors = FALSE)
 }
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
#' @return A `data.frame` with columns `test`, `statistic`, `p.value`.
#' @examples
#' # test variance homogeneity of birth weight across races (race)
#' birthwt <- MASS::birthwt
#' mc_var_test(birthwt$bwt, birthwt$race, method = "bartlett")
#' mc_var_test(birthwt$bwt, birthwt$race, method = "levene")
#' @export
mc_var_test <- function(x, g, method = c("bartlett", "fligner", "levene")) {
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
 p.value = as.numeric(a["g", "Pr(>F)"]),
 stringsAsFactors = FALSE))
 }
 data.frame(test = nm, statistic = as.numeric(t$statistic),
 p.value = as.numeric(t$p.value), stringsAsFactors = FALSE)
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
#' Understates error. For clearly skewed or small samples, consider mc_nonparam()'s
#' Wilcoxon test.
#' @param x Numeric vector.
#' @param y A second vector (for independent/paired); `NULL` for one-sample.
#' @param paired Logical; paired? default `FALSE`.
#' @param mu Hypothesised mean for one-sample; default `0`.
#' @param conf Confidence level; default `0.95`.
#' @param alternative `"two.sided"` (default) / `"less"` / `"greater"`.
#' @return A `data.frame` with columns `method`, `estimate`, `conf.low`,
#' `conf.high`, `statistic`, `p.value`.
#' @examples
#' #
#' birthwt <- MASS::birthwt
#' # one-sample: is the average birth weight equal to 3000 g?
#' mc_ttest(birthwt$bwt, mu = 3000)
#' # independent two-sample (welch correction): smoking vs non-smoking mothers' birth weight
#' mc_ttest(birthwt$bwt[birthwt$smoke == 1], birthwt$bwt[birthwt$smoke == 0])
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
 data.frame(method = tt$method, estimate = as.numeric(est),
 conf.low = as.numeric(tt$conf.int[1]),
 conf.high = as.numeric(tt$conf.int[2]),
 statistic = as.numeric(tt$statistic),
 p.value = as.numeric(tt$p.value), stringsAsFactors = FALSE)
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
#' And equal variance per group - verify first with mc_norm_test() and mc_var_test().
#' If variances are unequal or distributions skewed, consider Kruskal-Wallis
#' (see mc_nonparam).
#' @param formula A formula like `y ~ group`.
#' @param data A `data.frame` containing the variables.
#' @return A `data.frame` with columns `term`, `df`, `sumsq`, `meansq`,
#' `statistic`, `p.value`.
#' @examples
#' # do mothers of different races (race) differ in birth weight? (one-way aNOVA)
#' birthwt <- MASS::birthwt
#' mc_anova(bwt ~ factor(race), data = birthwt)
#' @export
mc_anova <- function(formula, data) {
 a <- stats::aov(formula, data = data)
 s <- summary(a)[[1]]
 data.frame(term = rownames(s), df = s[, "Df"], sumsq = s[, "Sum Sq"],
 meansq = s[, "Mean Sq"], statistic = s[, "F value"],
 p.value = s[, "Pr(>F)"], stringsAsFactors = FALSE)
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
#' @return A `data.frame` with columns `test`, `statistic`, `p.value`.
#' @examples
#' # birth weight by smoking status (two-group wilcoxon) and race (multi-group kruskal-wallis)
#' birthwt <- MASS::birthwt
#' mc_nonparam(birthwt$bwt, birthwt$smoke, method = "wilcoxon")
#' mc_nonparam(birthwt$bwt, factor(birthwt$race), method = "kruskal")
#' @export
mc_nonparam <- function(x, g, method = c("wilcoxon", "kruskal")) {
 method <- match.arg(method)
 g <- as.factor(g)
 if (method == "wilcoxon") {
 t <- stats::wilcox.test(x ~ g)
 nm <- "Wilcoxon rank-sum"
 } else {
 t <- stats::kruskal.test(x ~ g)
 nm <- "Kruskal-Wallis"
 }
 data.frame(test = nm, statistic = as.numeric(t$statistic),
 p.value = as.numeric(t$p.value), stringsAsFactors = FALSE)
}
