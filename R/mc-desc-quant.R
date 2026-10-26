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

#' One-shot quantitative description (flagship function)
#'
#' Corresponds to Chapter 3 of Part 2: produces about 20 descriptive
#' statistics in one call (n, n missing, mean, median, sd, variance, min,
#' max, range, IQR, SE, 95% CI, CV, skewness, kurtosis, Shapiro-Wilk
#' normality). Supports grouped computation.
#'
#' @param .data A data.frame.
#' @param vars Numeric variable names; NULL means all numeric columns.
#' @param group Optional grouping variable name.
#' @param conf Confidence level. Default 0.95.
#' @return A wide data.frame, one row per variable (one row per group if grouped).
#' @examples
#' mc_desc_quant(MASS::birthwt, vars = c("age", "lwt", "bwt"))
#' mc_desc_quant(MASS::birthwt, vars = c("age", "lwt"), group = "smoke")
#' @export
mc_desc_quant <- function(.data, vars = NULL, group = NULL, conf = 0.95) {
  if (!is.data.frame(.data)) stop("`.data` must be a data.frame.")
  if (is.null(vars)) vars <- names(.data)[sapply(.data, is.numeric)]
  stopifnot(all(vars %in% names(.data)))

  rows <- list()
  if (is.null(group)) {
    for (v in vars) rows[[length(rows) + 1]] <- c(variable = v, .desc_one(.data[[v]], conf))
  } else {
    gvec <- .data[[group]]
    for (lv in unique(gvec)) {
      sub <- .data[gvec == lv, vars, drop = FALSE]
      for (v in vars) {
        one <- .desc_one(sub[[v]], conf)
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
#' Corresponds to Chapter 3: Shapiro-Wilk (default, robust for n <= 5000)
#' or Kolmogorov-Smirnov.
#'
#' @param x A numeric vector.
#' @param method "shapiro" (default) or "ks".
#' @param ... Extra arguments passed to \code{ks.test}.
#' @return A data.frame: test / statistic / p.value.
#' @examples
#' mc_norm_test(MASS::birthwt$age)
#' @export
mc_norm_test <- function(x, method = c("shapiro", "ks"), ...) {
  method <- match.arg(method)
  x <- x[!is.na(x)]
  if (method == "shapiro") {
    t <- stats::shapiro.test(x)
    data.frame(test = "Shapiro-Wilk", statistic = as.numeric(t$statistic),
               p.value = as.numeric(t$p.value), stringsAsFactors = FALSE)
  } else {
    t <- stats::ks.test(x, "pnorm", mean(x), stats::sd(x), ...)
    data.frame(test = "Kolmogorov-Smirnov", statistic = as.numeric(t$statistic),
               p.value = as.numeric(t$p.value), stringsAsFactors = FALSE)
  }
}

#' Variance homogeneity test
#'
#' Corresponds to Chapter 3 ("three methods for variance homogeneity"):
#' Bartlett, Fligner-Killeen, and Levene (mean-deviation method).
#'
#' @param x A numeric vector.
#' @param g A grouping factor.
#' @param method "bartlett" (default), "fligner" or "levene".
#' @return A data.frame: test / statistic / p.value.
#' @examples
#' mc_var_test(MASS::birthwt$bwt, MASS::birthwt$smoke, method = "bartlett")
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

#' t-test (unified interface)
#'
#' Corresponds to Chapter 3: one-sample, independent two-sample, and
#' paired t-tests through a single entry point.
#'
#' @param x A numeric vector.
#' @param y A second vector (for independent/paired); NULL for one-sample.
#' @param paired Logical; paired? Default FALSE.
#' @param mu Hypothesised mean for one-sample test. Default 0.
#' @param conf Confidence level. Default 0.95.
#' @param alternative "two.sided" (default), "less" or "greater".
#' @return A data.frame: method / estimate / conf.low / conf.high / statistic / p.value
#' @examples
#' set.seed(1); mc_ttest(rnorm(30, mean = 1), mu = 0)
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
#' Corresponds to Chapter 3 ("ANOVA").
#'
#' @param formula A formula like \code{y ~ group}.
#' @param data A data.frame containing the variables.
#' @return A data.frame: term / df / sumsq / meansq / statistic / p.value
#' @examples
#' mc_anova(bwt ~ factor(race), data = MASS::birthwt)
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
#' Corresponds to Chapter 3 ("non-parametric tests between groups"):
#' Wilcoxon rank-sum (two groups) or Kruskal-Wallis (multiple groups).
#'
#' @param x A numeric vector.
#' @param g Grouping vector/factor.
#' @param method "wilcoxon" (default) or "kruskal".
#' @return A data.frame: test / statistic / p.value
#' @examples
#' mc_nonparam(MASS::birthwt$bwt, MASS::birthwt$smoke, method = "wilcoxon")
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
