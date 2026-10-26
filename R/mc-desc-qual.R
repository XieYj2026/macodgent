#' Contingency table and independence test (with odds ratio)
#'
#' Corresponds to Chapter 4: build a contingency table and run a chi-square
#' independence test; for a 2x2 table, also compute the odds ratio (OR) and its
#' 95% CI from Fisher's exact test.
#'
#' @param .data A data.frame.
#' @param row Row variable name (categorical).
#' @param col Column variable name (categorical).
#' @param correct Continuity correction for chi-square. Default TRUE.
#' @return A list: \code{table} (contingency table), \code{test} (test
#'   data.frame), \code{OR} (only for 2x2 tables).
#' @examples
#' mc_contingency(MASS::birthwt, row = "smoke", col = "low")
#' @export
mc_contingency <- function(.data, row, col, correct = TRUE) {
  tab <- table(.data[[row]], .data[[col]])
  ct <- stats::chisq.test(tab, correct = correct)
  res <- list(
    table = tab,
    test = data.frame(test = "Chi-square",
                      statistic = as.numeric(ct$statistic),
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
#' @param .data A data.frame.
#' @param vars Variable names; NULL means all numeric columns.
#' @param method "pearson" (default) or "spearman".
#' @param control Optional control variable names (for partial correlation).
#' @return A long data.frame: var1 / var2 / r / p.value
#' @examples
#' mc_corr(MASS::birthwt, vars = c("age", "lwt", "bwt"))
#' mc_corr(MASS::birthwt, vars = c("age", "lwt", "bwt"), control = "race")
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

#' Baseline characteristics table (Table One)
#'
#' Corresponds to Chapter 4 ("Table One"). Numeric variables are summarised as
#' mean +/- SD; categorical variables as counts and percentages; can be
#' stratified by a grouping variable.
#'
#' @param .data A data.frame.
#' @param vars Variable names to include.
#' @param group Optional grouping variable name.
#' @param digits Number of decimals. Default 2.
#' @return A data.frame: variable / level / stratum / value
#' @examples
#' mc_tableone(MASS::birthwt, vars = c("age", "lwt", "smoke", "low"), group = "smoke")
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
