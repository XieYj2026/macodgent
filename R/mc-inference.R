# =============================================================================
# Software name : macodgent medical data-science research toolkit
# Copyright holder : Dr.Xie
# Version : 1.0.2
# Language : R
# -----------------------------------------------------------------------------
# Module : inference layer (effect sizes, multiplicity, group comparison)
# Everything here is implemented with base R only. The point of this module is
# to deliver the quantities a medical manuscript actually needs -- effect size,
# its confidence interval, and multiplicity-adjusted p values -- instead of a
# bare p value.
# =============================================================================

# ---- internal helpers -------------------------------------------------------

# Evaluate an expression under a temporary seed, restoring the caller's RNG.
.mc_with_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) {
    old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)), add = TRUE)
  }
  set.seed(seed)
  force(expr)
}

# Hedges' small-sample correction for Cohen's d.
.mc_hedges_j <- function(n1, n2) 1 - 3 / (4 * (n1 + n2) - 9)

# Large-sample SE of Cohen's d (Hedges & Olkin), normal-approximate CI.
.mc_d_ci <- function(d, n1, n2, conf = 0.95) {
  se <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2 - 2)))
  zc <- stats::qnorm(1 - (1 - conf) / 2)
  c(d - zc * se, d + zc * se)
}

# Fisher-z CI for correlation-type effect sizes (r, rank-biserial).
.mc_r_ci <- function(r, n, conf = 0.95) {
  if (!is.finite(r) || abs(r) >= 1 || n <= 3) return(c(NA_real_, NA_real_))
  z <- atanh(r)
  se <- 1 / sqrt(n - 3)
  zc <- stats::qnorm(1 - (1 - conf) / 2)
  tanh(z + c(-1, 1) * zc * se)
}

# Cramer's V from a chi-square statistic.
.mc_cramers_v <- function(chi2, n, r, c) {
  k <- min(r, c) - 1L
  if (k < 1L || n <= 0) return(NA_real_)
  sqrt(chi2 / (n * k))
}

# Expected-count diagnostics for a contingency table.
.mc_table_diagnostics <- function(tab, expected = NULL) {
  if (is.null(expected)) {
    expected <- tryCatch(stats::chisq.test(tab, correct = FALSE)$expected,
                         error = function(e)
                           matrix(NA_real_, nrow(tab), ncol(tab)))
  }
  n_sparse <- sum(expected < 5, na.rm = TRUE)
  min_exp <- suppressWarnings(min(expected, na.rm = TRUE))
  if (!is.finite(min_exp)) min_exp <- NA_real_
  rec <- if (is.finite(min_exp) && n_sparse > 0L)
    "Expected counts below 5: prefer the Fisher exact test." else
    "Chi-square expected counts are adequate."
  data.frame(n = sum(tab), n_cells = length(tab),
             min_expected = min_exp, n_cells_lt5 = n_sparse,
             recommendation = rec, stringsAsFactors = FALSE)
}

# ---- exported: multiplicity -------------------------------------------------

#' Adjust p values for multiple testing
#'
#' Wraps [stats::p.adjust()] and returns a report-ready table that also flags
#' which findings survive the adjustment at a stated alpha.
#'
#' @details
#' Purpose: a p value is only interpretable together with the number of tests
#' performed. Baseline tables and correlation matrices routinely run dozens of
#' tests, so unadjusted p values overstate evidence.
#' Methods: \code{"holm"} (default) is uniformly more powerful than Bonferroni
#' and controls the family-wise error rate under any dependence;
#' \code{"bonferroni"} is the most conservative; \code{"BH"} controls the false
#' discovery rate and is appropriate for exploratory screening; \code{"BY"}
#' is the conservative FDR variant for arbitrary dependence.
#' Arguments: \code{p} is a numeric vector of raw p values; \code{method} is
#' passed to [stats::p.adjust()]; \code{alpha} is the significance level used
#' only to populate the \code{significant} flag.
#' Result: a data frame with the raw p, the adjusted p, the method used and a
#' logical significance flag. \code{NA} p values are preserved.
#' Note: the family of tests must be decided before looking at the data.
#' @param p Numeric vector of raw p values.
#' @param method Adjustment method: `"holm"` (default), `"bonferroni"`,
#'   `"BH"`, `"BY"`, or `"none"`.
#' @param alpha Significance level for the flag; default `0.05`.
#' @return A `data.frame` with columns `p.value`, `p.adjusted`, `method`,
#'   `significant`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' p <- mc_compare_groups(d, group = "treatment",
#'         vars = c("age", "bmi", "crp", "glucose"))$p.value
#' mc_adjust_p(p, method = "holm")
#' @export
mc_adjust_p <- function(p, method = c("holm", "bonferroni", "BH", "BY", "none"),
                        alpha = 0.05) {
  method <- match.arg(method)
  p <- as.numeric(p)
  if (any(!is.na(p) & (p < 0 | p > 1)))
    stop("All p values must lie in the 0-1 range.")
  padj <- stats::p.adjust(p, method = method)
  data.frame(p.value = p, p.adjusted = as.numeric(padj),
             method = method, significant = !is.na(padj) & padj < alpha,
             stringsAsFactors = FALSE)
}

# ---- exported: post-hoc comparisons ----------------------------------------

#' Tukey HSD post-hoc comparisons after one-way ANOVA
#'
#' A significant ANOVA only says that at least two groups differ. This function
#' produces the pairwise comparisons that a manuscript actually reports, with
#' multiplicity-adjusted p values and confidence intervals for each difference.
#'
#' @details
#' Purpose: after [mc_anova_oneway()] rejects the global null, identify which
#' specific pairs differ while controlling the family-wise error rate.
#' Method: Tukey's honestly significant difference (HSD) via
#' [stats::TukeyHSD()], based on the studentised range distribution. Unlike
#' repeated t-tests it controls the family-wise error rate exactly for balanced
#' designs and conservatively otherwise.
#' Arguments: pass either a fitted \code{aov}/\code{lm} object as \code{model},
#' or \code{data} plus \code{formula} for the function to fit the ANOVA first;
#' \code{conf} is the simultaneous confidence level.
#' Result: a tidy data frame with one row per pairwise contrast: the term,
#' the contrast label, the estimated difference, its simultaneous CI and the
#' adjusted p value.
#' Note: Tukey HSD assumes normality and homogeneous variances. When those
#' assumptions fail, use pairwise Wilcoxon tests with \code{mc_adjust_p()}
#' instead. A confidence interval crossing zero agrees with an adjusted
#' p value above 0.05.
#' @param model A fitted `aov` or `lm` object. If `NULL`, supply `data` and
#'   `formula`.
#' @param data A `data.frame` (used when `model` is `NULL`).
#' @param formula A formula such as `y ~ group` (used when `model` is `NULL`).
#' @param conf Simultaneous confidence level; default `0.95`.
#' @return A `data.frame` (class `mc_posthoc`) with columns `term`,
#'   `comparison`, `estimate`, `conf.low`, `conf.high`, `p.adjusted`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' mc_posthoc_tukey(data = d, formula = sbp_week12 ~ race)
#' @export
mc_posthoc_tukey <- function(model = NULL, data = NULL, formula = NULL,
                             conf = 0.95) {
  if (is.null(model)) {
    if (is.null(data) || is.null(formula))
      stop("Provide either a fitted model, or both `data` and `formula`.")
    model <- stats::aov(formula, data = data)
  } else if (inherits(model, "lm") && !inherits(model, "aov")) {
    mf <- stats::model.frame(model)
    model <- stats::aov(stats::formula(model), data = mf)
  }
  if (!inherits(model, "aov")) stop("`model` must be an aov or lm object.")
  th <- stats::TukeyHSD(model, conf.level = conf)
  rows <- list()
  for (nm in names(th)) {
    m <- as.data.frame(th[[nm]])
    rows[[length(rows) + 1]] <- data.frame(
      term = nm, comparison = rownames(m),
      estimate = m[["diff"]], conf.low = m[["lwr"]],
      conf.high = m[["upr"]], p.adjusted = m[["p adj"]],
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("mc_posthoc", class(out))
  out
}

# ---- exported: one-call group comparison -----------------------------------

#' Compare every variable between groups, with effect sizes and adjusted p values
#'
#' The inferential counterpart of a baseline table. For each variable it
#' selects an appropriate test from its type and distribution, then reports the
#' effect size with a confidence interval, the raw p value and a
#' multiplicity-adjusted p value across the whole family of tests.
#'
#' @details
#' Purpose: turn "compare the two arms" into one reproducible call that yields
#' a table which can go straight into a manuscript or supplement. This removes
#' the two most common reporting failures in medical papers: p values with no
#' effect size, and many tests with no multiplicity correction.
#'
#' Method (automatic test selection, recorded in the `note` column):
#' \itemize{
#' \item Numeric variable, two groups: Shapiro-Wilk on each group plus an
#' F test of equal variances. If both pass at `normal_alpha`, use Welch's
#' t-test with Hedges' g (small-sample corrected Cohen's d) and an
#' approximate CI; otherwise use the Wilcoxon rank-sum test with the
#' rank-biserial correlation and a Fisher-z CI.
#' \item Numeric variable, three or more groups: one-way ANOVA with eta-squared
#' (and omega-squared in `note`) when assumptions hold, otherwise
#' Kruskal-Wallis with epsilon-squared.
#' \item Categorical variable: chi-square with Cramer's V when all expected
#' cell counts are at least 5, otherwise a simulated Fisher exact test.
#' }
#'
#' Arguments: `data` is a data frame; `group` names the grouping column (two or
#' more levels); `vars` are the variables to compare (all others by default);
#' `p_adjust` controls multiplicity correction across `vars` (default Holm);
#' `conf` is the confidence level for effect sizes; `alpha` is the significance
#' level for the flag; `normal_alpha` is the threshold used by the assumption
#' checks; `seed` makes the simulated Fisher p value reproducible without
#' changing the caller's random-number state.
#'
#' Result: a data frame with one row per variable (or per contrast), carrying
#' `test`, `statistic`, `effect`, `estimate`, `conf.low`, `conf.high`,
#' `p.value`, `p.adjusted`, `significant` and a `note` explaining the choice.
#'
#' Note: automatic test selection is a convenience, not a substitute for
#' clinical judgement. Shapiro-Wilk becomes over-sensitive in large samples;
#' inspect `note` and the distribution before reporting. Confidence intervals
#' for d and rank-biserial are large-sample approximations. Effects are
#' oriented as level 2 minus level 1 (or the first-listed contrast).
#' Post-hoc pairwise contrasts for more than two groups are available from
#' [mc_posthoc_tukey()].
#'
#' @param data A `data.frame`.
#' @param group Name of the grouping column.
#' @param vars Variables to compare; `NULL` = every other column.
#' @param p_adjust Multiplicity method passed to [mc_adjust_p()].
#' @param conf Confidence level for effect-size intervals; default `0.95`.
#' @param alpha Significance level; default `0.05`.
#' @param normal_alpha Threshold for normality / equal-variance checks.
#' @param seed Seed for the simulated Fisher exact test; `NULL` = no seeding.
#' @return A `data.frame` (class `mc_group_comparison`) with one row per
#'   variable; see Details for the columns.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' res <- mc_compare_groups(d, group = "treatment",
#'         vars = c("age", "bmi", "crp", "smoke", "sex"))
#' res[, c("variable", "test", "effect", "estimate", "p.value", "p.adjusted")]
#' @export
mc_compare_groups <- function(data, group, vars = NULL,
                              p_adjust = c("holm", "bonferroni", "BH", "BY", "none"),
                              conf = 0.95, alpha = 0.05, normal_alpha = 0.05,
                              seed = 1L) {
  p_adjust <- match.arg(p_adjust)
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (length(group) != 1L || !group %in% names(data))
    stop("`group` must name a single column in `data`.")
  g_all <- droplevels(factor(data[[group]]))
  if (nlevels(g_all) < 2L) stop("`group` needs at least two non-missing levels.")
  if (is.null(vars)) vars <- setdiff(names(data), group)
  miss <- setdiff(vars, names(data))
  if (length(miss)) stop("Unknown column(s): ", paste(miss, collapse = ", "))

  new_row <- function(variable, n, group1, group2, test, statistic,
                      effect, estimate, ci, p.value, note) {
    data.frame(variable = variable, n = n, group1 = group1, group2 = group2,
               test = test, statistic = statistic, effect = effect,
               estimate = estimate, conf.low = ci[1], conf.high = ci[2],
               p.value = p.value, note = note, stringsAsFactors = FALSE)
  }

  rows <- list()
  for (v in vars) {
    x <- data[[v]]
    ok <- !is.na(x) & !is.na(g_all)
    xx <- x[ok]
    gg <- droplevels(g_all[ok])
    lv <- levels(gg)
    k <- nlevels(gg)
    if (k < 2L) next

    if (is.numeric(xx)) {
      groups <- lapply(lv, function(l) xx[gg == l])
      ns <- lengths(groups)
      norm_ok <- all(ns >= 3L) && all(ns <= 5000L) &&
        all(vapply(groups, function(a)
          isTRUE(tryCatch(stats::shapiro.test(a)$p.value > normal_alpha,
                          error = function(e) FALSE)), logical(1)))
      var_ok <- if (k == 2L) {
        isTRUE(tryCatch(stats::var.test(groups[[1]], groups[[2]])$p.value >
                          normal_alpha, error = function(e) FALSE))
      } else {
        isTRUE(tryCatch(stats::bartlett.test(xx ~ gg)$p.value > normal_alpha,
                        error = function(e) FALSE))
      }

      if (k == 2L) {
        a <- groups[[1]]; b <- groups[[2]]
        n1 <- length(a); n2 <- length(b); n <- n1 + n2
        if (norm_ok && var_ok) {
          tt <- stats::t.test(a, b, var.equal = FALSE, conf.level = conf)
          sp <- sqrt(((n1 - 1) * stats::var(a) + (n2 - 1) * stats::var(b)) /
                       (n1 + n2 - 2))
          d <- (mean(a) - mean(b)) / sp
          j <- .mc_hedges_j(n1, n2)
          est <- d * j
          ci <- .mc_d_ci(d, n1, n2, conf) * j
          note <- paste0("normal + equal variance; n=", n1, "/", n2)
          rows[[length(rows) + 1]] <- new_row(
            v, n, lv[1], lv[2], "Welch t-test",
            as.numeric(tt$statistic), "Hedges g", est, ci,
            as.numeric(tt$p.value), note)
        } else {
          u <- sum(rank(c(a, b))[seq_len(n1)]) - n1 * (n1 + 1) / 2
          r <- 2 * u / (n1 * n2) - 1
          wt <- stats::wilcox.test(a, b, exact = FALSE, conf.level = conf)
          ci <- .mc_r_ci(r, n, conf)
          note <- paste0(
            "non-normal or unequal variance; rank-biserial via Fisher z",
            if (abs(r) >= 1)
              "; effect at boundary, CI not estimable (check for outcome leakage)"
            else "")
          rows[[length(rows) + 1]] <- new_row(
            v, n, lv[1], lv[2], "Wilcoxon rank-sum", u,
            "rank-biserial r", r, ci, as.numeric(wt$p.value), note)
        }
      } else if (norm_ok && var_ok) {
        fit <- stats::aov(xx ~ gg)
        s <- summary(fit)[[1]]
        ssb <- s[1, "Sum Sq"]; ssw <- s[length(s[, "Sum Sq"]), "Sum Sq"]
        dfb <- s[1, "Df"]
        msw <- s[nrow(s), "Mean Sq"]
        eta2 <- ssb / (ssb + ssw)
        omega2 <- (ssb - dfb * msw) / (ssb + ssw + msw)
        note <- paste0("omega-squared=", signif(omega2, 3),
                       "; use mc_posthoc_tukey() for pairwise contrasts")
        rows[[length(rows) + 1]] <- new_row(
          v, length(xx), lv[1], paste(lv[-1], collapse = "/"),
          "One-way ANOVA", s[1, "F value"], "eta-squared", eta2,
          c(NA_real_, NA_real_), s[1, "Pr(>F)"], note)
      } else {
        kt <- stats::kruskal.test(xx ~ gg)
        h <- as.numeric(kt$statistic)
        eps2 <- (h - k + 1) / (length(xx) - k)
        note <- "rank-based; use pairwise Wilcoxon + mc_adjust_p() post-hoc"
        rows[[length(rows) + 1]] <- new_row(
          v, length(xx), lv[1], paste(lv[-1], collapse = "/"),
          "Kruskal-Wallis", h, "epsilon-squared", eps2,
          c(NA_real_, NA_real_), as.numeric(kt$p.value), note)
      }

    } else {
      tab <- table(gg, factor(xx))
      if (nrow(tab) < 2L || ncol(tab) < 2L) next
      n <- sum(tab)
      ct <- stats::chisq.test(tab, correct = FALSE)
      sparse <- any(ct$expected < 5)
      if (sparse) {
        ft <- .mc_with_seed(seed, stats::fisher.test(
          tab, simulate.p.value = TRUE, B = 2000))
        pval <- as.numeric(ft$p.value)
        test <- "Fisher exact (simulated, B=2000)"
        stat <- NA_real_
        v_est <- NA_real_
        note <- paste0("expected cell < 5; Cramer's V unavailable; ",
                       sum(ct$expected < 5), " sparse cells")
      } else {
        pval <- as.numeric(ct$p.value)
        test <- "Chi-square"
        stat <- as.numeric(ct$statistic)
        v_est <- .mc_cramers_v(stat, n, nrow(tab), ncol(tab))
        note <- "all expected counts >= 5"
      }
      rows[[length(rows) + 1]] <- new_row(
        v, n, lv[1], paste(lv[-1], collapse = "/"), test, stat,
        "Cramer's V", v_est, c(NA_real_, NA_real_), pval, note)
    }
  }

  if (!length(rows)) stop("No comparable variables with at least two groups.")
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$p.adjusted <- as.numeric(stats::p.adjust(out$p.value, method = p_adjust))
  out$significant <- !is.na(out$p.adjusted) & out$p.adjusted < alpha
  attr(out, "p_adjust") <- p_adjust
  attr(out, "alpha") <- alpha
  class(out) <- c("mc_group_comparison", class(out))
  out
}

#' Print a group-comparison table
#'
#' @param x An object of class `mc_group_comparison`.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.mc_group_comparison <- function(x, ...) {
  cat("Group comparison --", nrow(x), "tests, p adjusted by",
      attr(x, "p_adjust"), "\n")
  cols <- intersect(c("variable", "test", "effect", "estimate",
                      "conf.low", "conf.high", "p.value", "p.adjusted",
                      "significant"), names(x))
  d <- as.data.frame(x)[cols]
  for (nm in c("estimate", "conf.low", "conf.high"))
    if (nm %in% names(d)) d[[nm]] <- signif(d[[nm]], 3)
  for (nm in c("p.value", "p.adjusted"))
    if (nm %in% names(d)) d[[nm]] <- signif(d[[nm]], 3)
  print(d, row.names = FALSE)
  invisible(x)
}
