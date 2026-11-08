# =============================================================================
# Software name : macodgent medical data-science research toolkit
# Copyright holder : Dr.Xie
# Version : 1.0.2
# Language : R
# -----------------------------------------------------------------------------
# Module : logistic model report (full diagnostic suite)
# mc_logistic() deliberately stays a thin coefficient table. This module adds
# the model-level quantities a clinical paper must report: discrimination,
# calibration, collinearity, influential observations, events per variable and
# an explicit methodological verdict for each check.
# =============================================================================

# Internal: AUC from the Mann-Whitney identity, with average ranks for ties.
.mc_auc <- function(p, y) {
  ok <- !is.na(p) & !is.na(y)
  p <- as.numeric(p[ok]); y <- as.numeric(y[ok])
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  (sum(rank(p)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Full diagnostic report for a binary logistic regression
#'
#' Fits one binary logistic model and returns everything needed to report and
#' defend it: OR table, discrimination (AUC with CI), Hosmer-Lemeshow
#' calibration and calibration curve, collinearity (VIF), influential
#' observations, events per variable, global fit statistics, and an explicit
#' pass/fail verdict for each methodological check.
#'
#' @details
#' Purpose: [mc_logistic()] returns only a coefficient table. A published
#' logistic model additionally requires evidence that it discriminates,
#' calibrates, is not overfitted and is not driven by a handful of influential
#' records. This function assembles all of that in one reproducible call, using
#' base R plus the existing package dependencies (no pROC, car or rms).
#'
#' Methods:
#' \itemize{
#' \item Discrimination: AUC via the Mann-Whitney U identity, equivalent to
#' the concordance index. Its confidence interval is a stratified-free
#' percentile bootstrap over observation pairs (\code{auc_boot} resamples);
#' set \code{auc_boot = 0} to use the Hanley-McNeil normal approximation.
#' \item Calibration: Hosmer-Lemeshow test on prediction deciles with
#' \code{g - 2} degrees of freedom, plus a decile-level observed-versus-
#' expected table and calibration plot.
#' \item Collinearity: approximate VIF from the inverse correlation matrix of
#' the design-matrix columns (excluding the intercept). This is a screening
#' approximation; a singular matrix yields \code{NA}.
#' \item Influence: standardised deviance residuals and Cook's distance;
#' records with \code{|residual| > 2} or Cook's D \code{> 4/n} are flagged.
#' \item Overfitting: events per variable (EPV) = number of minority outcome
#' events / number of estimated coefficients excluding the intercept.
#' \item Global fit: likelihood-ratio chi-square against the null model,
#' Nagelkerke R-squared, Brier score, AIC and BIC.
#' }
#'
#' Arguments: `data` is a data frame; `formula` is `outcome ~ predictors`.
#' The outcome must be 0/1, logical, or a two-level factor (the second level is
#' taken as the event). `conf` sets all confidence levels; `auc_boot` is the
#' number of bootstrap resamples for the AUC CI; `seed` makes the bootstrap
#' reproducible without disturbing the caller's RNG state.
#'
#' Result: an object of class `mc_logistic_report` (a list) with
#' `coefficients`, `performance`, `checks`, `conclusion`, `calibration`,
#' `roc`, `vif`, `influential`, `model`, `plot` and `roc_plot`. Printing the
#' object shows a compact summary rather than dumping the whole list.
#'
#' Note: Hosmer-Lemeshow is sensitive to the number of groups and loses power
#' in small samples; a non-significant result is not proof of good calibration.
#' Inspect `plot` alongside the test. EPV < 10 signals overfitting risk and
#' cannot be repaired by reporting more diagnostics. The VIF here is a
#' design-matrix approximation, not the weighted GLM VIF of package `car`.
#'
#' @param data A `data.frame`.
#' @param formula A model formula, e.g. `low ~ age + smoke + lwt`.
#' @param conf Confidence level; default `0.95`.
#' @param auc_boot Bootstrap resamples for the AUC CI; `0` = analytic SE.
#' @param seed Random seed for the AUC bootstrap.
#' @return An object of class `mc_logistic_report`; see Details.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' rep <- mc_logistic_report(d,
#'   bp_control ~ age + bmi + treatment + diabetes, auc_boot = 0L)
#' rep$checks
#' rep$performance
#' @seealso \code{\link{mc_logistic}}, \code{\link{mc_stepwise}},
#'   \code{\link{mc_compare_groups}}.
#' @export
mc_logistic_report <- function(data, formula, conf = 0.95,
                               auc_boot = 500L, seed = 1L) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  fit <- eval(substitute(stats::glm(formula, data = DATA,
                                    family = stats::binomial),
                         list(DATA = data)))

  mf <- stats::model.frame(fit)
  yraw <- stats::model.response(mf)
  if (is.factor(yraw) || is.character(yraw)) {
    ylev <- levels(factor(yraw))
    if (length(ylev) != 2L)
      stop("The outcome must have exactly two levels for a binary model.")
    y <- as.integer(factor(yraw) == ylev[2])
  } else {
    y <- as.numeric(yraw)
    if (!all(y[!is.na(y)] %in% c(0, 1)))
      stop("The numeric outcome must be coded 0/1.")
  }
  p <- as.numeric(stats::fitted(fit))
  ok <- !is.na(p) & !is.na(y)
  p <- p[ok]; y <- y[ok]
  n <- length(y)
  n_event <- sum(y == 1); n_nonevent <- sum(y == 0)
  n_coef <- length(stats::coef(fit))
  epv <- min(n_event, n_nonevent) / max(1L, n_coef - 1L)

  ## ---- coefficient table ----
  td <- broom::tidy(fit, conf.int = TRUE, conf.level = conf,
                    exponentiate = TRUE)

  ## ---- discrimination ----
  auc <- .mc_auc(p, y)
  zc <- stats::qnorm(1 - (1 - conf) / 2)
  if (auc_boot > 0L) {
    bs <- .mc_with_seed(seed, replicate(auc_boot, {
      i <- sample.int(n, n, replace = TRUE)
      a <- .mc_auc(p[i], y[i])
      if (is.finite(a)) a else NA_real_
    }))
    bs <- bs[is.finite(bs)]
    auc_ci <- if (length(bs) >= 20L) {
      as.numeric(stats::quantile(bs, c((1 - conf) / 2, 1 - (1 - conf) / 2),
                                 na.rm = TRUE))
    } else c(NA_real_, NA_real_)
    auc_method <- paste0("percentile bootstrap (", length(bs), " resamples)")
  } else {
    # Hanley-McNeil standard error of the Mann-Whitney AUC.
    q1 <- auc / (2 - auc); q2 <- 2 * auc^2 / (1 + auc)
    se <- sqrt((auc * (1 - auc) + (n_event - 1) * (q1 - auc^2) +
                  (n_nonevent - 1) * (q2 - auc^2)) / (n_event * n_nonevent))
    auc_ci <- c(auc - zc * se, auc + zc * se)
    auc_method <- "Hanley-McNeil approximation"
  }

  ## ---- calibration ----
  br <- unique(stats::quantile(p, probs = seq(0, 1, length.out = 11)))
  grp <- cut(p, breaks = br, include.lowest = TRUE, labels = FALSE)
  ng <- as.numeric(tapply(y, grp, length))
  obs <- as.numeric(tapply(y, grp, sum))
  expc <- as.numeric(tapply(p, grp, sum))
  cal <- data.frame(decile = seq_along(ng), n = ng,
                    mean_predicted = expc / ng,
                    mean_observed = obs / ng,
                    observed_events = obs, expected_events = expc,
                    stringsAsFactors = FALSE)
  denom <- expc * (1 - expc / ng)
  usable <- is.finite(denom) & denom > 0
  g <- length(ng)
  hl <- if (sum(usable) >= 3L) {
    sum((obs[usable] - expc[usable])^2 / denom[usable])
  } else NA_real_
  hl_df <- max(1L, g - 2L)
  hl_p <- if (is.finite(hl)) stats::pchisq(hl, hl_df, lower.tail = FALSE) else NA_real_
  brier <- mean((p - y)^2)

  ## ---- collinearity ----
  X <- stats::model.matrix(fit)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (ncol(X) >= 2L) {
    vv <- tryCatch(diag(solve(stats::cor(X))), error = function(e)
      rep(NA_real_, ncol(X)))
    vif <- data.frame(term = colnames(X), VIF = as.numeric(vv),
                      stringsAsFactors = FALSE)
  } else {
    vif <- data.frame(term = colnames(X),
                      VIF = rep(1, max(0L, ncol(X))),
                      stringsAsFactors = FALSE)
  }
  max_vif <- suppressWarnings(max(vif$VIF, na.rm = TRUE))
  if (!is.finite(max_vif)) max_vif <- NA_real_

  ## ---- influence ----
  sr <- as.numeric(stats::rstandard(fit))
  cd <- as.numeric(stats::cooks.distance(fit))
  infl <- data.frame(row = seq_along(sr), std_residual = sr,
                     cooks_distance = cd,
                     flagged = abs(sr) > 2 | cd > 4 / nrow(mf),
                     stringsAsFactors = FALSE)
  n_flag <- sum(infl$flagged)

  ## ---- global fit ----
  lr <- fit$null.deviance - fit$deviance
  lr_df <- fit$df.null - fit$df.residual
  lr_p <- stats::pchisq(lr, lr_df, lower.tail = FALSE)
  nagel <- (1 - exp((fit$deviance - fit$null.deviance) / n)) /
    (1 - exp(-fit$null.deviance / n))
  perf <- data.frame(
    metric = c("n", "events", "non_events", "n_coefficients", "EPV",
               "AUC", "AUC_CI_low", "AUC_CI_high", "AUC_method",
               "Brier_score", "HosmerLemeshow_chi2", "HosmerLemeshow_df",
               "HosmerLemeshow_p", "Nagelkerke_R2",
               "LikelihoodRatio_chi2", "LikelihoodRatio_df",
               "LikelihoodRatio_p", "AIC", "BIC", "max_VIF",
               "influential_rows"),
    value = c(n, n_event, n_nonevent, n_coef, epv,
              auc, auc_ci[1], auc_ci[2], auc_method,
              brier, hl, hl_df, hl_p, nagel,
              lr, lr_df, lr_p, stats::AIC(fit), stats::BIC(fit),
              max_vif, n_flag),
    stringsAsFactors = FALSE)

  ## ---- ROC curve ----
  thr <- seq(0, 1, length.out = 101)
  roc <- data.frame(
    threshold = thr,
    tpr = vapply(thr, function(t) sum(y == 1 & p >= t) / max(1L, n_event),
                 numeric(1)),
    fpr = vapply(thr, function(t) sum(y == 0 & p >= t) / max(1L, n_nonevent),
                 numeric(1)),
    stringsAsFactors = FALSE)

  ## ---- plots ----
  cal_plot <- ggplot2::ggplot(cal, ggplot2::aes(
    x = mean_predicted, y = mean_observed, size = n)) +
    ggplot2::geom_point(color = "#1A439C", alpha = 0.8) +
    ggplot2::geom_line(color = "#1A439C", linewidth = 0.6) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                         color = "#C0392B", linewidth = 0.6) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "Mean predicted probability",
                  y = "Observed event proportion",
                  title = "Calibration by prediction decile") + mc_theme()
  roc_plot <- ggplot2::ggplot(roc, ggplot2::aes(x = fpr, y = tpr)) +
    ggplot2::geom_line(color = "#1C7331", linewidth = 0.8) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                         color = "#8A94A6", linewidth = 0.6) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "1 - Specificity", y = "Sensitivity",
                  title = sprintf("ROC (AUC = %.3f)", auc)) + mc_theme()

  ## ---- explicit methodological checks ----
  discr <- if (!is.finite(auc)) "unknown" else if (auc >= 0.8) "good" else
    if (auc >= 0.7) "acceptable" else "poor"
  checks <- data.frame(
    check = c("Discrimination", "Calibration", "Overfitting (EPV)",
              "Multicollinearity", "Influential observations"),
    value = c(auc, hl_p, epv, max_vif, n_flag),
    threshold = c("AUC >= 0.70", "HL p > 0.05", "EPV >= 10",
                  "max VIF < 5", "flagged = 0"),
    status = c(
      if (!is.finite(auc)) "UNDETERMINED" else if (auc >= 0.7) "PASS" else "REVIEW",
      if (!is.finite(hl_p)) "UNDETERMINED" else if (hl_p > 0.05) "PASS" else "REVIEW",
      if (epv >= 10) "PASS" else "REVIEW",
      if (!is.finite(max_vif)) "UNDETERMINED" else if (max_vif < 5) "PASS" else "REVIEW",
      if (n_flag == 0L) "PASS" else "REVIEW"),
    stringsAsFactors = FALSE)
  conclusion <- c(
    sprintf("Discrimination: AUC = %.3f (%.0f%% CI %.3f-%.3f), %s; %s.",
            auc, 100 * conf, auc_ci[1], auc_ci[2], discr, auc_method),
    sprintf("Calibration: Hosmer-Lemeshow chi2 = %.2f on %d df, p = %s (%s).",
            hl, hl_df, if (is.finite(hl_p)) format.pval(hl_p, digits = 3) else "NA",
            if (!is.finite(hl_p)) "not computable" else
              if (hl_p > 0.05) "no evidence of miscalibration" else "evidence of miscalibration"),
    sprintf("Overfitting: %d events and %d predictors give EPV = %.1f; %s.",
            min(n_event, n_nonevent), n_coef - 1L, epv,
            if (epv >= 10) "above the conventional threshold of 10" else
              "below 10, consider reducing predictors or penalisation"),
    sprintf("Collinearity: maximum approximate VIF = %s (%s).",
            if (is.finite(max_vif)) sprintf("%.2f", max_vif) else "NA",
            if (!is.finite(max_vif)) "design matrix singular" else
              if (max_vif < 5) "no strong collinearity" else "collinearity warrants review"),
    sprintf("Influence: %d of %d observations flagged (|std residual| > 2 or Cook's D > 4/n).",
            n_flag, nrow(mf)))

  out <- list(coefficients = td, performance = perf, checks = checks,
              conclusion = conclusion, calibration = cal, roc = roc,
              vif = vif, influential = infl, model = fit,
              plot = cal_plot, roc_plot = roc_plot)
  class(out) <- "mc_logistic_report"
  out
}

#' Print a compact logistic model report
#'
#' @param x An object of class `mc_logistic_report`.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.mc_logistic_report <- function(x, ...) {
  cat("Logistic model report\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  p <- x$performance
  getv <- function(k) p$value[p$metric == k]
  cat(sprintf("N = %s   Events = %s   Predictors = %s   EPV = %.1f\n",
              getv("n"), getv("events"),
              as.integer(as.numeric(getv("n_coefficients"))) - 1L,
              as.numeric(getv("EPV"))))
  cat(sprintf("AUC = %.3f (%.3f-%.3f)   Nagelkerke R2 = %.3f   Brier = %.3f\n",
              as.numeric(getv("AUC")), as.numeric(getv("AUC_CI_low")),
              as.numeric(getv("AUC_CI_high")),
              as.numeric(getv("Nagelkerke_R2")),
              as.numeric(getv("Brier_score"))))
  cat("\nMethodological checks\n")
  print(x$checks[, c("check", "value", "status")], row.names = FALSE,
        digits = 3)
  cat("\n")
  cat(paste0("- ", x$conclusion, collapse = "\n"), "\n")
  cat("\nCoefficient table (first rows)\n")
  print(utils::head(x$coefficients, 8L), row.names = FALSE, digits = 3)
  invisible(x)
}
