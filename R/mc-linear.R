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
# Module : linear regression (supplement to Chapter 5, common regression types)
# This file contains : mc_linear (simple/multiple linear regression, with
#              self-implemented diagnostics module)
# =============================================================================

#' Linear regression with self-implemented diagnostics
#'
#' Corresponds to Chapter 5 ("Linear regression", 5.2 & 5.6). Wraps
#' `stats::lm` and adds three *self-implemented* diagnostics: a bootstrap
#' confidence interval for coefficients, variance-inflation-factor (VIF)
#' multicollinearity checks via auxiliary regressions, and the four standard
#' regression diagnostic plots.
#'
#' @details
#' Purpose: fit simple or multiple linear regression and automatically produce a
#' reportable coefficient table and diagnostics.
#' Methods: ordinary least squares via `stats::lm`; on top of it, three
#' self-implemented diagnostics—
#' (1) *Bootstrap CI*: resample observations with replacement and refit
#' repeatedly, taking empirical quantiles of coefficients as the interval, with
#' no normality assumption (enable with `bootstrap > 0`);
#' (2) *VIF*: for each predictor, fit an auxiliary linear regression on the
#' remaining predictors and take 1/(1-R^2) as its variance inflation factor;
#' VIF > 5 usually signals multicollinearity;
#' (3) *Four diagnostic plots*: residuals vs fitted (non-linearity /
#' heteroscedasticity), normal Q-Q (residual normality), scale-location
#' (whether spread changes with fitted values), Cook's distance (influential
#' points).
#' Arguments: `data` is a data frame; `formula` is the model formula
#' (e.g. `y ~ x1 + x2`); `bootstrap` is the number of resamples (default 0 = no
#' bootstrap CI); `conf_level` is the confidence level, default 0.95.
#' Result: returns a list—`coef` (coefficient table: estimate/SE/t/p/normal CI,
#' plus boot lower/upper if bootstrap on), `performance` (R^2/adjusted R^2/F
#' statistic/AIC/BIC), `vif` (VIF per predictor), `diagnostics` (diagnostic
#' data frame), `diag_plots` (list of four ggplot objects), `plot` (main
#' residual diagnostic plot).
#' Note: VIF is computed only for numeric predictors; bootstrap significantly
#' increases runtime, so for very large samples try a small number (e.g. 200)
#' first. This function depends only on base R and ggplot2, no extra packages.
#' @param data A `data.frame`.
#' @param formula A model `formula`, e.g. `y ~ x1 + x2`.
#' @param bootstrap Number of bootstrap resamples; `0` = skip bootstrap CI.
#' @param conf_level Confidence level; default `0.95`.
#' @return A list with `coef`, `performance`, `vif`, `diagnostics`,
#'   `diag_plots`, `plot` (coefficient table, performance metrics, VIF,
#'   diagnostics, diagnostic plots and the main residual plot).
#' @examples
#' # Chapter 5, Sections 5.2/5.6: predict birth weight by mother's age, pregnancy weight, smoking
#' birthwt <- MASS::birthwt
#' res <- mc_linear(birthwt, bwt ~ age + lwt + smoke + factor(race),
#'                  bootstrap = 200)
#' res$coef                       # coefficient table (with bootstrap CI)
#' res$performance                # R^2 / adjusted R^2 / AIC / BIC
#' res$vif                        # variance inflation factor (multicollinearity)
#' print(res$plot)                # residuals vs fitted diagnostic plot
#' @export
mc_linear <- function(data, formula, bootstrap = 0L, conf_level = 0.95) {
  if (!inherits(data, "data.frame"))
    stop("`data` must be a data.frame.")
  if (!inherits(formula, "formula"))
    stop("`formula` must be a formula.")
  bootstrap <- as.integer(bootstrap)
  if (bootstrap < 0L) stop("`bootstrap` must be a non-negative integer.")

  fit <- stats::lm(formula, data = data)
  sm  <- summary(fit)
  cf  <- as.data.frame(sm$coefficients)
  names(cf) <- c("estimate", "std_error", "t_value", "p_value")
  cf$term <- rownames(cf)
  ci <- stats::confint(fit, level = conf_level)
  cf$ci_lower <- ci[, 1]
  cf$ci_upper <- ci[, 2]

  ## ---- self-implemented: bootstrap coefficient CI ----
  boot_ci <- NULL
  if (bootstrap > 0L) {
    n <- nrow(data)
    bmat <- matrix(NA_real_, bootstrap, length(coef(fit)))
    colnames(bmat) <- names(coef(fit))
    for (b in seq_len(bootstrap)) {
      idx <- sample.int(n, n, replace = TRUE)
      bf  <- tryCatch(stats::lm(formula, data = data[idx, , drop = FALSE]),
                      error = function(e) NULL)
      if (!is.null(bf)) bmat[b, names(coef(bf))] <- coef(bf)
    }
    a <- (1 - conf_level) / 2
    boot_ci <- apply(bmat, 2, function(v) {
      v <- v[!is.na(v)]
      if (length(v) < 2) c(NA_real_, NA_real_)
      else stats::quantile(v, c(a, 1 - a), na.rm = TRUE)
    })
    boot_ci <- as.data.frame(t(boot_ci))
    names(boot_ci) <- c("boot_lower", "boot_upper")
    cf <- cbind(cf, boot_ci[cf$term, , drop = FALSE])
  }
  cf <- cf[, c("term", "estimate", "std_error", "t_value", "p_value",
               "ci_lower", "ci_upper",
               if (bootstrap > 0L) c("boot_lower", "boot_upper") else NULL)]

  ## ---- self-implemented: VIF (auxiliary regression 1/(1-R^2)) ----
  mm <- stats::model.matrix(formula, data = data)
  preds <- colnames(mm)[colnames(mm) != "(Intercept)"]
  vif <- data.frame(term = preds, VIF = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_along(preds)) {
    others <- setdiff(preds, preds[i])
    if (length(others) == 0) { vif$VIF[i] <- 1; next }
    r2 <- tryCatch({
      ar <- stats::lm(mm[, preds[i]] ~ mm[, others])
      summary(ar)$r.squared
    }, error = function(e) NA_real_)
    vif$VIF[i] <- if (is.na(r2) || r2 >= 1) NA_real_ else 1 / (1 - r2)
  }

  ## ---- performance metrics ----
  performance <- data.frame(
    metric = c("r_squared", "adj_r_squared", "f_statistic", "aic", "bic"),
    value  = c(sm$r.squared, sm$adj.r.squared,
               sm$fstatistic["value"], stats::AIC(fit), stats::BIC(fit)),
    stringsAsFactors = FALSE)

  ## ---- diagnostic data ----
  fitted <- as.numeric(fit$fitted.values)
  resid  <- as.numeric(fit$residuals)
  stdres <- as.numeric(stats::rstandard(fit))
  cooks  <- as.numeric(stats::cooks.distance(fit))
  diag_df <- data.frame(fitted = fitted, resid = resid,
                        std_resid = stdres, cooks = cooks,
                        stringsAsFactors = FALSE)

  ## ---- four diagnostic plots (ggplot2) ----
  p1 <- ggplot2::ggplot(diag_df, ggplot2::aes(fitted, resid)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::geom_smooth(se = FALSE, color = "blue", size = 0.8) +
    ggplot2::labs(x = "Fitted values", y = "Residuals") + mc_theme()
  p2 <- ggplot2::ggplot(diag_df, ggplot2::aes(sample = std_resid)) +
    ggplot2::stat_qq(size = 1.2, alpha = 0.6) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red",
                         linetype = "dashed") +
    ggplot2::labs(x = "Theoretical quantiles", y = "Standardized residuals") +
    mc_theme()
  p3 <- ggplot2::ggplot(diag_df, ggplot2::aes(fitted, sqrt(abs(std_resid)))) +
    ggplot2::geom_point(size = 1.2, alpha = 0.6) +
    ggplot2::geom_smooth(se = FALSE, color = "blue", size = 0.8) +
    ggplot2::labs(x = "Fitted values", y = "Sqrt(|standardized residuals|)") +
    mc_theme()
  p4 <- ggplot2::ggplot(diag_df, ggplot2::aes(seq_along(cooks), cooks)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    ggplot2::labs(x = "Observation index", y = "Cook's distance") + mc_theme()

  list(coef = cf, performance = performance, vif = vif,
       diagnostics = diag_df, diag_plots = list(p1, p2, p3, p4),
       plot = p1, fit = fit)
}
