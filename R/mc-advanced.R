# =============================================================================
# Software name: macodgent medical-data-science research toolkit
# Copyright owner: Dr.Xie
# Version: 1.0.0
# Language: R
# Completion date: 2026-10-30
# Rights notice: This program is independently developed by the copyright owner
# and protected under the Computer Software Protection Regulation;
# without written permission, it may not be copied, distributed,
# rented, or used for commercial purposes.
# -----------------------------------------------------------------------------
# Module: Advanced modeling module (advanced)
# This file contains: mc_gam / mc_lasso / mc_ridge / mc_elastic_net / mc_rcs /
# mc_quantile_reg / mc_pca_regression (principal components regression) /
# mc_nnet_reg (neural-network regression)
# where.rcs_basis (restricted cubic spline basis) is an original
# implementation per Harrell, independent of the rms package;
# mc_rcs builds this spline internally and fits the prediction curve.
# =============================================================================

#' Generalized Additive Model (GAM)
#'
#' Generalized additive regression models (GAM) capture non-linear relationships while keeping interpretability.
#' Built on \code{\link[mgcv]{gam}}, it returns the
#' degrees of freedom and significance of each smooth term.
#'
#' @details
#' Purpose: Fit a generalized additive model (GAM) that, while keeping the
#' Interpretability of a linear model, lets predictors enter the model as smooth
#' (non-linear) terms, automatically capturing curved relationships.
#' Method: Built on mgcv's \code{gam()}, using splines (e.g. \code{s(x)}) for the
#' Smooth terms; it auto-selects the smoothing degrees of freedom and reports
#' Significance.
#' Arguments: \code{data} is a data frame; \code{formula} is of the form
#' \code{y ~ s(x1) + x2}, where \code{s()}-wrapped terms are smooth and the rest are linear.
#' Results: Returns a list with \code{model} (a gam object; use \code{plot()} to view
#' The smooth curves) and \code{terms} (edf, reference df, F value and p value of
#' Each smooth term, used to judge whether the variable has a non-linear effect).
#' Notes: If \code{edf} is near 1 and not significant, the variable is approximately
#' Linear and smoothing is unnecessary; a very large \code{edf} may overfit. Smooth
#' Terms need sufficient sample size for stable estimation. Depends on mgcv (Suggests).
#' To view the effect plot, use \code{plot(res$model, pages = 1, shade = TRUE)}.
#' @param data A `data.frame`.
#' @param formula A GAM formula such as \code{y ~ s(x1) + x2}.
#' @param ... Additional arguments passed to \code{\link[mgcv]{gam}}.
#' @return A list with \code{model} and \code{terms}.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' \donttest{
#' g <- mc_gam(d, sbp_week12 ~ s(sbp_baseline) + treatment)
#' g$smooths
#' g$performance
#' }
#' @seealso \code{\link{mc_rcs}} (restricted cubic spline with self-implemented basis).
#' @export
mc_gam <- function(data, formula, ...) {
 if (!requireNamespace("mgcv", quietly = TRUE)) {
 stop("`mc_gam` requires the mgcv package (shipped with R).")
 }
 fit <- eval(substitute(mgcv::gam(formula, data = DATA, ...),
 list(DATA = data)))
 s <- summary(fit)
 smooths <- tryCatch({
 t <- as.data.frame(s$s.table)
 t$term <- rownames(t)
 rownames(t) <- NULL
 t
 }, error = function(e) NULL)
 if (is.null(smooths) || !nrow(smooths))
 smooths <- data.frame(edf = NA_real_, Ref.df = NA_real_,
 F = NA_real_, `p-value` = NA_real_,
 term = "none", check.names = FALSE)
 names(smooths) <- c("edf", "ref_df", "F", "p.value", "term")
 parametric <- tryCatch({
 t <- as.data.frame(s$p.table)
 t$term <- rownames(t)
 rownames(t) <- NULL
 t
 }, error = function(e) NULL)
 n <- tryCatch(as.numeric(stats::nobs(fit)), error = function(e) nrow(data))
 perf <- data.frame(
 metric = c("AIC", "BIC", "deviance_explained", "adj_r_squared",
 "n", "total_edf"),
 value = c(stats::AIC(fit), stats::BIC(fit), s$dev.expl, s$r.sq,
 n, sum(smooths$edf, na.rm = TRUE)),
 stringsAsFactors = FALSE)
 p_sm <- suppressWarnings(min(smooths$p.value, na.rm = TRUE))
 edf_max <- suppressWarnings(max(smooths$edf, na.rm = TRUE))
 checks <- rbind(
 .mc_check_row("Smooth terms fitted", sum(!is.na(smooths$p.value)),
 ">= 1", sum(!is.na(smooths$p.value)) >= 1L),
 .mc_check_row("Non-linear signal", p_sm, "min smooth p < 0.05",
 if (is.finite(p_sm)) p_sm < 0.05 else NA,
 "edf near 1 with p > 0.05 suggests a linear term is sufficient"),
 .mc_check_row("Smooth complexity", edf_max, "max edf < n / 10",
 if (is.finite(edf_max)) edf_max < n / 10 else NA))
 list(model = fit, terms = broom::tidy(fit), smooths = smooths,
 parametric = parametric, performance = perf, checks = checks)
}

# Internal: unified implementation for elastic-net models.
.elastic_net <- function(data, x, y, alpha, nfolds = 10,
                        family = "gaussian", stability = 0L, seed = NULL) {
 if (!requireNamespace("glmnet", quietly = TRUE)) {
 stop("`mc_elastic_net` requires the glmnet package; please install.packages('glmnet').")
 }
 xm <- as.matrix(data[x])
 yv <- data[[y]]
 if (!is.null(seed)) set.seed(seed)
 cv <- glmnet::cv.glmnet(xm, yv, alpha = alpha, nfolds = nfolds,
 family = family)
 cf <- as.matrix(stats::coef(cv, s = "lambda.min"))
 df <- data.frame(variable = rownames(cf), coefficient = cf[, 1],
 stringsAsFactors = FALSE)
 df <- df[df$variable != "(Intercept)", ]
 rownames(df) <- NULL
 df$selected <- abs(df$coefficient) > 0
 cf1 <- as.matrix(stats::coef(cv, s = "lambda.1se"))
 cv_table <- data.frame(lambda = as.numeric(cv$lambda),
 cvm = as.numeric(cv$cvm), cvsd = as.numeric(cv$cvsd),
 cv_lo = as.numeric(cv$cvm) - as.numeric(cv$cvsd),
 cv_hi = as.numeric(cv$cvm) + as.numeric(cv$cvsd),
 stringsAsFactors = FALSE)
 performance <- data.frame(
 metric = c("lambda.min", "lambda.1se", "cv_mse_min",
 "n_nonzero_min", "n_nonzero_1se"),
 value = c(cv$lambda.min, cv$lambda.1se, min(cv$cvm, na.rm = TRUE),
 sum(abs(cf[, 1]) > 0) - 1L, sum(abs(cf1[, 1]) > 0) - 1L),
 stringsAsFactors = FALSE)
 stab <- NULL
 if (stability > 0L) {
 if (!is.null(seed)) set.seed(seed + 1L)
 cnt <- setNames(rep(0L, nrow(df)), df$variable)
 n <- nrow(xm)
 for (b in seq_len(stability)) {
 idx <- sample.int(n, n, replace = TRUE)
 c2 <- tryCatch(glmnet::cv.glmnet(xm[idx, , drop = FALSE], yv[idx],
 alpha = alpha, nfolds = nfolds,
 family = family), error = function(e) NULL)
 if (is.null(c2)) next
 cc <- as.matrix(stats::coef(c2, s = "lambda.min"))
 nz <- setdiff(rownames(cc)[abs(cc[, 1]) > 0], "(Intercept)")
 hit <- intersect(nz, names(cnt))
 cnt[hit] <- cnt[hit] + 1L
 }
 stab <- data.frame(variable = names(cnt),
 selection_frequency = as.numeric(cnt) / stability,
 stringsAsFactors = FALSE)
 }
 list(cv = cv, coefficients = df, lambda.min = cv$lambda.min,
 lambda.1se = cv$lambda.1se, alpha = alpha,
 cv_table = cv_table, performance = performance, stability = stab)
}

#' Lasso Regression
#'
#' Lasso regression (L1 regularization) for prediction and embedded variable selection.
#' Built on \pkg{glmnet} with \code{alpha = 1}. Returns the non-zero
#' coefficient table selected by 10-fold cross-validated lambda.
#'
#' @details
#' Purpose: Use Lasso (L1 regularization) for regression and variable selection; in
#' High-dimensional or collinear data it automatically shrinks the coefficients of
#' Unimportant variables to zero, achieving embedded feature screening.
#' Method: Based on the glmnet package; \code{alpha = 1} is Lasso; k-fold
#' Cross-validation (default 10 folds) is used to select the optimal penalty
#' Strength lambda.
#' Arguments: \code{data} is a data frame; \code{x} is a vector of predictor column
#' Names; \code{y} is the response column name; \code{nfolds} is the number of CV
#' Folds (default 10); \code{family} defaults to gaussian (continuous outcome), and
#' Can also be binomial, etc.
#' Results: Returns a list with \code{coefficients} (table of non-zero coefficients;
#' Variables shrunk to zero do not appear).
#' Notes: Lasso coefficients are biased (shrunken overall), so it suits prediction
#' And screening but not direct statistical inference. When several variables are
#' Highly correlated, Lasso randomly keeps one and zeroes the rest. Larger lambda
#' Means stronger shrinkage; tune with cv.glmnet's lambda.min / lambda.1se. Depends
#' On glmnet (Suggests).
#' @param data A `data.frame`.
#' @param x Character vector of predictor column names.
#' @param y Response variable name.
#' @param nfolds Cross-validation folds; default `10`.
#' @param family GLM family; `"gaussian"` (default) or `"binomial"`, etc.
#' @return A list with \code{cv}, \code{coefficients}, \code{lambda.min}.
#' @param stability Number of bootstrap refits for coefficient-selection frequency; `0` (default) skips it.
#' @param seed Random seed for cross-validation and stability resamples.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' \donttest{
#' xvars <- c("age", "bmi", "sbp_baseline", "ldl", "glucose", "creatinine")
#' mc_lasso(d, x = xvars, y = "sbp_week12", nfolds = 5, seed = 1)$coefficients
#' }
#' @seealso \code{\link{mc_ridge}}, \code{\link{mc_elastic_net}}.
#' @export
mc_lasso <- function(data, x, y, nfolds = 10, family = "gaussian",
               stability = 0L, seed = NULL) {
 .elastic_net(data, x, y, alpha = 1, nfolds = nfolds,
                family = family, stability = stability,
                seed = seed)
}

#' Ridge Regression
#'
#' Ridge regression (L2 regularization) handles multicollinearity by shrinking all coefficients proportionally.
#' Built on \pkg{glmnet} with \code{alpha = 0}. Returns the coefficient
#' table (ridge coefficients are typically non-sparse).
#'
#' @details
#' Purpose: Use ridge regression (L2 regularization) to handle multicollinearity and
#' Stabilize coefficient estimates, suited to cases where all variables matter and
#' You do not want to drop any.
#' Method: Based on glmnet; \code{alpha = 0} is ridge; it shrinks all coefficients
#' Proportionally (without zeroing) and keeps every variable.
#' Arguments: \code{data} is a data frame; \code{x} is a vector of predictor column
#' Names; \code{y} is the response column name; \code{nfolds} is the number of CV
#' Folds (default 10); \code{family} defaults to gaussian.
#' Results: Returns a list with \code{coefficients} (shrunk but all non-zero coefficients).
#' Notes: Unlike Lasso, ridge coefficients are never exactly zero, so it performs no
#' Variable selection, only shrinkage. It is robust to collinearity (spreading
#' Variance across correlated variables), but all variables remain in the
#' Interpretation. Larger lambda means stronger shrinkage. Depends on glmnet (Suggests).
#' @param data A `data.frame`.
#' @param x Character vector of predictor column names.
#' @param y Response variable name.
#' @param nfolds Cross-validation folds; default `10`.
#' @param family GLM family; `"gaussian"` (default) or `"binomial"`, etc.
#' @return A list with \code{cv}, \code{coefficients}, \code{lambda.min}.
#' @param stability Number of bootstrap refits for coefficient-selection frequency; `0` (default) skips it.
#' @param seed Random seed for cross-validation and stability resamples.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' \donttest{
#' xvars <- c("age", "bmi", "sbp_baseline", "ldl", "glucose", "creatinine")
#' mc_ridge(d, x = xvars, y = "sbp_week12", nfolds = 5, seed = 1)$coefficients
#' }
#' @seealso \code{\link{mc_lasso}}, \code{\link{mc_elastic_net}}.
#' @export
mc_ridge <- function(data, x, y, nfolds = 10, family = "gaussian",
               stability = 0L, seed = NULL) {
 .elastic_net(data, x, y, alpha = 0, nfolds = nfolds,
                family = family, stability = stability,
                seed = seed)
}

#' Elastic Net Regression
#'
#' Elastic net blends L1 and L2 penalties, balancing variable selection and shrinkage.
#' Built on \pkg{glmnet}; \code{alpha} tunes between
#' ridge (0) and lasso (1).
#'
#' @details
#' Purpose: Use elastic net to combine Lasso's sparsity with ridge's robustness to
#' Collinearity; it is a continuous interpolation between the two.
#' Method: Based on glmnet; \code{alpha} tunes the mixing ratio between 0 (pure
#' Ridge) and 1 (pure Lasso); \code{alpha = 0.5} is a common compromise.
#' Arguments: \code{data} is a data frame; \code{x} is a vector of predictor column
#' Names; \code{y} is the response column name; \code{alpha} is the mixing parameter
#' (default 0.5); \code{nfolds} is the number of CV folds; \code{family} defaults to gaussian.
#' Results: Returns a list with \code{coefficients} (coefficients shrunk per the
#' Chosen alpha, some possibly zero).
#' Notes: When predictors are grouped and highly correlated, elastic net is more
#' Robust than Lasso: it keeps or drops the whole group together rather than
#' Randomly keeping one. The closer alpha is to 1, the more it favours selection;
#' The closer to 0, the more it favours shrinkage. Depends on glmnet (Suggests).
#' @param data A `data.frame`.
#' @param x Character vector of predictor column names.
#' @param y Response variable name.
#' @param alpha Elastic-net mixing parameter; default `0.5`.
#' @param nfolds Cross-validation folds; default `10`.
#' @param family GLM family; `"gaussian"` (default) or `"binomial"`, etc.
#' @return A list with \code{cv}, \code{coefficients}, \code{lambda.min}, \code{alpha}.
#' @param stability Number of bootstrap refits for coefficient-selection frequency; `0` (default) skips it.
#' @param seed Random seed for cross-validation and stability resamples.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' \donttest{
#' xvars <- c("age", "bmi", "sbp_baseline", "ldl", "glucose", "creatinine")
#' mc_elastic_net(d, x = xvars, y = "sbp_week12",
#'                nfolds = 5, seed = 1)$coefficients
#' }
#' @seealso \code{\link{mc_lasso}}, \code{\link{mc_ridge}}.
#' @export
mc_elastic_net <- function(data, x, y, alpha = 0.5, nfolds = 10, family = "gaussian",
                           stability = 0L, seed = NULL) {
 .elastic_net(data, x, y, alpha = alpha, nfolds = nfolds,
                family = family, stability = stability,
                seed = seed)
}

# Internal: restricted cubic spline basis (Harrell's definition).
.rcs_basis <- function(x, knots = NULL, nk = 3) {
 x <- as.numeric(x)
 if (is.null(knots)) {
 qp <- c(0, (1:nk) / (nk + 1), 1)
 knots <- stats::quantile(x, probs = qp, na.rm = TRUE, names = FALSE)
 }
 knots <- sort(unique(knots))
 K <- length(knots)
 t1 <- knots[1]; tK <- knots[K]
 out <- data.frame(x = x); names(out) <- "rcs_lin"
 for (i in 2:(K - 1)) {
 ti <- knots[i]
 lambda_i <- (tK - ti) / (tK - t1)
 out[[paste0("rcs_", i)]] <- pmax(0, x - ti)^3 -
 lambda_i * pmax(0, x - t1)^3 - (1 - lambda_i) * pmax(0, x - tK)^3
 }
 attr(out, "knots") <- knots
 out
}

#' Restricted Cubic Spline Regression
#'
#' Restricted cubic splines model flexible non-linear dose-response relationships without external packages.
#' Without depending on \pkg{rms}, it builds a
#' restricted cubic spline basis and fits a linear / generalized linear model,
#' returning the prediction curve.
#'
#' @details
#' Purpose: Use restricted cubic splines (RCS) to describe the non-linear
#' Dose-response relationship between a predictor and the outcome, and to visualize
#' The curve; a common approach in epidemiology for presenting exposure-effect curves.
#' Method: A self-implemented spline basis (no extra dependencies) applies linear
#' Constraints at both ends of the predictor ("restricted") to avoid wild swings at
#' The boundaries; default interior knots \code{nk = 3}.
#' Arguments: \code{data} is a data frame; \code{x} is the predictor column name;
#' \code{y} is the response column name; \code{nk} is the number of knots (default 3,
#' Larger = more flexible); \code{knots} can manually set knot positions (default
#' Auto-selected by quantiles); \code{family} defaults to gaussian and also supports
#' Binomial (then it describes how the log-odds change with x).
#' Results: Returns a list with \code{knots} (knot positions), \code{coefficients}
#' (spline basis coefficients), and \code{plot} (prediction-curve ggplot object).
#' Notes: More knots make the curve wigglier and more prone to overfitting; usually
#' 3-5 knots suffice. RCS does not assume monotonicity and can show U-shaped or
#' Inverted-U relationships. View the result with \code{print(res$plot)}.
#' @param data A `data.frame`.
#' @param x Numeric predictor name.
#' @param y Response variable name.
#' @param nk Number of interior knots; default `3`.
#' @param knots Full knot vector including boundaries; `NULL` = automatic quantile knots.
#' @param family GLM family; `"gaussian"` (default) or `"binomial"`, etc.
#' @return A list with \code{model}, \code{coefficients}, \code{plot}, \code{knots}.
#' @param conf Confidence level for the curve band and coefficient intervals; default `0.95`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' r <- mc_rcs(d, x = "sbp_baseline", y = "sbp_week12", nk = 4)
#' r$knots
#' r$nonlinearity_p
#' head(r$curve)
#' @seealso \code{\link{mc_gam}} (smooth terms are also used for non-linear relationship modeling).
#' @export
mc_rcs <- function(data, x, y, nk = 3, knots = NULL, family = "gaussian",
                   conf = 0.95) {
  bx <- .rcs_basis(data[[x]], knots = knots, nk = nk)
  df <- data.frame(bx, y = data[[y]])
  f <- stats::as.formula(paste("y ~", paste(names(bx), collapse = " + ")))
  is_lm <- is.character(family) && length(family) == 1L &&
    identical(family, "gaussian")
  fit <- if (is_lm) stats::lm(f, data = df) else
    stats::glm(f, data = df, family = family)
  grid <- data.frame(x = seq(min(data[[x]], na.rm = TRUE),
                             max(data[[x]], na.rm = TRUE), length.out = 100))
  gx <- .rcs_basis(grid$x, knots = attr(bx, "knots"), nk = nk)
  zc <- stats::qnorm(1 - (1 - conf) / 2)
  if (is_lm) {
    pred <- stats::predict(fit, newdata = gx, se.fit = TRUE)
    pc <- data.frame(x = grid$x, fit = as.numeric(pred$fit),
                     lo = as.numeric(pred$fit) - zc * pred$se.fit,
                     hi = as.numeric(pred$fit) + zc * pred$se.fit)
  } else {
    ## predict.glm(se.fit = TRUE) returns se.fit on the LINK scale even when
    ## type = "response". Build the interval on the link scale, then invert.
    pred <- stats::predict(fit, newdata = gx, se.fit = TRUE, type = "link")
    inv <- stats::family(fit)$linkinv
    eta <- as.numeric(pred$fit)
    pc <- data.frame(x = grid$x, fit = inv(eta),
                     lo = inv(eta - zc * pred$se.fit),
                     hi = inv(eta + zc * pred$se.fit))
  }
  plt <- ggplot2::ggplot(pc, ggplot2::aes(x = x, y = fit)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2) +
    ggplot2::scale_x_continuous(breaks = .nice_breaks(pc$x, 6)) +
    ggplot2::labs(x = x, y = paste("predicted", y)) +
    mc_theme()
  td <- broom::tidy(fit, conf.int = TRUE, conf.level = conf,
                    exponentiate = !is_lm)
  ## Formal non-linearity test: linear term vs the full spline basis.
  nl_p <- tryCatch({
    f0 <- stats::as.formula("y ~ rcs_lin")
    m0 <- if (is_lm) stats::lm(f0, data = df) else
      stats::glm(f0, data = df, family = family)
    if (is_lm) {
      a <- stats::anova(m0, fit)
      as.numeric(a[2, "Pr(>F)"])
    } else {
      a <- stats::anova(m0, fit, test = "Chisq")
      as.numeric(a[2, "Pr(>Chi)"])
    }
  }, error = function(e) NA_real_)
  list(model = fit, coefficients = td, plot = plt,
       knots = attr(bx, "knots"), curve = pc,
       nonlinearity_p = nl_p,
       checks = .mc_check_row("Non-linearity", nl_p,
                              "spline vs linear p < alpha",
                              if (is.na(nl_p)) NA else nl_p < 0.05,
                              "likelihood-ratio / F test against a linear term"))
}

#' Quantile Regression
#'
#' Corresponds to "Principles and Application Scenarios of Quantile
#' Regression" Built on
#' \code{\link[quantreg]{rq}}, it returns coefficients for multiple quantiles
#' At once.
#'
#' @details
#' Purpose: Use quantile regression to describe how a predictor affects different
#' Quantiles of the outcome (e.g. median, upper/lower quartiles); especially useful
#' Under heteroscedasticity or when focusing on extreme quantiles (e.g. low income,
#' High death risk).
#' Method: Based on quantreg's \code{rq()}, minimizing weighted absolute residuals;
#' It can estimate coefficients for multiple quantiles at once.
#' Arguments: \code{data} is a data frame; \code{formula} is of the form \code{y ~ x};
#' \code{taus} is the vector of quantiles to estimate, default \code{c(0.25, 0.5, 0.75)}.
#' Results: Returns a data frame of coefficients at each quantile (including intercept
#' And predictor slope), allowing comparison of slopes across quantiles.
#' Notes: If the predictor slope differs markedly across quantiles (e.g. under
#' Heteroscedasticity), the conditional distribution shape changes with x -- information
#' Invisible to ordinary mean regression (OLS). Quantile regression is robust to
#' Outliers. Depends on quantreg (Suggests).
#' @param data A `data.frame`.
#' @param formula A regression formula.
#' @param taus Vector of quantile points; default `c(0.25, 0.5, 0.75)`.
#' @return A `data.frame` with columns `term`, `estimate`, `std.error`, `statistic`,
#' `p.value`, `tau` (stacked over taus).
#' @param conf Confidence level for coefficient intervals; default `0.95`.
#' @param se Standard-error method passed to `quantreg::summary.rq()`; `"rank"` (default), `"nid"` fallback.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' \donttest{
#' mc_quantile_reg(d, sbp_week12 ~ sbp_baseline + age,
#'                 taus = c(0.25, 0.5, 0.75))
#' }
#' @export
mc_quantile_reg <- function(data, formula, taus = c(0.25, 0.5, 0.75),
                            conf = 0.95, se = "nid") {
 if (!requireNamespace("quantreg", quietly = TRUE)) {
 stop("`mc_quantile_reg` requires the quantreg package; please install.packages('quantreg').")
 }
 zc <- stats::qnorm(1 - (1 - conf) / 2)
 rows <- list()
 for (t in taus) {
 fit <- quantreg::rq(formula, data = data, tau = t)
 ## "rank" SE fails on non-unique solutions; try robust alternatives in order.
 sm <- NULL
 for (candidate in unique(c(se, "nid", "iid", "ker"))) {
 sm <- suppressWarnings(tryCatch(summary(fit, se = candidate),
 error = function(e) NULL))
 if (!is.null(sm) && !is.null(sm$coefficients) &&
 ncol(sm$coefficients) >= 4L) break
 sm <- NULL
 }
 if (!is.null(sm)) {
 co <- as.data.frame(sm$coefficients)
 est <- co[, 1]; serr <- co[, 2]; tv <- co[, 3]; pv <- co[, 4]
 } else {
 est <- stats::coef(fit)
 serr <- rep(NA_real_, length(est))
 tv <- rep(NA_real_, length(est)); pv <- rep(NA_real_, length(est))
 }
 rows[[length(rows) + 1L]] <- data.frame(
 term = names(stats::coef(fit)), estimate = as.numeric(est),
 std_error = as.numeric(serr), statistic = as.numeric(tv),
 p.value = as.numeric(pv),
 conf.low = as.numeric(est) - zc * as.numeric(serr),
 conf.high = as.numeric(est) + zc * as.numeric(serr),
 tau = rep(t, length(est)), stringsAsFactors = FALSE)
 }
 out <- do.call(rbind, rows)
 rownames(out) <- NULL
 ## Approximate across-quantile slope test (independent-SE approximation).
 terms <- setdiff(unique(out$term), "(Intercept)")
 if (length(taus) >= 2L && length(terms)) {
 i1 <- which(out$tau == min(taus)); i2 <- which(out$tau == max(taus))
 tr <- do.call(rbind, lapply(terms, function(v) {
 a <- out[i1, ][out$term[i1] == v, ][1, ]
 b <- out[i2, ][out$term[i2] == v, ][1, ]
 se2 <- sqrt(a$std_error^2 + b$std_error^2)
 z <- (b$estimate - a$estimate) / se2
 data.frame(term = v, tau_low = min(taus), tau_high = max(taus),
 estimate_low = a$estimate, estimate_high = b$estimate,
 difference = b$estimate - a$estimate, std_error = se2,
 z = z, p.value = 2 * stats::pnorm(-abs(z)),
 stringsAsFactors = FALSE)
 }))
 } else tr <- NULL
 attr(out, "trend") <- tr
 out
}

#' Principal Components Regression
#'
#' The whole pipeline is *self-implemented*:
#' PCA via \code{stats::prcomp}, an in-house k-fold CV to pick the number of
#' Components, and a back-transform of coefficients to the original scale.
#'
#' @details
#' Purpose: When predictors are highly correlated (multicollinear), use principal
#' Components regression for stable estimation.
#' Method (self-implemented): (1) standardize predictors (center + scale); (2) PCA
#' Yields uncorrelated principal components; (3) if the number of components is not
#' Given, an in-house k-fold CV selects the optimum by prediction MSE; (4) fit linear
#' Regression on the selected components; (5) back-transform component coefficients to
#' The original variable scale (using the rotation matrix and standardization parameters).
#' Arguments: \code{data} is a data frame; \code{formula} is of the form
#' \code{y ~ x1 + x2 +...}; \code{ncomp} is the number of components
#' (\code{NULL} = auto-selected by CV); \code{k} is the number of CV folds (default 5);
#' \code{seed} is the CV random seed.
#' Results: Returns a list -- \code{best_ncomp} (chosen components), \code{cv_table}
#' (CV-MSE for each candidate), \code{coef_original} (original-scale coefficients +
#' Intercept), \code{performance} (R-squared / adjusted R-squared / RMSE),
#' \code{plot} (CV-MSE vs number of components), \code{fit} (final lm object).
#' Notes: PCR standardizes predictors, so coefficients relate to the standardized
#' Scale; categorical variables should be dummy-coded first. This function depends
#' Only on base R, no extra package needed.
#' @param data A `data.frame`.
#' @param formula A model `formula`, e.g. `y ~ x1 + x2`.
#' @param ncomp Number of principal components; `NULL` = auto-select by CV.
#' @param k CV folds; default `5`.
#' @param seed Random seed for CV.
#' @return A list with chosen components, CV table, original-scale coefficients,
#' Performance and the CV plot.
#' @param conf Confidence level for original-scale coefficient intervals; default `0.95`.
#' @param validation `"none"` (default) or `"holdout"` for an untouched test set.
#' @param prop Training proportion for holdout validation; default `0.7`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' pcr <- mc_pca_regression(d,
#'   sbp_week12 ~ age + bmi + sbp_baseline + ldl + glucose, k = 5, seed = 1)
#' pcr$best_ncomp
#' pcr$coef_original
#' @export
mc_pca_regression <- function(data, formula, ncomp = NULL, k = 5L,
                              seed = NULL, conf = 0.95,
                              validation = c("none", "holdout"),
                              prop = 0.7) {
 if (!inherits(data, "data.frame") || !inherits(formula, "formula"))
 stop("`data` must be a data.frame and `formula` a formula.")
 validation <- match.arg(validation)
 mf <- stats::model.frame(formula, data = data)
 y <- as.numeric(stats::model.response(mf))
 X <- as.matrix(mf[, -1, drop = FALSE])
 if (!is.numeric(X)) stop("All predictors must be numeric for PCR.")
 p <- ncol(X)
 Xs <- scale(X)
 mu <- attr(Xs, "scaled:center")
 sg <- attr(Xs, "scaled:scale")
 pca <- stats::prcomp(Xs, center = FALSE, scale. = FALSE)
 Z <- pca$x
 maxc <- min(p, length(y) - 2L)
 cv_mse <- NULL
 if (is.null(ncomp)) {
 if (!is.null(seed)) set.seed(seed)
 idx <- sample(rep(1:k, length.out = length(y)))
 cv_mse <- numeric(maxc)
 for (c in 1:maxc) {
 errs <- numeric(length(y))
 for (fold in 1:k) {
 te <- which(idx == fold)
 tr <- setdiff(seq_along(y), te)
 m <- stats::lm(y[tr] ~ Z[tr, 1:c])
 pr <- as.numeric(cbind(1, Z[te, 1:c]) %*% coef(m))
 errs[te] <- (y[te] - pr)^2
 }
 cv_mse[c] <- mean(errs, na.rm = TRUE)
 }
 ncomp <- which.min(cv_mse)
 } else {
 ncomp <- min(as.integer(ncomp), maxc)
 }
 fit <- stats::lm(y ~ Z[, 1:ncomp])
 beta <- as.numeric(coef(fit))
 rot <- pca$rotation[, 1:ncomp, drop = FALSE]
 ## Back-transform PC coefficients: rotation (p x ncomp) times PC slopes.
 gamma <- as.numeric(rot %*% beta[-1])
 coef_orig <- gamma / sg
 intercept <- beta[1] - sum(mu * coef_orig)
 ## Analytic standard errors on the original scale: rotate the PC
 ## coefficient covariance back through the loading matrix.
 zc <- stats::qnorm(1 - (1 - conf) / 2)
 V <- tryCatch(stats::vcov(fit)[-1, -1, drop = FALSE],
 error = function(e) NULL)
 se_orig <- if (!is.null(V))
 sqrt(pmax(0, diag(rot %*% V %*% t(rot)))) / sg else
 rep(NA_real_, length(coef_orig))
 tstat <- coef_orig / se_orig
 coef_original <- data.frame(
 term = c("(Intercept)", colnames(X)),
 estimate = c(intercept, coef_orig),
 std.error = c(NA_real_, se_orig),
 conf.low = c(NA_real_, coef_orig - zc * se_orig),
 conf.high = c(NA_real_, coef_orig + zc * se_orig),
 p.value = c(NA_real_, 2 * stats::pnorm(-abs(tstat))),
 stringsAsFactors = FALSE)
 pred <- as.numeric(cbind(1, Z[, 1:ncomp]) %*% coef(fit))
 perf <- data.frame(
 metric = c("r_squared", "adj_r_squared", "rmse"),
 value = c(summary(fit)$r.squared, summary(fit)$adj.r.squared,
 sqrt(mean((y - pred)^2))), stringsAsFactors = FALSE)
 ## Honest external validation on an untouched holdout set.
 external <- NULL
 if (validation == "holdout" && length(y) >= 20L) {
 if (!is.null(seed)) set.seed(seed + 1L)
 n <- length(y)
 tr <- sample.int(n, max(10L, floor(prop * n)))
 te <- setdiff(seq_len(n), tr)
 Xtr <- scale(X[tr, , drop = FALSE])
 mu_tr <- attr(Xtr, "scaled:center"); sg_tr <- attr(Xtr, "scaled:scale")
 Xte <- sweep(sweep(X[te, , drop = FALSE], 2L, mu_tr, "-"),
 2L, sg_tr, "/")
 pca_tr <- stats::prcomp(Xtr, center = FALSE, scale. = FALSE)
 Rtr <- pca_tr$rotation[, seq_len(min(ncomp, ncol(pca_tr$x))), drop = FALSE]
 Ztr <- pca_tr$x[, seq_len(ncol(Rtr)), drop = FALSE]
 Zte <- as.matrix(Xte) %*% Rtr
 m_tr <- stats::lm(y[tr] ~ Ztr)
 pr <- as.numeric(cbind(1, Zte) %*% coef(m_tr))
 sse <- sum((y[te] - pr)^2, na.rm = TRUE)
 external <- data.frame(
 metric = c("n_train", "n_test", "external_rmse", "external_r2"),
 value = c(length(tr), length(te), sqrt(sse / length(te)),
 1 - sse / sum((y[te] - mean(y[te]))^2)),
 stringsAsFactors = FALSE)
 perf <- rbind(perf, external)
 }
 cv_table <- if (is.null(cv_mse)) NULL else
 data.frame(ncomp = seq_along(cv_mse), cv_mse = cv_mse,
 stringsAsFactors = FALSE)
 plot <- NULL
 if (!is.null(cv_table)) {
 plot <- ggplot2::ggplot(cv_table, ggplot2::aes(ncomp, cv_mse)) +
 ggplot2::geom_line(color = "steelblue") +
 ggplot2::geom_point(size = 2, color = "steelblue") +
 ggplot2::geom_vline(xintercept = ncomp, linetype = "dashed", color = "red") +
 ggplot2::labs(x = "Number of principal components",
 y = "Cross-validated MSE") + mc_theme()
 }
 list(best_ncomp = ncomp, cv_table = cv_table,
 coef_original = coef_original, performance = perf,
 external = external, plot = plot, fit = fit)
}

#' Neural-network regression
#'
#' A single-hidden-layer
#' Feed-forward network for regression via \code{nnet::nnet()} (with `linout = TRUE`),
#' Plus a predicted-vs-observed plot and RMSE.
#'
#' @details
#' Purpose: Use a neural network to fit the non-linear relationship between predictors
#' And a continuous outcome.
#' Dependencies: Requires the nnet package (in Suggests); a clear message is given if
#' Not installed.
#' Method: \code{nnet::nnet(y ~., linout = TRUE, size = hidden)}; the number of hidden
#' Nodes is controlled by \code{hidden}; data are auto-standardized to speed up convergence.
#' Arguments: \code{data} is a data frame; \code{formula} is of the form
#' \code{y ~ x1 + x2 +...}; \code{hidden} is the number of hidden nodes (default 5);
#' \code{maxit} is the maximum iterations (default 200); \code{seed} is the random seed.
#' Results: Returns a list -- \code{fit} (nnet object), \code{predictions} (fitted
#' Values), \code{rmse}, \code{plot} (predicted vs observed scatter).
#' Notes: Neural nets are scale-sensitive; this function standardizes predictors
#' Internally. Results are sensitive to random initialization, so set \code{seed} for
#' Reproducibility. Too many nodes overfit; compare structures with
#' \code{mc_model_compare()} or cross-validation.
#' @param data A `data.frame`.
#' @param formula A model `formula`.
#' @param hidden Hidden-layer nodes; default `5`.
#' @param maxit Maximum iterations; default `200`.
#' @param seed Random seed.
#' @return A list with the fitted net, predictions, RMSE and the fit plot.
#' @param validation `"holdout"` (default), `"cv"` or `"none"`. Out-of-sample error replaces the optimistic in-sample RMSE.
#' @param prop Training-set proportion for `validation = "holdout"`; default `0.7`.
#' @param folds Number of folds for `validation = "cv"`; default `5`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' \donttest{
#' nn <- mc_nnet_reg(d,
#'   sbp_week12 ~ age + bmi + sbp_baseline + treatment, hidden = 5, seed = 7)
#' nn$validation_rmse
#' }
#' @export
mc_nnet_reg <- function(data, formula, hidden = 5, maxit = 200, seed = NULL,
                        validation = c("holdout", "cv", "none"),
                        prop = 0.7, folds = 5L) {
  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("`mc_nnet_reg` requires the nnet package; please install.packages('nnet').")
  }
  validation <- match.arg(validation)
  mf <- stats::model.frame(formula, data = data)
  y <- as.numeric(stats::model.response(mf))
  n <- length(y)
  fit_one <- function(d) suppressWarnings(
    nnet::nnet(formula, data = d, size = hidden, linout = TRUE,
               maxit = maxit, trace = FALSE))
  if (!is.null(seed)) set.seed(seed)
  fit <- tryCatch(fit_one(data), error = function(e) NULL)
  if (is.null(fit)) stop("The neural network failed to converge; reduce `hidden` or increase `maxit`.")
  pred <- as.numeric(stats::predict(fit, newdata = data))
  train_rmse <- sqrt(mean((y - pred)^2, na.rm = TRUE))
  val_rmse <- NA_real_; val_r2 <- NA_real_
  if (validation == "holdout" && n >= 10L) {
    if (!is.null(seed)) set.seed(seed + 1L)
    tr <- sample.int(n, max(5L, floor(prop * n)))
    te <- setdiff(seq_len(n), tr)
    m <- tryCatch(fit_one(data[tr, , drop = FALSE]), error = function(e) NULL)
    if (!is.null(m)) {
      yt <- y[te]
      pt <- as.numeric(stats::predict(m, newdata = data[te, , drop = FALSE]))
      sse <- sum((yt - pt)^2, na.rm = TRUE)
      val_rmse <- sqrt(mean((yt - pt)^2, na.rm = TRUE))
      val_r2 <- 1 - sse / sum((yt - mean(yt))^2)
    }
  } else if (validation == "cv" && n >= 10L) {
    if (!is.null(seed)) set.seed(seed + 1L)
    idx <- sample(rep(seq_len(folds), length.out = n))
    errs <- rep(NA_real_, n)
    for (k in seq_len(folds)) {
      te <- which(idx == k); tr <- setdiff(seq_len(n), te)
      m <- tryCatch(fit_one(data[tr, , drop = FALSE]), error = function(e) NULL)
      if (!is.null(m))
        errs[te] <- y[te] - as.numeric(
          stats::predict(m, newdata = data[te, , drop = FALSE]))
    }
    mse <- mean(errs^2, na.rm = TRUE)
    val_rmse <- sqrt(mse)
    val_r2 <- 1 - mse / stats::var(y)
  }
  ## RMSE now reports out-of-sample error whenever validation is possible.
  rmse <- if (is.na(val_rmse)) train_rmse else val_rmse
  performance <- data.frame(
    metric = c("train_rmse", "validation_rmse", "validation_r2", "n"),
    value = c(train_rmse, val_rmse, val_r2, n),
    stringsAsFactors = FALSE)
  plot <- ggplot2::ggplot(data.frame(obs = y, pred = pred),
                          ggplot2::aes(obs, pred)) +
    ggplot2::geom_point(size = 1.3, alpha = 0.6) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red",
                         linetype = "dashed") +
    ggplot2::labs(x = "Observed", y = "Predicted") + mc_theme()
  list(fit = fit, predictions = pred, rmse = rmse,
       performance = performance, validation = validation, plot = plot)
}
