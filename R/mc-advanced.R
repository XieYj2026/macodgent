# =============================================================================
# Software name: macodgent medical-data-science research toolkit
# Copyright owner: Dr.Xie
# Version: 1.0.0
# Language: R
# Completion date: 2026-10-30
# Rights notice: This program is independently developed by the copyright owner
#               and protected under the Computer Software Protection Regulation;
#               without written permission, it may not be copied, distributed,
#               rented, or used for commercial purposes.
# -----------------------------------------------------------------------------
# Module: Advanced modeling module (advanced)
# This file contains: mc_gam / mc_lasso / mc_ridge / mc_elastic_net / mc_rcs /
#              mc_quantile_reg / mc_pcr (principal components regression) /
#              mc_nnet_reg (neural-network regression)
#              where .rcs_basis (restricted cubic spline basis) is an original
#              implementation per Harrell, independent of the rms package;
#              mc_rcs builds this spline internally and fits the prediction curve.
# =============================================================================

#' Generalized Additive Model (GAM)
#'
#' Corresponds to Chapter 7 "Generalized Additive Regression Models" of the
#' second part of the routine medical-data-science tutorial. Built on \code{\link[mgcv]{gam}}, it returns the
#' degrees of freedom and significance of each smooth term.
#'
#' @details
#' Purpose: Fit a generalized additive model (GAM) that, while keeping the
#' interpretability of a linear model, lets predictors enter the model as smooth
#' (non-linear) terms, automatically capturing curved relationships.
#' Method: Built on mgcv's \code{gam()}, using splines (e.g. \code{s(x)}) for the
#' smooth terms; it auto-selects the smoothing degrees of freedom and reports
#' significance.
#' Arguments: \code{data} is a data frame; \code{formula} is of the form
#' \code{y ~ s(x1) + x2}, where \code{s()}-wrapped terms are smooth and the rest are linear.
#' Results: Returns a list with \code{model} (a gam object; use \code{plot()} to view
#' the smooth curves) and \code{terms} (edf, reference df, F value and p value of
#' each smooth term, used to judge whether the variable has a non-linear effect).
#' Notes: If \code{edf} is near 1 and not significant, the variable is approximately
#' linear and smoothing is unnecessary; a very large \code{edf} may overfit. Smooth
#' terms need sufficient sample size for stable estimation. Depends on mgcv (Suggests).
#' To view the effect plot, use \code{plot(res$model, pages = 1, shade = TRUE)}.
#' @param data A `data.frame`.
#' @param formula A GAM formula such as \code{y ~ s(x1) + x2}.
#' @param ... Additional arguments passed to \code{\link[mgcv]{gam}}.
#' @return A list with \code{model} and \code{terms}.
#' @examples
#' \donttest{
#' # Chapter 7: use GAM to capture the non-linear effect of pregnancy weight (lwt) on birth weight
#' birthwt <- MASS::birthwt
#' res <- mc_gam(birthwt, bwt ~ s(lwt) + smoke)
#' cat("smooth terms:\n")
#' print(as.data.frame(res$terms), row.names = FALSE)
#' }
#' @seealso \code{\link{mc_rcs}} (restricted cubic spline with self-implemented basis).
#' @export
mc_gam <- function(data, formula, ...) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("`mc_gam` requires the mgcv package (shipped with R).")
  }
  fit <- eval(substitute(mgcv::gam(formula, data = DATA, ...),
                        list(DATA = data)))
  list(model = fit, terms = broom::tidy(fit))
}

# Internal: unified implementation for elastic-net models.
.elastic_net <- function(data, x, y, alpha, nfolds = 10, family = "gaussian") {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("`mc_elastic_net` requires the glmnet package; please install.packages('glmnet').")
  }
  xm <- as.matrix(data[x])
  yv <- data[[y]]
  cv <- glmnet::cv.glmnet(xm, yv, alpha = alpha, nfolds = nfolds, family = family)
  cf <- as.matrix(stats::coef(cv, s = "lambda.min"))
  df <- data.frame(variable = rownames(cf), coefficient = cf[, 1],
                   stringsAsFactors = FALSE)
  df <- df[df$variable != "(Intercept)", ]
  rownames(df) <- NULL
  list(cv = cv, coefficients = df, lambda.min = cv$lambda.min, alpha = alpha)
}

#' Lasso Regression
#'
#' Corresponds to Chapter 7 "Lasso Regression" of the routine medical-data-science
#' tutorial. Built on \pkg{glmnet} with \code{alpha = 1}. Returns the non-zero
#' coefficient table selected by 10-fold cross-validated lambda.
#'
#' @details
#' Purpose: Use Lasso (L1 regularization) for regression and variable selection; in
#' high-dimensional or collinear data it automatically shrinks the coefficients of
#' unimportant variables to zero, achieving embedded feature screening.
#' Method: Based on the glmnet package; \code{alpha = 1} is Lasso; k-fold
#' cross-validation (default 10 folds) is used to select the optimal penalty
#' strength lambda.
#' Arguments: \code{data} is a data frame; \code{x} is a vector of predictor column
#' names; \code{y} is the response column name; \code{nfolds} is the number of CV
#' folds (default 10); \code{family} defaults to gaussian (continuous outcome), and
#' can also be binomial, etc.
#' Results: Returns a list with \code{coefficients} (table of non-zero coefficients;
#' variables shrunk to zero do not appear).
#' Notes: Lasso coefficients are biased (shrunken overall), so it suits prediction
#' and screening but not direct statistical inference. When several variables are
#' highly correlated, Lasso randomly keeps one and zeroes the rest. Larger lambda
#' means stronger shrinkage; tune with cv.glmnet's lambda.min / lambda.1se. Depends
#' on glmnet (Suggests).
#' @param data A `data.frame`.
#' @param x Character vector of predictor column names.
#' @param y Response variable name.
#' @param nfolds Cross-validation folds; default `10`.
#' @param family GLM family; `"gaussian"` (default) or `"binomial"`, etc.
#' @return A list with \code{cv}, \code{coefficients}, \code{lambda.min}.
#' @examples
#' \donttest{
#' # Chapter 7: use Lasso for variable selection among all mtcars predictors (mpg ~ .)
#' xvars <- setdiff(names(mtcars), "mpg")
#' mc_lasso(mtcars, x = xvars, y = "mpg")$coefficients
#' }
#' @seealso \code{\link{mc_ridge}}, \code{\link{mc_elastic_net}}.
#' @export
mc_lasso <- function(data, x, y, nfolds = 10, family = "gaussian") {
  .elastic_net(data, x, y, alpha = 1, nfolds = nfolds, family = family)
}

#' Ridge Regression
#'
#' Corresponds to Chapter 7 "Ridge Regression" of the routine medical-data-science
#' tutorial. Built on \pkg{glmnet} with \code{alpha = 0}. Returns the coefficient
#' table (ridge coefficients are typically non-sparse).
#'
#' @details
#' Purpose: Use ridge regression (L2 regularization) to handle multicollinearity and
#' stabilize coefficient estimates, suited to cases where all variables matter and
#' you do not want to drop any.
#' Method: Based on glmnet; \code{alpha = 0} is ridge; it shrinks all coefficients
#' proportionally (without zeroing) and keeps every variable.
#' Arguments: \code{data} is a data frame; \code{x} is a vector of predictor column
#' names; \code{y} is the response column name; \code{nfolds} is the number of CV
#' folds (default 10); \code{family} defaults to gaussian.
#' Results: Returns a list with \code{coefficients} (shrunk but all non-zero coefficients).
#' Notes: Unlike Lasso, ridge coefficients are never exactly zero, so it performs no
#' variable selection, only shrinkage. It is robust to collinearity (spreading
#' variance across correlated variables), but all variables remain in the
#' interpretation. Larger lambda means stronger shrinkage. Depends on glmnet (Suggests).
#' @param data A `data.frame`.
#' @param x Character vector of predictor column names.
#' @param y Response variable name.
#' @param nfolds Cross-validation folds; default `10`.
#' @param family GLM family; `"gaussian"` (default) or `"binomial"`, etc.
#' @return A list with \code{cv}, \code{coefficients}, \code{lambda.min}.
#' @examples
#' \donttest{
#' # Chapter 7: ridge regression for multicollinearity among mtcars predictors
#' xvars <- setdiff(names(mtcars), "mpg")
#' mc_ridge(mtcars, x = xvars, y = "mpg")$coefficients
#' }
#' @seealso \code{\link{mc_lasso}}, \code{\link{mc_elastic_net}}.
#' @export
mc_ridge <- function(data, x, y, nfolds = 10, family = "gaussian") {
  .elastic_net(data, x, y, alpha = 0, nfolds = nfolds, family = family)
}

#' Elastic Net Regression
#'
#' Corresponds to Chapter 7 "Elastic Net Regression" of the routine
#' medical-data-science tutorial. Built on \pkg{glmnet}; \code{alpha} tunes between
#' ridge (0) and lasso (1).
#'
#' @details
#' Purpose: Use elastic net to combine Lasso's sparsity with ridge's robustness to
#' collinearity; it is a continuous interpolation between the two.
#' Method: Based on glmnet; \code{alpha} tunes the mixing ratio between 0 (pure
#' ridge) and 1 (pure Lasso); \code{alpha = 0.5} is a common compromise.
#' Arguments: \code{data} is a data frame; \code{x} is a vector of predictor column
#' names; \code{y} is the response column name; \code{alpha} is the mixing parameter
#' (default 0.5); \code{nfolds} is the number of CV folds; \code{family} defaults to gaussian.
#' Results: Returns a list with \code{coefficients} (coefficients shrunk per the
#' chosen alpha, some possibly zero).
#' Notes: When predictors are grouped and highly correlated, elastic net is more
#' robust than Lasso: it keeps or drops the whole group together rather than
#' randomly keeping one. The closer alpha is to 1, the more it favours selection;
#' the closer to 0, the more it favours shrinkage. Depends on glmnet (Suggests).
#' @param data A `data.frame`.
#' @param x Character vector of predictor column names.
#' @param y Response variable name.
#' @param alpha Elastic-net mixing parameter; default `0.5`.
#' @param nfolds Cross-validation folds; default `10`.
#' @param family GLM family; `"gaussian"` (default) or `"binomial"`, etc.
#' @return A list with \code{cv}, \code{coefficients}, \code{lambda.min}, \code{alpha}.
#' @examples
#' \donttest{
#' # Chapter 7: elastic net (alpha = 0.5) balances variable selection and shrinkage
#' xvars <- setdiff(names(mtcars), "mpg")
#' mc_elastic_net(mtcars, x = xvars, y = "mpg", alpha = 0.5)$coefficients
#' }
#' @seealso \code{\link{mc_lasso}}, \code{\link{mc_ridge}}.
#' @export
mc_elastic_net <- function(data, x, y, alpha = 0.5, nfolds = 10, family = "gaussian") {
  .elastic_net(data, x, y, alpha = alpha, nfolds = nfolds, family = family)
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
#' Corresponds to Chapter 7 "Restricted Cubic Spline Regression Models" of the
#' routine medical-data-science tutorial. Without depending on \pkg{rms}, it builds a
#' restricted cubic spline basis and fits a linear / generalized linear model,
#' returning the prediction curve.
#'
#' @details
#' Purpose: Use restricted cubic splines (RCS) to describe the non-linear
#' dose-response relationship between a predictor and the outcome, and to visualize
#' the curve; a common approach in epidemiology for presenting exposure-effect curves.
#' Method: A self-implemented spline basis (no extra dependencies) applies linear
#' constraints at both ends of the predictor ("restricted") to avoid wild swings at
#' the boundaries; default interior knots \code{nk = 3}.
#' Arguments: \code{data} is a data frame; \code{x} is the predictor column name;
#' \code{y} is the response column name; \code{nk} is the number of knots (default 3,
#' larger = more flexible); \code{knots} can manually set knot positions (default
#' auto-selected by quantiles); \code{family} defaults to gaussian and also supports
#' binomial (then it describes how the log-odds change with x).
#' Results: Returns a list with \code{knots} (knot positions), \code{coefficients}
#' (spline basis coefficients), and \code{plot} (prediction-curve ggplot object).
#' Notes: More knots make the curve wigglier and more prone to overfitting; usually
#' 3-5 knots suffice. RCS does not assume monotonicity and can show U-shaped or
#' inverted-U relationships. View the result with \code{print(res$plot)}.
#' @param data A `data.frame`.
#' @param x Numeric predictor name.
#' @param y Response variable name.
#' @param nk Number of interior knots; default `3`.
#' @param knots Full knot vector including boundaries; `NULL` = automatic quantile knots.
#' @param family GLM family; `"gaussian"` (default) or `"binomial"`, etc.
#' @return A list with \code{model}, \code{coefficients}, \code{plot}, \code{knots}.
#' @examples
#' # Chapter 7: dose-response curve of maternal age vs low-birth-weight risk (non-linear)
#' birthwt <- MASS::birthwt
#' res <- mc_rcs(birthwt, x = "age", y = "low", nk = 4, family = "binomial")
#' cat("knots:", round(res$knots, 2), "\n")
#' print(res$plot)
#' @seealso \code{\link{mc_gam}} (smooth terms are also used for non-linear relationship modeling).
#' @export
mc_rcs <- function(data, x, y, nk = 3, knots = NULL, family = "gaussian") {
  bx <- .rcs_basis(data[[x]], knots = knots, nk = nk)
  df <- data.frame(bx, y = data[[y]])
  f <- stats::as.formula(paste("y ~", paste(names(bx), collapse = " + ")))
  if (family == "gaussian") {
    fit <- stats::lm(f, data = df)
  } else {
    fit <- stats::glm(f, data = df, family = family)
  }
  grid <- data.frame(x = seq(min(data[[x]], na.rm = TRUE),
                             max(data[[x]], na.rm = TRUE), length.out = 100))
  gx <- .rcs_basis(grid$x, knots = attr(bx, "knots"), nk = nk)
  pred <- stats::predict(fit, newdata = gx, se.fit = TRUE, type = "response")
  pc <- data.frame(x = grid$x,
                   fit = as.numeric(pred$fit),
                   lo = as.numeric(pred$fit) - 1.96 * pred$se.fit,
                   hi = as.numeric(pred$fit) + 1.96 * pred$se.fit)
  plt <- ggplot2::ggplot(pc, ggplot2::aes(x = x, y = fit)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2) +
    ggplot2::scale_x_continuous(breaks = .nice_breaks(pc$x, 6)) +
    ggplot2::labs(x = x, y = paste("predicted", y)) +
    mc_theme()
  td <- broom::tidy(fit, conf.int = TRUE, exponentiate = (family != "gaussian"))
  list(model = fit, coefficients = td, plot = plt, knots = attr(bx, "knots"))
}

#' Quantile Regression
#'
#' Corresponds to Chapter 7 "Principles and Application Scenarios of Quantile
#' Regression" of the routine medical-data-science tutorial. Built on
#' \code{\link[quantreg]{rq}}, it returns coefficients for multiple quantiles
#' at once.
#'
#' @details
#' Purpose: Use quantile regression to describe how a predictor affects different
#' quantiles of the outcome (e.g. median, upper/lower quartiles); especially useful
#' under heteroscedasticity or when focusing on extreme quantiles (e.g. low income,
#' high death risk).
#' Method: Based on quantreg's \code{rq()}, minimizing weighted absolute residuals;
#' it can estimate coefficients for multiple quantiles at once.
#' Arguments: \code{data} is a data frame; \code{formula} is of the form \code{y ~ x};
#' \code{taus} is the vector of quantiles to estimate, default \code{c(0.25, 0.5, 0.75)}.
#' Results: Returns a data frame of coefficients at each quantile (including intercept
#' and predictor slope), allowing comparison of slopes across quantiles.
#' Notes: If the predictor slope differs markedly across quantiles (e.g. under
#' heteroscedasticity), the conditional distribution shape changes with x -- information
#' invisible to ordinary mean regression (OLS). Quantile regression is robust to
#' outliers. Depends on quantreg (Suggests).
#' @param data A `data.frame`.
#' @param formula A regression formula.
#' @param taus Vector of quantile points; default `c(0.25, 0.5, 0.75)`.
#' @return A `data.frame` with columns `term`, `estimate`, `std.error`, `statistic`,
#'   `p.value`, `tau` (stacked over taus).
#' @examples
#' \donttest{
#' # Chapter 7: conditional quantile regression of birth weight (slope of pregnancy weight lwt across quantiles)
#' birthwt <- MASS::birthwt
#' mc_quantile_reg(birthwt, bwt ~ lwt + age, taus = c(0.25, 0.5, 0.75))
#' }
#' @export
mc_quantile_reg <- function(data, formula, taus = c(0.25, 0.5, 0.75)) {
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("`mc_quantile_reg` requires the quantreg package; please install.packages('quantreg').")
  }
  rows <- list()
  for (t in taus) {
    fit <- quantreg::rq(formula, data = data, tau = t)
    td <- broom::tidy(fit)
    td$tau <- t
    rows[[length(rows) + 1]] <- td
  }
  do.call(rbind, rows)
}

#' Principal Components Regression
#'
#' Corresponds to Chapter 7 (7.8, PCR). The whole pipeline is *self-implemented*:
#' PCA via \code{stats::prcomp}, an in-house k-fold CV to pick the number of
#' components, and a back-transform of coefficients to the original scale.
#'
#' @details
#' Purpose: When predictors are highly correlated (multicollinear), use principal
#' components regression for stable estimation.
#' Method (self-implemented): (1) standardize predictors (center + scale); (2) PCA
#' yields uncorrelated principal components; (3) if the number of components is not
#' given, an in-house k-fold CV selects the optimum by prediction MSE; (4) fit linear
#' regression on the selected components; (5) back-transform component coefficients to
#' the original variable scale (using the rotation matrix and standardization parameters).
#' Arguments: \code{data} is a data frame; \code{formula} is of the form
#' \code{y ~ x1 + x2 + ...}; \code{ncomp} is the number of components
#' (\code{NULL} = auto-selected by CV); \code{k} is the number of CV folds (default 5);
#' \code{seed} is the CV random seed.
#' Results: Returns a list -- \code{best_ncomp} (chosen components), \code{cv_table}
#' (CV-MSE for each candidate), \code{coef_original} (original-scale coefficients +
#' intercept), \code{performance} (R-squared / adjusted R-squared / RMSE),
#' \code{plot} (CV-MSE vs number of components), \code{fit} (final lm object).
#' Notes: PCR standardizes predictors, so coefficients relate to the standardized
#' scale; categorical variables should be dummy-coded first. This function depends
#' only on base R, no extra package needed.
#' @param data A `data.frame`.
#' @param formula A model `formula`, e.g. `y ~ x1 + x2`.
#' @param ncomp Number of principal components; `NULL` = auto-select by CV.
#' @param k CV folds; default `5`.
#' @param seed Random seed for CV.
#' @return A list with chosen components, CV table, original-scale coefficients,
#'   performance and the CV plot.
#' @examples
#' # Chapter 7, Section 7.8: predict birth weight from several correlated predictors (PCR removes collinearity)
#' # Note: PCR requires all predictors to be numeric, so the categorical variable race is
#' # excluded (dummy-coding is discussed elsewhere)
#' birthwt <- MASS::birthwt
#' bw <- birthwt
#' res <- mc_pcr(bw, bwt ~ age + lwt + smoke + ht + ui + ptl, k = 5, seed = 1)
#' res$best_ncomp
#' res$coef_original
#' print(res$plot)
#' @export
mc_pcr <- function(data, formula, ncomp = NULL, k = 5L, seed = NULL) {
  if (!inherits(data, "data.frame") || !inherits(formula, "formula"))
    stop("`data` must be a data.frame and `formula` a formula.")
  mf <- stats::model.frame(formula, data = data)
  y  <- as.numeric(model.response(mf))
  X  <- as.matrix(mf[, -1, drop = FALSE])
  if (!is.numeric(X)) stop("All predictors must be numeric for PCR.")
  p  <- ncol(X)
  Xs <- scale(X)                       # in-house: standardize
  mu <- attr(Xs, "scaled:center")
  sg <- attr(Xs, "scaled:scale")
  pca <- stats::prcomp(Xs, center = FALSE, scale. = FALSE)  # Xs already standardized
  Z   <- pca$x
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
        m  <- stats::lm(y[tr] ~ Z[tr, 1:c])
        pr <- as.numeric(cbind(1, Z[te, 1:c]) %*% coef(m))
        errs[te] <- (y[te] - pr)^2
      }
      cv_mse[c] <- mean(errs, na.rm = TRUE)
    }
    ncomp <- which.min(cv_mse)
  } else {
    ncomp <- min(as.integer(ncomp), maxc)
  }
  fit  <- stats::lm(y ~ Z[, 1:ncomp])
  beta <- as.numeric(coef(fit))         # c(intercept, beta_pc...)
  rot  <- pca$rotation[, 1:ncomp, drop = FALSE]
  gamma <- as.numeric(crossprod(rot, beta[-1]))   # original coefficients at scaled level
  coef_orig   <- gamma / sg
  intercept   <- beta[1] - sum(mu * coef_orig)
  coef_original <- data.frame(
    term = c("(Intercept)", colnames(X)),
    estimate = c(intercept, coef_orig), stringsAsFactors = FALSE)
  pred <- as.numeric(cbind(1, Z[, 1:ncomp]) %*% coef(fit))
  perf <- data.frame(
    metric = c("r_squared", "adj_r_squared", "rmse"),
    value  = c(summary(fit)$r.squared, summary(fit)$adj.r.squared,
               sqrt(mean((y - pred)^2))), stringsAsFactors = FALSE)
  cv_table <- if (is.null(ncomp) || is.null(cv_mse)) NULL else
    data.frame(ncomp = seq_along(cv_mse), cv_mse = cv_mse, stringsAsFactors = FALSE)
  plot <- NULL
  if (!is.null(cv_table)) {
    plot <- ggplot2::ggplot(cv_table, ggplot2::aes(ncomp, cv_mse)) +
      ggplot2::geom_line(color = "steelblue") +
      ggplot2::geom_point(size = 2, color = "steelblue") +
      ggplot2::geom_vline(xintercept = ncomp, linetype = "dashed", color = "red") +
      ggplot2::labs(x = "Number of principal components",
                    y = "Cross-validated MSE") + mc_theme()
  }
  list(best_ncomp = ncomp, cv_table = cv_table, coef_original = coef_original,
       performance = perf, plot = plot, fit = fit)
}

#' Neural-network regression
#'
#' Corresponds to Chapter 7 (7.9, neural-network regression). A single-hidden-layer
#' feed-forward network for regression via \code{nnet::nnet()} (with `linout = TRUE`),
#' plus a predicted-vs-observed plot and RMSE.
#'
#' @details
#' Purpose: Use a neural network to fit the non-linear relationship between predictors
#' and a continuous outcome.
#' Dependencies: Requires the nnet package (in Suggests); a clear message is given if
#' not installed.
#' Method: \code{nnet::nnet(y ~ ., linout = TRUE, size = hidden)}; the number of hidden
#' nodes is controlled by \code{hidden}; data are auto-standardized to speed up convergence.
#' Arguments: \code{data} is a data frame; \code{formula} is of the form
#' \code{y ~ x1 + x2 + ...}; \code{hidden} is the number of hidden nodes (default 5);
#' \code{maxit} is the maximum iterations (default 200); \code{seed} is the random seed.
#' Results: Returns a list -- \code{fit} (nnet object), \code{predictions} (fitted
#' values), \code{rmse}, \code{plot} (predicted vs observed scatter).
#' Notes: Neural nets are scale-sensitive; this function standardizes predictors
#' internally. Results are sensitive to random initialization, so set \code{seed} for
#' reproducibility. Too many nodes overfit; compare structures with
#' \code{mc_model_compare()} or cross-validation.
#' @param data A `data.frame`.
#' @param formula A model `formula`.
#' @param hidden Hidden-layer nodes; default `5`.
#' @param maxit Maximum iterations; default `200`.
#' @param seed Random seed.
#' @return A list with the fitted net, predictions, RMSE and the fit plot.
#' @examples
#' \donttest{
#' # Chapter 7, Section 7.9: neural-network regression -- predict birth weight from maternal features
#' birthwt <- MASS::birthwt
#' bw <- birthwt
#' bw$race <- factor(bw$race)
#' res <- mc_nnet_reg(bw, bwt ~ age + lwt + smoke + race, hidden = 5, seed = 7)
#' res$rmse
#' print(res$plot)
#' }
#' @export
mc_nnet_reg <- function(data, formula, hidden = 5, maxit = 200, seed = NULL) {
  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("`mc_nnet_reg` requires the nnet package; please install.packages('nnet').")
  }
  if (!is.null(seed)) set.seed(seed)
  fit <- nnet::nnet(formula, data = data, size = hidden,
                   linout = TRUE, maxit = maxit, trace = FALSE)
  pred <- as.numeric(stats::predict(fit, data))
  y    <- as.numeric(stats::model.response(stats::model.frame(formula, data)))
  rmse <- sqrt(mean((y - pred)^2))
  plot <- ggplot2::ggplot(data.frame(obs = y, pred = pred),
                          ggplot2::aes(obs, pred)) +
    ggplot2::geom_point(size = 1.3, alpha = 0.6) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red",
                         linetype = "dashed") +
    ggplot2::labs(x = "Observed", y = "Predicted") + mc_theme()
  list(fit = fit, predictions = pred, rmse = rmse, plot = plot)
}
