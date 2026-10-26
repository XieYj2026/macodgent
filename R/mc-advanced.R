#' Generalized Additive Model (GAM)
#'
#' Corresponds to Chapter 7 "Generalized Additive Regression Models" of the
#' second part of the column. Built on \code{\link[mgcv]{gam}}, it returns the
#' degrees of freedom and significance of each smooth term.
#'
#' @param data A data.frame.
#' @param formula A GAM formula such as \code{y ~ s(x1) + x2}.
#' @param ... Additional arguments passed to \code{\link[mgcv]{gam}}.
#' @return A list with \code{model} (the fitted gam object) and \code{terms}
#'   (a tidy table of smooth terms).
#' @examples
#' \dontrun{
#' mc_gam(MASS::birthwt, bwt ~ s(age) + smoke)
#' }
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
#' Corresponds to Chapter 7 "Lasso Regression" of the second part of the
#' column. Built on \pkg{glmnet} with \code{alpha = 1}.
#'
#' @param data A data.frame.
#' @param x Character vector of predictor column names.
#' @param y Name of the response variable.
#' @param nfolds Number of cross-validation folds, default 10.
#' @param family "gaussian" (default) or "binomial", etc.
#' @return A list with \code{cv}, \code{coefficients} (non-zero coefficient
#'   table), and \code{lambda.min}.
#' @examples
#' \dontrun{
#' mc_lasso(mtcars, x = c("hp", "wt", "disp"), y = "mpg")
#' }
#' @export
mc_lasso <- function(data, x, y, nfolds = 10, family = "gaussian") {
  .elastic_net(data, x, y, alpha = 1, nfolds = nfolds, family = family)
}

#' Ridge Regression
#'
#' Corresponds to Chapter 7 "Ridge Regression" of the second part of the
#' column. Built on \pkg{glmnet} with \code{alpha = 0}.
#'
#' @param data A data.frame.
#' @param x Character vector of predictor column names.
#' @param y Name of the response variable.
#' @param nfolds Number of cross-validation folds, default 10.
#' @param family "gaussian" (default) or "binomial", etc.
#' @return A list with \code{cv}, \code{coefficients} (coefficient table), and
#'   \code{lambda.min}.
#' @examples
#' \dontrun{
#' mc_ridge(mtcars, x = c("hp", "wt", "disp"), y = "mpg")
#' }
#' @export
mc_ridge <- function(data, x, y, nfolds = 10, family = "gaussian") {
  .elastic_net(data, x, y, alpha = 0, nfolds = nfolds, family = family)
}

#' Elastic Net Regression
#'
#' Corresponds to Chapter 7 "Elastic Net Regression" of the second part of the
#' column. Built on \pkg{glmnet}; \code{alpha} tunes between ridge (0) and
#' lasso (1).
#'
#' @param data A data.frame.
#' @param x Character vector of predictor column names.
#' @param y Name of the response variable.
#' @param alpha Elastic-net mixing parameter, default 0.5.
#' @param nfolds Number of cross-validation folds, default 10.
#' @param family "gaussian" (default) or "binomial", etc.
#' @return A list with \code{cv}, \code{coefficients} (coefficient table),
#'   \code{lambda.min}, and \code{alpha}.
#' @examples
#' \dontrun{
#' mc_elastic_net(mtcars, x = c("hp", "wt", "disp"), y = "mpg", alpha = 0.5)
#' }
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

#' Restricted Cubic Spline Regression (Self-implemented Spline Basis)
#'
#' Corresponds to Chapter 7 "Restricted Cubic Spline Regression Models" of the
#' second part of the column. Without depending on \pkg{rms}, it builds a
#' restricted cubic spline basis and fits a linear / generalized linear model,
#' returning the prediction curve.
#'
#' @param data A data.frame.
#' @param x Numeric predictor name.
#' @param y Response variable name.
#' @param nk Number of interior knots, default 3.
#' @param knots Full knot vector (including boundaries); NULL uses automatic
#'   quantile-based knots.
#' @param family "gaussian" (default) or "binomial", etc.
#' @return A list with \code{model}, \code{coefficients} (coefficient table),
#'   \code{plot} (prediction curve), and \code{knots}.
#' @examples
#' mc_rcs(MASS::birthwt, x = "age", y = "bwt", nk = 3)
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
    ggplot2::labs(x = x, y = paste("predicted", y)) +
    mc_theme()
  td <- broom::tidy(fit, conf.int = TRUE, exponentiate = (family != "gaussian"))
  list(model = fit, coefficients = td, plot = plt, knots = attr(bx, "knots"))
}

#' Quantile Regression
#'
#' Corresponds to Chapter 7 "Principles and Application Scenarios of Quantile
#' Regression" of the second part of the column. Built on
#' \code{\link[quantreg]{rq}}, it returns coefficients for multiple quantiles
#' at once.
#'
#' @param data A data.frame.
#' @param formula A regression formula.
#' @param taus Vector of quantile points, default c(0.25, 0.5, 0.75).
#' @return A data.frame with columns term / estimate / std.error / statistic /
#'   p.value / tau.
#' @examples
#' \dontrun{
#' mc_quantile_reg(MASS::birthwt, bwt ~ age + lwt, taus = c(0.25, 0.5, 0.75))
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
