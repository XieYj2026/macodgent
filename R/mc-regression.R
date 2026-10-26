#' Binary logistic regression
#'
#' Corresponds to Chapter 5: based on \code{glm(family = binomial)}, returns a
#' tidy coefficient table (with OR and its 95% CI).
#'
#' @param data A data.frame.
#' @param formula Model formula, e.g. \code{y ~ x1 + x2}.
#' @param conf Confidence level. Default 0.95.
#' @return A data.frame: term / estimate(OR) / std.error / statistic / p.value /
#'   conf.low / conf.high
#' @examples
#' bw <- MASS::birthwt
#' bw$low <- factor(bw$low)
#' mc_logistic(bw, low ~ age + smoke)
#' @export
mc_logistic <- function(data, formula, conf = 0.95) {
  fit <- eval(substitute(stats::glm(formula, data = DATA,
                                    family = stats::binomial),
                         list(DATA = data)))
  td <- broom::tidy(fit, conf.int = TRUE, conf.level = conf, exponentiate = TRUE)
  td
}

#' Univariate screening then multivariable logistic regression (flagship)
#'
#' Corresponds to Chapter 5 (extended reading): run a univariate logistic
#' regression for each candidate variable, screen significant ones at \code{alpha},
#' then build a multivariable model from the selected variables.
#'
#' @param data A data.frame.
#' @param y Binary outcome variable name (character).
#' @param x_vars Candidate predictor names (character vector).
#' @param alpha Univariate screening significance threshold. Default 0.05.
#' @param conf Multivariable model confidence level. Default 0.95.
#' @return A list: \code{univariate} (univariate table), \code{selected}
#'   (selected variables), \code{multivariable} (multivariable coefficient table).
#' @examples
#' bw <- MASS::birthwt
#' mc_logistic_screen(bw, y = "low",
#'   x_vars = c("age", "lwt", "race", "smoke", "ht", "ui"))
#' @export
mc_logistic_screen <- function(data, y, x_vars, alpha = 0.05, conf = 0.95) {
  univ <- list()
  sig <- c()
  for (xv in x_vars) {
    f <- stats::as.formula(paste(y, "~", xv))
    m <- eval(substitute(stats::glm(f, data = DATA, family = stats::binomial),
                        list(DATA = data)))
    td <- broom::tidy(m, conf.int = TRUE, exponentiate = TRUE)
    rows_x <- td[td$term != "(Intercept)", ]
    if (nrow(rows_x) == 0) next
    p_min <- min(rows_x$p.value)
    best <- rows_x[which.min(rows_x$p.value), ]
    univ[[xv]] <- data.frame(
      variable = xv, OR = best$estimate,
      CI.low = best$conf.low, CI.high = best$conf.high,
      p.value = p_min, stringsAsFactors = FALSE)
    if (p_min < alpha) sig <- c(sig, xv)
  }
  univ_df <- do.call(rbind, univ)
  if (length(sig) > 0) {
    f2 <- stats::as.formula(paste(y, "~", paste(sig, collapse = " + ")))
    m2 <- eval(substitute(stats::glm(f2, data = DATA, family = stats::binomial),
                          list(DATA = data)))
    mult <- broom::tidy(m2, conf.int = TRUE, conf.level = conf, exponentiate = TRUE)
    mult$variable <- mult$term
  } else {
    mult <- NULL
    message("Univariate screening selected no variables; multivariable model not built.")
  }
  list(univariate = univ_df, selected = sig, multivariable = mult)
}

#' Stepwise variable selection
#'
#' Corresponds to Chapters 5 and 7: stepwise selection based on AIC
#' (default \code{k = 2}), returning a coefficient table (OR / regression
#' coefficients).
#'
#' @param data A data.frame.
#' @param formula Full model formula.
#' @param direction "both" (default), "backward" or "forward".
#' @param k AIC penalty. Default 2.
#' @return A list: \code{model} (final glm object), \code{coefficients}
#'   (coefficient table), \code{formula} (final formula).
#' @examples
#' bw <- MASS::birthwt
#' mc_stepwise(bw, low ~ age + lwt + race + smoke + ht + ui)
#' @export
mc_stepwise <- function(data, formula, direction = c("both", "backward", "forward"),
                        k = 2) {
  direction <- match.arg(direction)
  fit <- eval(substitute(stats::glm(formula, data = DATA,
                                   family = stats::binomial),
                         list(DATA = data)))
  step_fit <- stats::step(fit, direction = direction, k = k, trace = 0)
  td <- broom::tidy(step_fit, conf.int = TRUE, exponentiate = TRUE)
  td$variable <- td$term
  list(model = step_fit, coefficients = td, formula = formula(step_fit))
}

#' Ordinal logistic regression
#'
#' Corresponds to Chapter 5 ("ordinal logistic regression"). Based on
#' \pkg{MASS}'s \code{polr}.
#'
#' @param data A data.frame.
#' @param formula Model formula; the outcome is an ordered factor.
#' @return A data.frame: variable / Value / Std.Error / t.value / OR / p.value
#' @examples
#' bw <- MASS::birthwt
#' bw$race_ord <- ordered(bw$race)
#' mc_ordinal_logistic(bw, race_ord ~ age + smoke)
#' @export
mc_ordinal_logistic <- function(data, formula) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("`mc_ordinal_logistic` requires the MASS package (shipped with R).")
  }
  fit <- eval(substitute(MASS::polr(formula, data = DATA, Hess = TRUE),
                        list(DATA = data)))
  co <- summary(fit)$coefficients
  cf <- as.data.frame(co)
  cf$variable <- rownames(cf)
  cf$OR <- exp(cf$Value)
  cf$p.value <- 2 * stats::pnorm(-abs(cf$`t value`))
  cf
}
