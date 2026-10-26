#' Kaplan-Meier survival curve and log-rank test
#'
#' Corresponds to Chapter 6: survival estimation, between-group log-rank test,
#' and visualisation.
#'
#' @param data A data.frame.
#' @param time Time variable name.
#' @param status Outcome variable name (0 = censored, 1 = event).
#' @param group Optional grouping variable name.
#' @return A list: \code{fit} (survfit object), \code{logrank} (test data.frame),
#'   \code{plot} (ggplot KM plot).
#' @examples
#' mc_survfit(survival::lung, time = "time", status = "status", group = "sex")
#' @export
mc_survfit <- function(data, time, status, group = NULL) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("The survival package is required.")
  sv <- survival::Surv(data[[time]], data[[status]])
  if (is.null(group)) {
    fit <- survival::survfit(sv ~ 1)
    lr_df <- NULL
    td <- broom::tidy(fit)
    td$strata <- "All"
  } else {
    g <- data[[group]]
    tmp <- data.frame(sv = sv, g = g)
    fit <- survival::survfit(sv ~ g, data = tmp)
    lr <- survival::survdiff(sv ~ g, data = tmp)
    lr_df <- data.frame(test = "log-rank",
                        chisq = as.numeric(lr$chisq),
                        p.value = as.numeric(lr$pvalue),
                        stringsAsFactors = FALSE)
    td <- broom::tidy(fit)
  }
  plt <- ggplot2::ggplot(td, ggplot2::aes(x = time, y = estimate, color = strata)) +
    ggplot2::geom_step() +
    ggplot2::labs(x = "Time", y = "Survival probability", color = "Group") +
    mc_theme()
  if (!is.null(lr_df)) {
    plt <- plt + ggplot2::ggtitle(paste0("log-rank p = ",
                                         format.pval(lr_df$p.value, digits = 3)))
  }
  list(fit = fit, logrank = lr_df, plot = plt)
}

#' Cox proportional hazards model (with PH test and forest plot)
#'
#' Corresponds to Chapter 6: Cox model, proportional-hazards assumption test
#' (Schoenfeld residuals), and a forest plot.
#'
#' @param data A data.frame.
#' @param formula A Cox model formula, e.g. \code{Surv(time, status) ~ x1 + x2}.
#' @param conf Confidence level. Default 0.95.
#' @return A list: \code{model} (coxph object), \code{coefficients} (HR table),
#'   \code{ph_test} (PH test table), \code{forest} (ggplot forest plot).
#' @examples
#' fit <- mc_cox(survival::lung, survival::Surv(time, status) ~ age + sex + ph.ecog)
#' fit$forest
#' @export
mc_cox <- function(data, formula, conf = 0.95) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("The survival package is required.")
  fit <- eval(substitute(survival::coxph(formula, data = DATA),
                        list(DATA = data)))
  td <- broom::tidy(fit, conf.int = TRUE, conf.level = conf, exponentiate = TRUE)
  td$variable <- td$term
  zph <- survival::cox.zph(fit)
  zph_df <- as.data.frame(zph$table)
  fp <- ggplot2::ggplot(td, ggplot2::aes(x = variable, y = estimate)) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = conf.low, ymax = conf.high)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#C0392B") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "", y = "Hazard Ratio (95% CI)") +
    mc_theme()
  list(model = fit, coefficients = td, ph_test = zph_df, forest = fp)
}

#' Competing-risk regression (Fine-Gray model)
#'
#' Corresponds to Chapter 8 ("competing-risk survival analysis"). Based on
#' \pkg{cmprsk}'s \code{crr}. \code{status}: 0 = censored, >0 = different event
#' types; \code{cause} selects the event of interest.
#'
#' @param data A data.frame.
#' @param time Time variable name.
#' @param status Competing-risk status variable name.
#' @param group Covariate (grouping/exposure) variable name.
#' @param cause Event type of interest. Default 1.
#' @return A list: \code{model} (crr object), \code{coefficients}
#'   (subdistribution HR table).
#' @examples
#' \dontrun{
#' mc_competing_risk(df, time = "ftime", status = "fstatus", group = "trt", cause = 1)
#' }
#' @export
mc_competing_risk <- function(data, time, status, group, cause = 1) {
  if (!requireNamespace("cmprsk", quietly = TRUE)) {
    stop("`mc_competing_risk` requires the cmprsk package; please install.packages('cmprsk').")
  }
  ftime <- data[[time]]
  fstatus <- data[[status]]
  cov <- stats::model.matrix(~ data[[group]])[, -1, drop = FALSE]
  fit <- cmprsk::crr(ftime, fstatus, cov1 = cov, failcode = cause)
  cf <- as.data.frame(fit$coef)
  coef_df <- data.frame(
    variable = rownames(cf), HR = cf$`exp(coef)`,
    coef = cf$coef, p.value = cf$`p-value`, stringsAsFactors = FALSE)
  list(model = fit, coefficients = coef_df)
}

#' Multi-state survival model
#'
#' Corresponds to Chapter 8 ("multi-state survival model"). Requires \code{data}
#' to be in \pkg{mstate}'s expanded long format (with Tstart, Tstop, status,
#' from, to and covariates); \code{trans} is the transition matrix. Returns the
#' Cox fit and \code{msfit} results.
#'
#' @param data mstate expanded long-format data.frame.
#' @param trans Transition matrix.
#' @param formula A formula like \code{Surv(time, status) ~ cov} (time/status are
#'   columns in data).
#' @param ... Extra arguments passed to \code{survival::coxph}.
#' @return A list: \code{model} (coxph object), \code{msfit} (mstate::msfit object).
#' @examples
#' \dontrun{
#' mc_multistate(mstate_data, trans = trans_mat, formula = Surv(time, status) ~ age)
#' }
#' @export
mc_multistate <- function(data, trans, formula, ...) {
  if (!requireNamespace("mstate", quietly = TRUE)) {
    stop("`mc_multistate` requires the mstate package; please install.packages('mstate').")
  }
  fit <- eval(substitute(survival::coxph(formula, data = DATA, ...),
                        list(DATA = data)))
  ms <- mstate::msfit(fit, trans = trans)
  list(model = fit, msfit = ms)
}
