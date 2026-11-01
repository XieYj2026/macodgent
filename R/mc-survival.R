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
# Module: Survival analysis module (survival)
# This file contains: mc_survfit / mc_cox / mc_competing_risk / mc_multistate /
#              mc_ph_test (Cox proportional-hazards assumption test)
#              mc_competing_risk rebuilds the subdistribution hazard-ratio table and
#              p values from crr coefficients; mc_multistate parses strata from the
#              formula and builds multi-state predictions on its own.
# =============================================================================

#' KM survival curve & log-rank test
#'
#' Corresponds to Chapter 6: survival estimation, between-group log-rank test,
#' and visualisation.
#'
#' @details
#' Purpose: Estimate the Kaplan-Meier survival curve and compare between-group
#' survival differences with the log-rank test; the most basic exploratory tool in
#' survival analysis.
#' Method: The KM method estimates the survival function from the product of
#' conditional survival probabilities; with grouping it automatically performs a
#' log-rank test (whether the group survival curves are overall equal).
#' Arguments: \code{data} is a data frame; \code{time} is the follow-up time column;
#' \code{status} is the event indicator column (usually 1 = event, 0 = censored);
#' \code{group} is an optional grouping column.
#' Results: Returns a list with \code{survfit} (survival object), \code{logrank}
#' (chi-square, df, p value of the log-rank test, present when grouped), and
#' \code{plot} (KM curve drawn by ggplot, with confidence band and group legend).
#' Notes: Censored data (\code{status = 0}) must not be counted as events; the
#' log-rank test is insensitive to late differences (use a weighted test to emphasize
#' distant differences). View the result directly with \code{print(res$plot)}.
#' @param data A `data.frame`.
#' @param time Time variable name.
#' @param status Outcome variable name; 0 = censored, 1 = event.
#' @param group Optional grouping variable name.
#' @return A list with \code{fit} (survfit object), \code{logrank} (test data frame:
#'   test/chisq/p.value; `NULL` when ungrouped), \code{plot} (ggplot KM plot).
#' @examples
#' # Chapter 6: demonstrate KM curve and log-rank test with the ovarian cancer prognostic dataset
#' ovarian <- survival::ovarian
#' ov <- ovarian
#' ov$status <- ifelse(ov$fustat == 2, 0, 1)   # 2 = censored -> 0, 1 = death -> 1
#' res <- mc_survfit(ov, time = "futime", status = "status", group = "rx")
#' cat("---- log-rank test ----\n"); print(res$logrank, row.names = FALSE)
#' print(res$plot)
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
    ggplot2::scale_x_continuous(breaks = .nice_breaks(td$time, 6)) +
    ggplot2::labs(x = "Time", y = "Survival probability", color = "Group") +
    mc_theme()
  if (!is.null(lr_df)) {
    plt <- plt + ggplot2::ggtitle(paste0("log-rank p = ",
                                         format.pval(lr_df$p.value, digits = 3)))
  }
  list(fit = fit, logrank = lr_df, plot = plt)
}

#' Cox PH model (PH test & forest plot)
#'
#' Corresponds to Chapter 6: Cox model, proportional-hazards assumption test
#' (Schoenfeld residuals), and a forest plot.
#'
#' @details
#' Purpose: Fit a Cox proportional-hazards model to assess the concurrent effect of
#' several covariates on the hazard of an event, yielding hazard ratios (HR).
#' Method: A semi-parametric model that makes no assumption about the baseline hazard
#' shape and only estimates covariate hazard ratios; by default it tests the
#' proportional-hazards (PH) assumption.
#' Arguments: \code{data} is a data frame; \code{formula} is of the form
#' \code{Surv(time, status) ~ covariates}, requiring \code{survival::Surv()} to build
#' the survival object; \code{conf} is the HR confidence level.
#' Results: Returns a list with \code{coefficients} (per-variable HR, SE, z, p value
#' and CI), \code{ph_test} (PH-assumption test result), and \code{forest} (HR forest
#' plot, ggplot object).
#' Notes: HR > 1 means the event hazard rises as that covariate increases. The PH
#' assumption is the core premise of Cox -- if \code{ph_test} shows a significant
#' violation (e.g. a very small p for a variable), consider stratified Cox,
#' time-varying effects, or a parametric model. View the forest plot with
#' \code{print(res$forest)}.
#' @param data A `data.frame`.
#' @param formula A Cox model formula, e.g. \code{Surv(time, status) ~ x1 + x2}.
#' @param conf Confidence level; default `0.95`.
#' @return A list with \code{model} (coxph object), \code{coefficients} (HR table:
#'   term/estimate/conf.low/conf.high/p.value), \code{ph_test} (PH test table),
#'   \code{forest} (ggplot forest plot).
#' @examples
#' # Chapter 6: Cox proportional-hazards model -- ovarian-cancer survival prognosis
#' ovarian <- survival::ovarian
#' ov <- ovarian
#' ov$status <- ifelse(ov$fustat == 2, 0, 1)
#' res <- mc_cox(ov, survival::Surv(futime, status) ~ age + rx + resid.ds + ecog.ps)
#' cat("---- HR table / coefficients ----\n"); print(as.data.frame(res$coefficients), row.names = FALSE)
#' cat("\n---- PH assumption test ----\n"); print(res$ph_test)
#' print(res$forest)
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

#' Competing-risk regression (Fine-Gray)
#'
#' Corresponds to Chapter 8 ("competing-risk survival analysis"). Based on
#' \pkg{cmprsk}'s \code{crr}. \code{status}: 0 = censored, >0 = different event
#' types; \code{cause} selects the event of interest.
#'
#' @details
#' Purpose: Perform regression in the presence of competing events, avoiding
#' overestimation of the target-event risk by mistakenly treating competing events as
#' censored. Typical scenario: a relapse analysis where death is the competing event.
#' Method: Based on cmprsk's Fine-Gray subdistribution-hazard model, which treats
#' competing events as a special kind of censoring and directly estimates the
#' cumulative incidence of the target cause.
#' Arguments: \code{data} is a data frame; \code{time} is the time column;
#' \code{status} is the status column (0 = censored, >0 = different event types);
#' \code{group} is the grouping/exposure variable; \code{cause} specifies the target
#' event code (default 1).
#' Results: Returns a list with \code{coefficients} (subdistribution hazard ratio SHR,
#' SE, z, p value per variable).
#' Notes: The \code{status} coding differs from \code{mc_survfit()}/\code{mc_cox()}:
#' here 0 is censored and each other positive integer represents one event type;
#' arrange the data accordingly. If there is only one event type with no competition,
#' an ordinary Cox model suffices. Depends on cmprsk (Suggests).
#' @param data A `data.frame`.
#' @param time Time variable name.
#' @param status Competing-risk status variable name.
#' @param group Covariate (grouping/exposure) variable name.
#' @param cause Event type of interest; default `1`.
#' @return A list with \code{model} (crr object), \code{coefficients} (subdistribution
#'   HR table: variable / HR / coef / p.value).
#' @examples
#' \donttest{
#' # Chapter 8, Section 8.2: competing risks of leukemia relapse (1) and death (2)
#' # Event of interest cause = 1 (relapse); death is a competing event and must not be treated as censored
#' set.seed(126)
#' n <- 300
#' trt     <- sample(c(0,1), n, TRUE)
#' ftime   <- rexp(n, .05)
#' fstatus <- sample(0:2, n, TRUE, c(.4,.35,.25))   # 0 censored, 1 relapse, 2 death (competing)
#' df <- data.frame(ftime = ftime, fstatus = fstatus, trt = trt)
#' res <- mc_competing_risk(df, time = "ftime", status = "fstatus", group = "trt", cause = 1)
#' print(res$coefficients)
#' }
#' @export
mc_competing_risk <- function(data, time, status, group, cause = 1) {
  if (!requireNamespace("cmprsk", quietly = TRUE)) {
    stop("`mc_competing_risk` requires the cmprsk package; please install.packages('cmprsk').")
  }
  ftime <- data[[time]]
  fstatus <- data[[status]]
  cov <- as.matrix(data[[group]])
  colnames(cov) <- group
  fit <- cmprsk::crr(ftime, fstatus, cov1 = cov, failcode = cause)
  cf <- as.numeric(fit$coef)
  names(cf) <- names(fit$coef)
  se <- sqrt(diag(as.matrix(fit$var)))
  z <- cf / se
  pval <- 2 * stats::pnorm(-abs(z))
  coef_df <- data.frame(
    variable = names(cf),
    coef = cf,
    HR = exp(cf),
    se = se,
    p.value = pval, stringsAsFactors = FALSE)
  list(model = fit, coefficients = coef_df)
}

#' Multi-state survival model
#'
#' Corresponds to Chapter 8 ("multi-state survival model"). Requires \code{data}
#' to be in \pkg{mstate}'s expanded long format (with Tstart, Tstop, status,
#' from, to and covariates); \code{trans} is the transition matrix.
#'
#' @details
#' Purpose: Fit a multi-state survival model to study the transition risk between
#' several health states (e.g. healthy -> complication -> death), which is closer to
#' the clinical process than a single endpoint.
#' Method: Based on the mstate package; first expand follow-up data into a from-to
#' transition long table, then model each transition separately with \code{coxph()}.
#' Arguments: \code{data} is the already-expanded multi-state long table (with columns
#' id, Tstart, Tstop, status, from, to, ...); \code{trans} is the transition matrix
#' (built by \code{mstate::transMat()}, defining which transitions between states are
#' allowed); \code{formula} contains the survival object and covariates, and may use
#' \code{strata(trans)} to stratify by transition type.
#' Results: Returns a list with \code{model} (coxph fit per transition) and
#' \code{msfit} (an msfit object for predicting cumulative transition probabilities).
#' Notes: The input must already be in mstate's required long format (not the raw wide
#' table), and the \code{trans} codes must match the from/to columns and
#' \code{strata(trans)}; a wrong construction leads to "all observations in one group"
#' or no estimable transitions. Depends on mstate (Suggests).
#' @param data mstate expanded long-format `data.frame`.
#' @param trans Transition matrix.
#' @param formula A formula, e.g. \code{Surv(time, status) ~ cov} (time/status are columns in data).
#' @param ... Extra arguments passed to \code{survival::coxph}.
#' @return A list with \code{model} (coxph object), \code{msfit} (mstate::msfit object).
#' @examples
#' \donttest{
#' # Chapter 8: three-state model (Healthy H -> Complication I / direct Death D)
#' set.seed(127)
#' library(mstate)
#' n <- 150
#' age <- rnorm(n, 55, 8)
#' t12 <- rexp(n, .15); t23 <- rexp(n, .20); t13 <- rexp(n, .10)
#' rows <- list()
#' for (i in 1:n) {
#'   if (runif(1) < .6) {
#'     rows[[length(rows)+1]] <- data.frame(id=i, Tstart=0,      Tstop=t12[i],        status=1, from=1, to=2, age=age[i])
#'     rows[[length(rows)+1]] <- data.frame(id=i, Tstart=t12[i], Tstop=t12[i]+t23[i], status=1, from=2, to=3, age=age[i])
#'   } else {
#'     rows[[length(rows)+1]] <- data.frame(id=i, Tstart=0, Tstop=t13[i], status=1, from=1, to=3, age=age[i])
#'   }
#' }
#' ms <- do.call(rbind, rows)
#' trans <- transMat(list(c(2,3), c(3), c()), names = c("H","I","D"))
#' ms$trans <- with(ms, ifelse(from == 1 & to == 2, 1, ifelse(from == 1 & to == 3, 2, 3)))
#' fit <- mc_multistate(ms, trans = trans, formula = survival::Surv(Tstart, Tstop, status) ~ age + strata(trans))
#' cat("coxph coefficients:\n")
#' print(summary(fit$model)$coefficients)
#' cat("\nmsfit class:", class(fit$msfit)[1], "\n")
#' }
#' @export
mc_multistate <- function(data, trans, formula, ...) {
  if (!requireNamespace("mstate", quietly = TRUE)) {
    stop("`mc_multistate` requires the mstate package; please install.packages('mstate').")
  }
  fit <- eval(substitute(survival::coxph(formula, data = DATA, ...),
                        list(DATA = data)))
  nd <- data
  fstr <- paste(deparse(formula), collapse = "")
  m <- regmatches(fstr, gregexpr("strata\\(([^)]+)\\)", fstr))[[1]]
  if (length(m)) {
    svar <- sub("strata\\(([^)]+)\\)", "\\1", m[1])
    nd$strata <- data[[svar]]
  }
  ms <- mstate::msfit(fit, newdata = nd, trans = trans)
  list(model = fit, msfit = ms)
}

#' Cox PH-assumption test
#'
#' Corresponds to Chapter 6 (6.5, testing the PH assumption). Runs
#' \code{survival::cox.zph()} on a Cox model and adds a *self-built* ggplot of
#' scaled Schoenfeld residuals versus transformed time.
#'
#' @details
#' Purpose: Test the core assumption of the Cox proportional-hazards model -- whether
#' each covariate's hazard ratio is constant over time.
#' Dependencies: Requires the survival package (already in Imports).
#' Method: Scaled Schoenfeld residuals should show no trend over time;
#' \code{survival::cox.zph()} smooths each curve and runs a correlation test; p > 0.05
#' is usually taken as satisfying the PH assumption.
#' Arguments: \code{fit} is a fitted \code{coxph} object (e.g. returned by
#' \code{mc_cox()}); or pass \code{data} and \code{formula} for the function to fit first.
#' Results: Returns a list -- \code{test} (test table: variable / chi-square / df / p
#' value, including a Global row), \code{plot} (ggplot of scaled Schoenfeld residuals
#' vs time; the \code{which} variable can be specified, defaulting to the first covariate).
#' Notes: With small sample size or few events the PH test has low power; if the PH
#' assumption is violated, consider stratified Cox, time-dependent terms, or an
#' alternative model such as \code{mc_competing_risk()}.
#' @param fit A fitted `coxph` object.
#' @param data A data frame (used only when \code{fit} is not supplied).
#' @param formula A survival formula (used only when \code{fit} is not supplied).
#' @param which Covariate name for the residual plot (default the first).
#' @param transform Time transformation; default `"km"` (Kaplan-Meier).
#' @return A list with the test table and residual plot.
#' @examples
#' # Chapter 6, Section 6.5: test the proportional-hazards assumption of an ovarian-cancer Cox model
#' ovarian <- survival::ovarian
#' ov <- ovarian
#' ov$status <- ifelse(ov$fustat == 2, 0, 1)
#' fit <- survival::coxph(survival::Surv(futime, status) ~ age + rx + resid.ds, data = ov)
#' r <- mc_ph_test(fit)
#' print(r$test)        # PH test table for each covariate
#' print(r$plot)        # scaled Schoenfeld residuals vs time
#' @export
mc_ph_test <- function(fit = NULL, data = NULL, formula = NULL,
                       which = NULL, transform = "km") {
  if (is.null(fit)) {
    if (is.null(data) || is.null(formula))
      stop("Provide either a fitted `coxph` object or both `data` and `formula`.")
    if (!requireNamespace("survival", quietly = TRUE))
      stop("`mc_ph_test` requires the survival package.")
    fit <- survival::coxph(formula, data = data)
  }
  if (!inherits(fit, "coxph"))
    stop("`fit` must be a coxph object.")
  zp <- survival::cox.zph(fit, transform = transform)
  tb <- as.data.frame(zp$table)
  names(tb) <- c("chisq", "df", "p_value")
  tb$variable <- rownames(tb)
  vnames <- colnames(zp$y)   # covariate names (columns of cox.zph$y = variables, rows = event times)
  w <- if (is.null(which)) vnames[1] else which
  if (!w %in% vnames) stop("`which` must be one of: ", paste(vnames, collapse = ", "))
  yy <- as.numeric(zp$y[, w, drop = TRUE])
  tt <- as.numeric(zp$x)
  dfp <- data.frame(time = tt, sresid = yy)
  plot <- ggplot2::ggplot(dfp, ggplot2::aes(time, sresid)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.6) +
    ggplot2::geom_smooth(se = TRUE, color = "steelblue", size = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::labs(x = paste("Transformed time (", transform, ")", sep = ""),
                  y = paste("Scaled Schoenfeld residual (", w, ")", sep = "")) +
    mc_theme()
  list(test = tb[, c("variable", "chisq", "df", "p_value")], plot = plot)
}
