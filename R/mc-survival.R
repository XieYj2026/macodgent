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
# Module: Survival analysis module (survival)
# This file contains: mc_kaplan_meier / mc_cox / mc_competing_risk / mc_multistate /
# mc_test_proportional_hazards (Cox proportional-hazards assumption test)
# mc_competing_risk rebuilds the subdistribution hazard-ratio table and
# p values from crr coefficients; mc_multistate parses strata from the
# formula and builds multi-state predictions on its own.
# =============================================================================

#' KM survival curve & log-rank test
#'
#' Survival estimation, between-group log-rank test,
#' And visualisation.
#'
#' @details
#' Purpose: Estimate the Kaplan-Meier survival curve and compare between-group
#' Survival differences with the log-rank test; the most basic exploratory tool in
#' Survival analysis.
#' Method: The KM method estimates the survival function from the product of
#' Conditional survival probabilities; with grouping it automatically performs a
#' Log-rank test (whether the group survival curves are overall equal).
#' Arguments: \code{data} is a data frame; \code{time} is the follow-up time column;
#' \code{status} is the event indicator column (usually 1 = event, 0 = censored);
#' \code{group} is an optional grouping column.
#' Results: Returns a list with \code{survfit} (survival object), \code{logrank}
#' (chi-square, df, p value of the log-rank test, present when grouped), and
#' \code{plot} (KM curve drawn by ggplot, with confidence band and group legend).
#' Notes: Censored data (\code{status = 0}) must not be counted as events; the
#' Log-rank test is insensitive to late differences (use a weighted test to emphasize
#' Distant differences). View the result directly with \code{print(res$plot)}.
#' @param data A `data.frame`.
#' @param time Time variable name.
#' @param status Outcome variable name; 0 = censored, 1 = event.
#' @param group Optional grouping variable name.
#' @param times Landmark time points for survival estimates and numbers at
#'   risk; `NULL` = the quartiles of follow-up time.
#' @return A list with \code{fit} (survfit object), \code{logrank} (test data
#'   frame: test/chisq/df/p.value; `NULL` when ungrouped), \code{plot} (ggplot
#'   KM plot), \code{median_survival} (median with 95% CI per stratum),
#'   \code{survival_at} (landmark survival with CI) and \code{risk_table}
#'   (numbers at risk at those landmarks).
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' km <- mc_kaplan_meier(d, time = "followup_days", status = "cv_event",
#'                       group = "treatment")
#' km$logrank
#' km$median_survival
#' @export
mc_kaplan_meier <- function(data, time, status, group = NULL, times = NULL) {
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
 df = length(lr$n) - 1L,
 p.value = as.numeric(lr$pvalue),
 stringsAsFactors = FALSE)
 td <- broom::tidy(fit)
 }
 plt <- ggplot2::ggplot(td, ggplot2::aes(x = time, y = estimate, color = strata)) +
 ggplot2::geom_step() +
 ggplot2::scale_x_continuous(breaks =.nice_breaks(td$time, 6)) +
 ggplot2::labs(x = "Time", y = "Survival probability", color = "Group") +
 mc_theme()
 if (!is.null(lr_df)) {
 plt <- plt + ggplot2::ggtitle(paste0("log-rank p = ",
 format.pval(lr_df$p.value, digits = 3)))
 }
 ## Median survival with its CI, landmark survival probabilities and the
 ## numbers at risk are the numbers a KM figure must carry in a manuscript.
 st <- summary(fit)$table
 if (is.null(dim(st))) {
 lcl <- grep("LCL", names(st), value = TRUE)[1]
 ucl <- grep("UCL", names(st), value = TRUE)[1]
 med <- data.frame(strata = "All", median = unname(st["median"]),
 conf.low = unname(st[lcl]), conf.high = unname(st[ucl]),
 stringsAsFactors = FALSE)
 } else {
 lcl <- grep("LCL", colnames(st), value = TRUE)[1]
 ucl <- grep("UCL", colnames(st), value = TRUE)[1]
 med <- data.frame(strata = rownames(st), median = st[, "median"],
 conf.low = st[, lcl], conf.high = st[, ucl],
 stringsAsFactors = FALSE)
 }
 if (is.null(times)) {
 times <- as.numeric(stats::quantile(data[[time]],
 c(0.25, 0.5, 0.75), na.rm = TRUE))
 }
 times <- sort(unique(as.numeric(times)))
 sm <- summary(fit, times = times)
 strata_col <- if (!is.null(sm$strata)) as.character(sm$strata) else
 rep("All", length(sm$time))
 landmark <- data.frame(
 strata = strata_col, time = as.numeric(sm$time),
 n.risk = as.numeric(sm$n.risk),
 survival = as.numeric(sm$surv),
 conf.low = as.numeric(sm$lower), conf.high = as.numeric(sm$upper),
 stringsAsFactors = FALSE)
 list(fit = fit, logrank = lr_df, plot = plt,
 median_survival = med, survival_at = landmark,
 risk_table = landmark[, c("strata", "time", "n.risk")])
}

#' Cox PH model (PH test & forest plot)
#'
#' Cox model, proportional-hazards assumption test
#' (Schoenfeld residuals), and a forest plot.
#'
#' @details
#' Purpose: Fit a Cox proportional-hazards model to assess the concurrent effect of
#' Several covariates on the hazard of an event, yielding hazard ratios (HR).
#' Method: A semi-parametric model that makes no assumption about the baseline hazard
#' Shape and only estimates covariate hazard ratios; by default it tests the
#' Proportional-hazards (PH) assumption.
#' Arguments: \code{data} is a data frame; \code{formula} is of the form
#' \code{Surv(time, status) ~ covariates}, requiring \code{survival::Surv()} to build
#' The survival object; \code{conf} is the HR confidence level.
#' Results: Returns a list with \code{coefficients} (per-variable HR, SE, z, p value
#' And CI), \code{ph_test} (PH-assumption test result), and \code{forest} (HR forest
#' Plot, ggplot object).
#' Notes: HR > 1 means the event hazard rises as that covariate increases. The PH
#' Assumption is the core premise of Cox -- if \code{ph_test} shows a significant
#' Violation (e.g. a very small p for a variable), consider stratified Cox,
#' Time-varying effects, or a parametric model. View the forest plot with
#' \code{print(res$forest)}.
#' @param data A `data.frame`.
#' @param formula A Cox model formula, e.g. \code{Surv(time, status) ~ x1 + x2}.
#' @param conf Confidence level; default `0.95`.
#' @return A list with \code{model} (coxph object), \code{coefficients} (HR
#'   table: term/estimate/conf.low/conf.high/p.value), \code{ph_test} (PH test
#'   table), \code{forest} (ggplot forest plot), \code{performance} (n, events,
#'   concordance, likelihood-ratio test, AIC, R-squared) and \code{checks}
#'   (explicit PH-assumption and global-fit verdicts).
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' cx <- mc_cox(d,
#'   survival::Surv(followup_days, cv_event) ~ age + bmi + treatment + diabetes)
#' cx$coefficients
#' cx$performance
#' @export
mc_cox <- function(data, formula, conf = 0.95) {
 if (!requireNamespace("survival", quietly = TRUE)) stop("The survival package is required.")
 fit <- eval(substitute(survival::coxph(formula, data = DATA),
 list(DATA = data)))
 td <- broom::tidy(fit, conf.int = TRUE, conf.level = conf, exponentiate = TRUE)
 td$variable <- td$term
 zph <- survival::cox.zph(fit)
 zph_df <- as.data.frame(zph$table)
 smc <- summary(fit)
 perf <- data.frame(
 metric = c("n", "events", "concordance", "concordance_se",
 "LRT_chi2", "LRT_df", "LRT_p", "AIC", "R2"),
 value = c(as.numeric(smc$n), as.numeric(smc$nevent),
 unname(smc$concordance[1]), unname(smc$concordance[2]),
 unname(smc$logtest[1]), unname(smc$logtest[2]),
 unname(smc$logtest[3]), stats::AIC(fit), unname(smc$rsq[1])),
 stringsAsFactors = FALSE)
 pcol <- grep("^p$", colnames(zph_df), value = TRUE)[1]
 viol <- if (is.na(pcol)) character(0) else
 rownames(zph_df)[!is.na(zph_df[[pcol]]) & zph_df[[pcol]] < 0.05 &
 rownames(zph_df) != "GLOBAL"]
 checks <- data.frame(
 check = c("Proportional hazards", "Global fit"),
 value = c(length(viol), unname(smc$logtest[3])),
 status = c(if (length(viol) == 0L) "PASS" else "REVIEW",
 if (!is.na(smc$logtest[3]) && smc$logtest[3] < 0.05) "PASS" else "REVIEW"),
 detail = c(if (length(viol) == 0L)
 "no covariate violates PH at alpha = 0.05" else
 paste("PH violation:", paste(viol, collapse = ", ")),
 "likelihood-ratio test against the null model"),
 stringsAsFactors = FALSE)
 fp <- ggplot2::ggplot(td, ggplot2::aes(x = variable, y = estimate)) +
 ggplot2::geom_pointrange(ggplot2::aes(ymin = conf.low, ymax = conf.high)) +
 ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#C0392B") +
 ggplot2::coord_flip() +
 ggplot2::labs(x = "", y = "Hazard Ratio (95% CI)") +
 mc_theme()
 list(model = fit, coefficients = td, ph_test = zph_df, forest = fp,
 performance = perf, checks = checks)
}

#' Competing-risk regression (Fine-Gray)
#'
#' Based on
#' \pkg{cmprsk}'s \code{crr}. \code{status}: 0 = censored, >0 = different event
#' Types; \code{cause} selects the event of interest.
#'
#' @details
#' Purpose: Perform regression in the presence of competing events, avoiding
#' Overestimation of the target-event risk by mistakenly treating competing events as
#' Censored. Typical scenario: a relapse analysis where death is the competing event.
#' Method: Based on cmprsk's Fine-Gray subdistribution-hazard model, which treats
#' Competing events as a special kind of censoring and directly estimates the
#' Cumulative incidence of the target cause.
#' Arguments: \code{data} is a data frame; \code{time} is the time column;
#' \code{status} is the status column (0 = censored, >0 = different event types);
#' \code{group} is the grouping/exposure variable; \code{cause} specifies the target
#' Event code (default 1).
#' Results: Returns a list with \code{coefficients} (subdistribution hazard ratio SHR,
#' SE, z, p value per variable).
#' Notes: The \code{status} coding differs from \code{mc_kaplan_meier()}/\code{mc_cox()}:
#' Here 0 is censored and each other positive integer represents one event type;
#' Arrange the data accordingly. If there is only one event type with no competition,
#' An ordinary Cox model suffices. Depends on cmprsk (Suggests).
#' @param data A `data.frame`.
#' @param time Time variable name.
#' @param status Competing-risk status variable name.
#' @param group Covariate (grouping/exposure) variable name.
#' @param cause Event type of interest; default `1`.
#' @return A list with \code{model} (crr object), \code{coefficients} (subdistribution
#' HR table: variable / HR / coef / p.value).
#' @param conf Confidence level for the subdistribution HR interval; default `0.95`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' \donttest{
#' cr <- mc_competing_risk(d, time = "followup_days", status = "comp_event",
#'                         group = "treatment", cause = 1)
#' cr$coefficients
#' cr$checks
#' }
#' @export
mc_competing_risk <- function(data, time, status, group, cause = 1,
                              conf = 0.95) {
 if (!requireNamespace("cmprsk", quietly = TRUE)) {
 stop("`mc_competing_risk` requires the cmprsk package; please install.packages('cmprsk').")
 }
  ftime <- data[[time]]
  fstatus <- data[[status]]
  # crr() passes cov1 to Fortran, so it must be numeric. as.matrix() on a
  # factor column yields a character matrix and fails with "NA/NaN/Inf in
  # foreign function call"; model.matrix() dummy-codes factors and keeps
  # numeric grouping variables unchanged.
  cov <- stats::model.matrix(~ ., data = as.data.frame(data[group]))[, -1L,
                                                                     drop = FALSE]
  if (!ncol(cov)) stop("`group` must provide at least one usable covariate.")
  fit <- cmprsk::crr(ftime, fstatus, cov1 = cov, failcode = cause)
 cf <- as.numeric(fit$coef)
 names(cf) <- names(fit$coef)
 se <- sqrt(diag(as.matrix(fit$var)))
 z <- cf / se
 pval <- 2 * stats::pnorm(-abs(z))
 zc <- stats::qnorm(1 - (1 - conf) / 2)
 coef_df <- data.frame(
 variable = names(cf), coef = cf, se = se, z = z,
 HR = exp(cf), HR_lower = exp(cf - zc * se),
 HR_upper = exp(cf + zc * se),
 p.value = pval, stringsAsFactors = FALSE)
 ## Cumulative incidence function for the event of interest.
 grp <- if (length(group) == 1L &&
 length(unique(data[[group[1L]]])) <= 5L)
 as.factor(data[[group[1L]]]) else NULL
 cif <- tryCatch({
 ci <- if (is.null(grp)) cmprsk::cuminc(ftime, fstatus, cencode = 0) else
 cmprsk::cuminc(ftime, fstatus, grp, cencode = 0)
 keys <- setdiff(names(ci)[grepl(paste0(" ", cause, "$"), names(ci))],
 "Tests")
 do.call(rbind, lapply(keys, function(k) {
 zz <- ci[[k]]
 data.frame(group = sub(paste0(" ", cause, "$"), "", k),
 time = as.numeric(zz$time), cif = as.numeric(zz$est),
 var = as.numeric(zz$var), stringsAsFactors = FALSE)
 }))
 }, error = function(e) NULL)
 tab <- table(fstatus)
 checks <- data.frame(
 metric = c("n", "n_events_cause", "n_competing", "n_censored",
 "events_per_variable", "n_coefficients"),
 value = c(length(ftime),
 sum(fstatus == cause), sum(fstatus > 0 & fstatus != cause),
 sum(fstatus == 0),
 sum(fstatus == cause) / max(1L, length(cf)), length(cf)),
 stringsAsFactors = FALSE)
 list(model = fit, coefficients = coef_df, cif = cif,
 cause_counts = tab, checks = checks)
}

#' Multi-state survival model
#'
#' Requires \code{data}
#' To be in \pkg{mstate}'s expanded long format (with Tstart, Tstop, status,
#' From, to and covariates); \code{trans} is the transition matrix.
#'
#' @details
#' Purpose: Fit a multi-state survival model to study the transition risk between
#' Several health states (e.g. healthy -> complication -> death), which is closer to
#' The clinical process than a single endpoint.
#' Method: Based on the mstate package; first expand follow-up data into a from-to
#' Transition long table, then model each transition separately with \code{coxph()}.
#' Arguments: \code{data} is the already-expanded multi-state long table (with columns
#' Id, Tstart, Tstop, status, from, to,...); \code{trans} is the transition matrix
#' (built by \code{mstate::transMat()}, defining which transitions between states are
#' Allowed); \code{formula} contains the survival object and covariates, and may use
#' \code{strata(trans)} to stratify by transition type.
#' Results: Returns a list with \code{model} (coxph fit per transition) and
#' \code{msfit} (an msfit object for predicting cumulative transition probabilities).
#' Notes: The input must already be in mstate's required long format (not the raw wide
#' Table), and the \code{trans} codes must match the from/to columns and
#' \code{strata(trans)}; a wrong construction leads to "all observations in one group"
#' Or no estimable transitions. Depends on mstate (Suggests).
#' @param data mstate expanded long-format `data.frame`.
#' @param trans Transition matrix.
#' @param formula A formula, e.g. \code{Surv(time, status) ~ cov} (time/status are columns in data).
#' @param ... Extra arguments passed to \code{survival::coxph}.
#' @return A list with \code{model} (coxph object), \code{msfit} (mstate::msfit object).
#' @param times Prediction time point(s) passed to `mstate::probtrans()`; default `0`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' \donttest{
#' ms_d <- d[1:200, ]                  # simulated 200-subject subset
#' set.seed(127)
#' n <- nrow(ms_d)
#' age <- ms_d$age
#' t12 <- rexp(n, .15); t23 <- rexp(n, .20); t13 <- rexp(n, .10)
#' rows <- list()
#' for (i in seq_len(n)) {
#'   if (runif(1) < .6) {
#'     rows[[length(rows) + 1]] <- data.frame(id = i, tstart = 0,
#'       tstop = t12[i], status = 1, from = 1, to = 2, age = age[i])
#'     rows[[length(rows) + 1]] <- data.frame(id = i, tstart = t12[i],
#'       tstop = t12[i] + t23[i], status = 1, from = 2, to = 3, age = age[i])
#'   } else {
#'     rows[[length(rows) + 1]] <- data.frame(id = i, tstart = 0,
#'       tstop = t13[i], status = 1, from = 1, to = 3, age = age[i])
#'   }
#' }
#' ms <- do.call(rbind, rows)
#' trans <- mstate::transMat(list(c(2, 3), c(3), c()),
#'                           names = c("h", "i", "d"))
#' ms$trans <- with(ms, ifelse(from == 1 & to == 2, 1,
#'                      ifelse(from == 1 & to == 3, 2, 3)))
#' fit <- mc_multistate(ms, trans = trans,
#'   formula = survival::Surv(tstart, tstop, status) ~
#'     age + survival::strata(trans))
#' head(fit$probability_table)
#' }
#' @export
mc_multistate <- function(data, trans, formula, ..., times = 0) {
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
 ## State-occupation probabilities; probtrans needs a strata column in
 ## newdata, so it is guarded and reported as unavailable when absent.
 tp <- tryCatch(mstate::probtrans(ms, predt = times),
 error = function(e) NULL)
 prob_table <- NULL
 if (!is.null(tp)) {
 first <- if (is.data.frame(tp)) tp else tp[[1L]]
 prob_table <- tryCatch(as.data.frame(first), error = function(e) NULL)
 }
 n_trans <- tryCatch(sum(!is.na(ms$Haz$trans)),
 error = function(e) NA_real_)
 checks <- data.frame(
 metric = c("n_rows", "n_transitions", "n_coefficients",
 "probabilities_available"),
 value = c(nrow(data), n_trans, length(stats::coef(fit)),
 as.numeric(!is.null(prob_table))),
 stringsAsFactors = FALSE)
 list(model = fit, msfit = ms,
 transition_probabilities = tp, probability_table = prob_table,
 checks = checks)
}

#' Cox PH-assumption test
#'
#' Runs
#' \code{survival::cox.zph()} on a Cox model and adds a *self-built* ggplot of
#' Scaled Schoenfeld residuals versus transformed time.
#'
#' @details
#' Purpose: Test the core assumption of the Cox proportional-hazards model -- whether
#' Each covariate's hazard ratio is constant over time.
#' Dependencies: Requires the survival package (already in Imports).
#' Method: Scaled Schoenfeld residuals should show no trend over time;
#' \code{survival::cox.zph()} smooths each curve and runs a correlation test; p > 0.05
#' Is usually taken as satisfying the PH assumption.
#' Arguments: \code{fit} is a fitted \code{coxph} object (e.g. returned by
#' \code{mc_cox()}); or pass \code{data} and \code{formula} for the function to fit first.
#' Results: Returns a list -- \code{test} (test table: variable / chi-square / df / p
#' Value, including a Global row), \code{plot} (ggplot of scaled Schoenfeld residuals
#' Vs time; the \code{which} variable can be specified, defaulting to the first covariate).
#' Notes: With small sample size or few events the PH test has low power; if the PH
#' Assumption is violated, consider stratified Cox, time-dependent terms, or an
#' Alternative model such as \code{mc_competing_risk()}.
#' @param fit A fitted `coxph` object.
#' @param data A data frame (used only when \code{fit} is not supplied).
#' @param formula A survival formula (used only when \code{fit} is not supplied).
#' @param which Covariate name for the residual plot (default the first).
#' @param transform Time transformation; default `"km"` (Kaplan-Meier).
#' @return A list with the test table and residual plot.
#' @param alpha Significance level for the PASS / REVIEW verdict; default `0.05`.
#' @param p_adjust Multiplicity correction for per-covariate tests; default `"BH"`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' cx <- mc_cox(d,
#'   survival::Surv(followup_days, cv_event) ~ age + bmi + treatment)
#' mc_test_proportional_hazards(fit = cx$model)
#' @export
mc_test_proportional_hazards <- function(fit = NULL, data = NULL, formula = NULL,
 which = NULL, transform = "km", alpha = 0.05,
 p_adjust = "BH") {
 if (is.null(fit)) {
 if (is.null(data) || is.null(formula))
 stop("Provide either a fitted `coxph` object or both `data` and `formula`.")
 fit <- survival::coxph(formula, data = data)
 }
 if (!inherits(fit, "coxph"))
 stop("`fit` must be a coxph object.")
 zp <- survival::cox.zph(fit, transform = transform)
 tb <- as.data.frame(zp$table)
 names(tb) <- c("chisq", "df", "p_value")
 tb$variable <- rownames(tb)
 vnames <- colnames(zp$y)
 w <- if (is.null(which)) vnames[1] else which
 if (!w %in% vnames) stop("`which` must be one of: ", paste(vnames, collapse = ", "))
 ## Multiplicity correction across the covariate-specific tests only.
 per <- tb$variable != "GLOBAL"
 tb$p.adjusted <- NA_real_
 tb$p.adjusted[per] <- stats::p.adjust(tb$p_value[per], method = p_adjust)
 global_p <- if ("GLOBAL" %in% tb$variable)
 tb$p_value[tb$variable == "GLOBAL"][1] else NA_real_
 verdict <- .mc_check_row("Proportional hazards (global)", global_p,
 "global PH p > alpha",
 if (is.na(global_p)) NA else global_p > alpha,
 "Schoenfeld residuals; per-covariate p values in `test`")
 checks <- rbind(verdict,
 .mc_check_row("Per-covariate PH", suppressWarnings(min(tb$p.adjusted[per], na.rm = TRUE)),
 "min adjusted p > alpha",
 if (all(is.na(tb$p.adjusted[per]))) NA else
 min(tb$p.adjusted[per], na.rm = TRUE) > alpha))
 yy <- as.numeric(zp$y[, w, drop = TRUE])
 tt <- as.numeric(zp$x)
 dfp <- data.frame(time = tt, sresid = yy)
 plot <- ggplot2::ggplot(dfp, ggplot2::aes(time, sresid)) +
 ggplot2::geom_point(size = 1.2, alpha = 0.6) +
 ggplot2::geom_smooth(se = TRUE, color = "steelblue", linewidth = 0.8) +
 ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
 ggplot2::labs(x = paste("Transformed time (", transform, ")", sep = ""),
 y = paste("Scaled Schoenfeld residual (", w, ")", sep = "")) +
 mc_theme()
 list(test = tb[, c("variable", "chisq", "df", "p_value", "p.adjusted")],
 global_p = global_p, verdict = verdict, checks = checks,
 p_adjust = p_adjust, plot = plot)
}
