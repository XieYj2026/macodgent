# =============================================================================
# Software name : macodgent medical data-science research toolkit
# Copyright owner : Dr.Xie
# Version : 1.0.0
# Language : R
# Completion date : 2026-10-30
# Copyright notice : This program is independently developed by the copyright
# owner and protected under the Software Copyright Protection
# Regulation; without written permission, it may not be
# copied, distributed, rented, or used for commercial purposes.
# -----------------------------------------------------------------------------
# Module : Regression modelling module
# This file contains: mc_logistic / mc_logistic_screen / mc_stepwise / mc_ordinal_logistic
# mc_clogistic (conditional logistic / matched case-control) /
# mc_model_compare (original model comparison)
# mc_logistic_screen (univariate screening then multivariable
# modelling) is an original combined workflow.
# =============================================================================

#' Binary logistic regression
#'
#' Based on \code{glm(family = binomial)}, returns a
#' Tidy coefficient table (with OR and its 95% CI, already exponentiated).
#'
#' @details
#' Purpose: Fit a binary logistic regression to obtain the log-odds ratio and OR
#' Of each predictor for the outcome.
#' Method: A binomial logit model estimated by maximum likelihood; factor/
#' Character variables are automatically handled as categorical.
#' Arguments: `data` is a data.frame; `formula` is of the form
#' Outcome ~ predictors (the outcome must be 0/1 or a two-level factor);
#' `conf` is the confidence level for the OR, default 0.95.
#' Result: Returns a tidy coefficient table with estimate (on the logit scale),
#' Standard error, z value, p value, and the exponentiated OR with its
#' Confidence interval, for direct interpretation.
#' Note: OR > 1 means the odds of the outcome increase by one unit of the
#' Predictor (or relative to the reference level); a CI crossing 1 indicates the
#' Predictor is not significant at alpha = 0.05. Check for severe
#' Multicollinearity and empty cells before modelling. For non-binary outcomes,
#' Use a multinomial or ordinal model instead.
#' @param data A \code{data.frame}.
#' @param formula A model formula, e.g. \code{y ~ x1 + x2}.
#' @param conf Confidence level; default \code{0.95}.
#' @return A tidy \code{data.frame} with columns \code{term} / \code{estimate}
#' (OR) / \code{std.error} / \code{statistic} / \code{p.value} / \code{conf.low}
#' / \code{conf.high}; \code{estimate} is the exponentiated OR.
#' @param diagnostics If `TRUE`, return the full diagnostic bundle from [mc_logistic_report()] instead of only the coefficient table.
#' @param auc_boot Bootstrap resamples for the AUC CI; `0` uses the Hanley-McNeil approximation.
#' @param seed Random seed used when `diagnostics = TRUE`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' mc_logistic(d, bp_control ~ age + bmi + treatment)
#' mc_logistic(d, bp_control ~ age + bmi + treatment,
#'             diagnostics = TRUE, auc_boot = 0L)$checks
#' @export
mc_logistic <- function(data, formula, conf = 0.95, diagnostics = FALSE,
                        auc_boot = 0L, seed = 1L) {
  fit <- eval(substitute(stats::glm(formula, data = DATA,
                                    family = stats::binomial),
                         list(DATA = data)))
  td <- broom::tidy(fit, conf.int = TRUE, conf.level = conf,
                    exponentiate = TRUE)
  if (!isTRUE(diagnostics)) return(td)
  rep <- mc_logistic_report(data, formula, conf = conf,
                            auc_boot = auc_boot, seed = seed)
  list(coefficients = td, model = fit, report = rep,
       checks = rep$checks, performance = rep$performance,
       diagnostics = rep$diagnostics)
}

#' Univariate screen -> multivariable logistic regression
#'
#' Run a univariate logistic
#' Regression for each candidate variable, screen significant ones at \code{alpha},
#' Then build a multivariable model from the selected variables.
#'
#' @details
#' Purpose: One-step implementation of the classic epidemiological workflow
#' "univariate screening then multivariable modelling", removing the need for
#' Manual stepwise work.
#' Method: First run a univariate logistic regression for each candidate
#' `x_vars`, screen significant variables at the `alpha` threshold (default 0.05),
#' Then fit a multivariable model with the selected variables.
#' Arguments: `data` is a data.frame; `y` is the binary outcome column name
#' (character); `x_vars` is a vector of candidate predictor column names;
#' `alpha` is the significance threshold for univariate screening; `conf` is the
#' Confidence level.
#' Result: Returns a list with `univariate` (univariate results per variable),
#' `selected` (names of selected variables), and `multivariable` (the tidy
#' Multivariable coefficient table).
#' Note: This is convenient but not a gold standard  -  stepwise/screening can
#' Introduce selection bias and overfitting, and univariate significance does not
#' Imply multivariable significance (variables may mask or amplify each other).
#' Report both univariate and multivariable results and state the screening rule.
#' It is still recommended to validate the final model with mc_logistic().
#' @param data A \code{data.frame}.
#' @param y Binary outcome variable name (character).
#' @param x_vars Candidate predictor names (character vector).
#' @param alpha Univariate screening significance threshold; default \code{0.05}.
#' @param conf Multivariable confidence level; default \code{0.95}.
#' @return A list with \code{univariate} (univariate table), \code{selected}
#' (selected variables), and \code{multivariable} (multivariable coefficient
#' Table; \code{NULL} if no variables were selected).
#' @param p_adjust Multiplicity correction for the univariate p values; `"none"` (default), `"holm"`, `"BH"`, etc.
#' @param bootstrap Number of resamples for variable-selection stability; `0` (default) skips it.
#' @param seed Random seed for the stability resamples.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' scr <- mc_logistic_screen(d, y = "bp_control",
#'   x_vars = c("age", "bmi", "ldl", "glucose", "smoke", "diabetes"),
#'   p_adjust = "BH")
#' scr$univariate
#' scr$selected
#' @export
mc_logistic_screen <- function(data, y, x_vars, alpha = 0.05, conf = 0.95,
                               p_adjust = "none", bootstrap = 0L,
                               seed = NULL) {
 univ <- list()
 for (xv in x_vars) {
 f <- stats::as.formula(paste(y, "~", xv))
 m <- tryCatch(stats::glm(f, data = data, family = stats::binomial),
 error = function(e) NULL)
 if (is.null(m)) next
 td <- broom::tidy(m, conf.int = TRUE, conf.level = conf, exponentiate = TRUE)
 rows_x <- td[td$term != "(Intercept)", ]
 if (nrow(rows_x) == 0) next
 best <- rows_x[which.min(rows_x$p.value), ]
 univ[[xv]] <- data.frame(
 variable = xv, OR = best$estimate,
 CI.low = best$conf.low, CI.high = best$conf.high,
 p.value = as.numeric(best$p.value), stringsAsFactors = FALSE)
 }
 if (!length(univ)) stop("No candidate variable could be fitted.")
 univ_df <- do.call(rbind, univ)
 rownames(univ_df) <- NULL
 univ_df$p.adjusted <- as.numeric(stats::p.adjust(univ_df$p.value,
 method = p_adjust))
 univ_df$selected <- univ_df$p.adjusted < alpha
 ## Screening now uses the adjusted p value; with p_adjust = "none" this is
 ## exactly the original behaviour.
 sig <- univ_df$variable[univ_df$selected]
 if (length(sig) > 0) {
 f2 <- stats::as.formula(paste(y, "~", paste(sig, collapse = " + ")))
 m2 <- stats::glm(f2, data = data, family = stats::binomial)
 mult <- broom::tidy(m2, conf.int = TRUE, conf.level = conf,
 exponentiate = TRUE)
 mult$variable <- mult$term
 } else {
 mult <- NULL
 message("Univariate screening selected no variables; multivariable model not built.")
 }
 stability <- NULL
 if (bootstrap > 0L) {
 if (!is.null(seed)) set.seed(seed)
 cnt <- setNames(rep(0L, length(x_vars)), x_vars)
 n <- nrow(data)
 for (b in seq_len(bootstrap)) {
 d <- data[sample.int(n, n, replace = TRUE), , drop = FALSE]
 for (xv in x_vars) {
 mm <- tryCatch(stats::glm(stats::as.formula(paste(y, "~", xv)),
 data = d, family = stats::binomial),
 error = function(e) NULL)
 if (is.null(mm)) next
 tt <- broom::tidy(mm)
 pp <- tt$p.value[tt$term != "(Intercept)"]
 if (length(pp) && min(pp, na.rm = TRUE) < alpha) cnt[xv] <- cnt[xv] + 1L
 }
 }
 stability <- data.frame(variable = x_vars,
 selection_frequency = as.numeric(cnt) / bootstrap,
 stringsAsFactors = FALSE)
 }
 list(univariate = univ_df, selected = sig, multivariable = mult,
 stability = stability, p_adjust = p_adjust)
}

#' Stepwise variable selection
#'
#' Stepwise selection based on AIC
#' (default \code{k = 2}), returning a coefficient table (OR / regression
#' Coefficients).
#'
#' @details
#' Purpose: Perform stepwise variable selection based on an information criterion
#' (AIC), automatically adding/removing variables from the full model to obtain a
#' Simpler regression model.
#' Method: Aims to minimise AIC and supports "both" (default, bidirectional, can
#' Enter or leave), "backward" (drop from the full model), and "forward" (add
#' From the null model); `k` is the penalty (default 2 gives AIC; using log(n)
#' Gives BIC).
#' Arguments: `data` is a data.frame; `formula` contains all candidate
#' Predictors; `direction` selects the direction; `k` is the AIC/BIC penalty.
#' Result: Returns a list with the final `coefficients` (coefficient table) and
#' `formula` (the formula of finally selected variables).
#' Note: Stepwise regression is convenient for exploration but has statistical
#' Flaws  -  p values are underestimated, confidence intervals are too narrow, and
#' Overfitting is likely; results also depend on arbitrary entry/removal
#' Thresholds. Treat it as a "draft" for variable screening; important
#' Conclusions still need domain knowledge and a pre-registered analysis plan.
#' @param data A \code{data.frame}.
#' @param formula Full model formula.
#' @param direction \code{"both"} (default), \code{"backward"}, or \code{"forward"}.
#' @param k AIC penalty; default \code{2}.
#' @return A list with \code{model} (final glm object), \code{coefficients}
#' (coefficient table), and \code{formula} (final formula).
#' @param family GLM family object or string; default `binomial()`. Previously hard-coded to logistic.
#' @param bootstrap Resamples for selection stability; `0` (default) skips it.
#' @param seed Random seed for the stability resamples.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' mc_stepwise(d, bp_control ~ age + bmi + ldl + treatment,
#'             direction = "both")
#' @export
mc_stepwise <- function(data, formula, direction = c("both", "backward", "forward"),
 k = 2, family = stats::binomial(), bootstrap = 0L, seed = NULL) {
 direction <- match.arg(direction)
 fit <- eval(substitute(stats::glm(formula, data = DATA, family = FAM),
 list(DATA = data, FAM = family)))
 if (direction == "forward") {
 null_fit <- stats::update(fit, . ~ 1)
 step_fit <- stats::step(null_fit,
 scope = list(lower = ~ 1, upper = formula),
 direction = "forward", k = k, trace = 0)
 } else {
 step_fit <- stats::step(fit, direction = direction, k = k, trace = 0)
 }
 fam_obj <- if (inherits(family, "family")) family else
 tryCatch(stats::binomial(), error = function(e) stats::gaussian())
 expo <- identical(fam_obj$family, "binomial") &&
 identical(fam_obj$link, "logit")
 td <- broom::tidy(step_fit, conf.int = TRUE, exponentiate = expo)
 td$variable <- td$term
 final_terms <- setdiff(all.vars(stats::formula(step_fit)[[3]]), ".")
 stability <- NULL
 if (bootstrap > 0L) {
 if (!is.null(seed)) set.seed(seed)
 cnt <- setNames(rep(0L, length(final_terms)), final_terms)
 n <- nrow(data)
 all_terms <- all.vars(formula[[3]])
 for (b in seq_len(bootstrap)) {
 d <- data[sample.int(n, n, replace = TRUE), , drop = FALSE]
 m <- tryCatch(stats::glm(formula, data = d, family = family),
 error = function(e) NULL)
 if (is.null(m)) next
 s <- tryCatch(stats::step(m, direction = "both", k = k, trace = 0),
 error = function(e) NULL)
 if (is.null(s)) next
 kept <- all.vars(stats::formula(s)[[3]])
 cnt[intersect(kept, names(cnt))] <- cnt[intersect(kept, names(cnt))] + 1L
 }
 stability <- data.frame(variable = final_terms,
 selection_frequency = as.numeric(cnt) / bootstrap,
 stringsAsFactors = FALSE)
 }
 list(model = step_fit, coefficients = td,
 formula = stats::formula(step_fit), family = fam_obj,
 stability = stability)
}

#' Ordinal logistic regression
#'
#' Based on
#' \pkg{MASS}'s \code{polr}; the outcome must be an ordered factor. Returns a
#' Coefficient table with OR.
#'
#' @details
#' Purpose: Fit an ordinal logistic regression (proportional-odds model),
#' Suitable when the outcome is an ordered multicategory variable
#' (e.g. mild/moderate/severe).
#' Method: Based on MASS's polr(), which assumes the predictor effects are the
#' Same across all cut points (proportional-odds assumption).
#' Arguments: `data` is a data.frame; `formula` is of the form
#' Ordered_outcome ~ predictors, where the outcome must be an ordered factor
#' (ordered()).
#' Result: Returns a coefficient table with coefficients (on the logit scale),
#' Standard errors, and the exponentiated OR (odds ratio).
#' Note: If the proportional-odds assumption fails (check via MASS::polr-based
#' Tests or large stratum-specific OR differences), interpret results cautiously;
#' Consider a partial proportional-odds model or a multinomial logistic model.
#' Always convert the outcome to an ordered factor first, otherwise the model is
#' Meaningless. Depends on MASS (Suggests).
#' @param data A \code{data.frame}.
#' @param formula Model formula; the outcome is an ordered factor.
#' @return A \code{data.frame} of coefficients with columns \code{variable} /
#' \code{Value} / \code{Std.Error} / \code{t.value} / \code{OR} /
#' \code{p.value}; \code{OR = exp(Value)}.
#' @param conf Confidence level for the OR interval; default `0.95`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' \donttest{
#' mc_ordinal_logistic(d, bp_grade ~ age + bmi + treatment)
#' }
#' @export
mc_ordinal_logistic <- function(data, formula, conf = 0.95) {
 if (!requireNamespace("MASS", quietly = TRUE)) {
 stop("`mc_ordinal_logistic` requires the MASS package (shipped with R).")
 }
 fit <- eval(substitute(MASS::polr(formula, data = DATA, Hess = TRUE),
 list(DATA = data)))
 co <- summary(fit)$coefficients
 cf <- as.data.frame(co)
 names(cf) <- c("estimate", "std_error", "t_value")
 cf$variable <- rownames(cf)
 zc <- stats::qnorm(1 - (1 - conf) / 2)
 cf$OR <- exp(cf$estimate)
 cf$OR_lower <- exp(cf$estimate - zc * cf$std_error)
 cf$OR_upper <- exp(cf$estimate + zc * cf$std_error)
 cf$p.value <- 2 * stats::pnorm(-abs(cf$t_value))
 cf$is_intercept <- grepl("\\|", cf$variable)
 ## Approximate proportional-odds check: refit a binary logit at every
 ## cut point and measure how far the shared slopes drift apart.
 mf <- stats::model.frame(formula, data = data)
 y <- stats::model.response(mf)
 lev <- levels(factor(y))
 po_stat <- NA_real_
 if (length(lev) >= 3L) {
 d <- mf[, -1, drop = FALSE]
 splits <- lapply(seq_len(length(lev) - 1L), function(i) {
 d$.y <- as.integer(as.integer(factor(y)) > i)
 m <- tryCatch(stats::glm(.y ~ ., data = d,
 family = stats::binomial),
 error = function(e) NULL)
 if (is.null(m)) NULL else stats::coef(m)
 })
 splits <- Filter(Negate(is.null), splits)
 if (length(splits) >= 2L) {
 common <- Reduce(intersect, lapply(splits, names))
 common <- setdiff(common, "(Intercept)")
 if (length(common)) {
 mm <- do.call(rbind, lapply(splits, function(v) v[common]))
 drift <- apply(mm, 2L, function(v) diff(range(v, na.rm = TRUE)))
 po_stat <- max(drift / pmax(0.05, cf$std_error[match(common, cf$variable)]),
 na.rm = TRUE)
 }
 }
 }
 checks <- .mc_check_row("Proportional odds (approximate)", po_stat,
 "max slope drift < 1.96 SE",
 if (is.na(po_stat)) NA else po_stat < 1.96,
 "binary logits at each cut point; screening approximation, not the Brant test")
 attr(cf, "checks") <- checks
 attr(cf, "fit") <- fit
 cf
}

#' Conditional logistic regression
#'
#' Suitable for
#' 1:M matched case-control studies; fits \code{survival::clogit()} which
#' Conditions on each matched set, yielding a conditional odds ratio.
#'
#' @details
#' Purpose: Analyse the association between exposure and outcome in matched
#' Case-control data, controlling for confounding from the matching variables.
#' Dependency: Requires the survival package (Suggests); a clear message is given
#' If not installed.
#' Method: Specify the matched stratum in the formula with
#' \code{strata(matching variable)}, and \code{survival::clogit()} performs a
#' Conditional logistic likelihood estimation within each stratum.
#' Arguments: \code{data} is a data.frame; \code{formula} is the model formula
#' With a 0/1 response (1 = case), and must include a
#' \code{strata(matching variable)} term.
#' Result: Returns a coefficient table with estimate, standard error, z value,
#' P value, OR (exp(estimate)), and the OR 95% confidence interval. OR > 1
#' Indicates the exposure is a risk factor.
#' Note: Each matched stratum should contain both cases and controls; if a
#' Stratum has only one status, it contributes nothing to the conditional
#' Likelihood. Model comparison can be done with mc_model_compare().
#' @param data A \code{data.frame}.
#' @param formula A formula with a \code{strata()} term.
#' @return A coefficient table with estimate / standard error / z / p / OR / OR
#' Confidence interval.
#' @param conf Confidence level for the OR interval; default `0.95`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' m <- d[1:600, ]                      # 300 simulated matched pairs
#' mc_clogistic(m,
#'   case_status ~ smoke + prior_cv + survival::strata(match_id))
#' @export
mc_clogistic <- function(data, formula, conf = 0.95) {
 if (!requireNamespace("survival", quietly = TRUE)) {
 stop("`mc_clogistic` requires the survival package; please install.packages('survival').")
 }
 # survival::clogit internally calls coxph() via eval(call, parent.frame()),
 # so it fails unless survival is on the search path.
 e <- new.env(parent = asNamespace("survival"))
 e$formula <- formula
 e$data <- data
 fit <- eval(call("clogit", formula = formula, data = data), envir = e)
 beta <- as.numeric(fit$coefficients)
 nam <- names(fit$coefficients)
 se <- sqrt(diag(as.matrix(fit$var)))
 z <- beta / se
 p <- 2 * stats::pnorm(-abs(z))
 zc <- stats::qnorm(1 - (1 - conf) / 2)
 ## Vectorised interval (the previous c(-1, 1) recycling was wrong for
 ## models with more than two coefficients).
 lo <- exp(beta - zc * se)
 hi <- exp(beta + zc * se)
 out <- data.frame(variable = nam, estimate = beta, std_error = se, z = z,
 p.value = p, OR = exp(beta), OR_lower = lo, OR_upper = hi,
 stringsAsFactors = FALSE)
 fstr <- paste(deparse(formula), collapse = " ")
 sm <- regmatches(fstr, regexpr("strata\\(([^)]+)\\)", fstr))
 strat_var <- if (length(sm)) sub("^strata\\(", "",
 sub("\\)$", "", sm[1])) else NA_character_
 n_strata <- if (!is.na(strat_var) && strat_var %in% names(data))
 length(unique(data[[strat_var]])) else NA_integer_
 scalar <- function(expr) {
 v <- tryCatch(as.numeric(expr), error = function(e) NA_real_)
 if (length(v) == 1L) v else NA_real_
 }
 events <- scalar(fit$nevent)
 n_obs <- scalar(if (is.null(fit$n)) nrow(data) else fit$n)
 lrt_p <- tryCatch({
 a <- stats::anova(fit)
 pp <- a[nrow(a), ]
 as.numeric(pp[grep("Pr", names(pp))[1]])
 }, error = function(e) NA_real_)
 checks <- data.frame(
 metric = c("n", "events", "n_strata", "logLik", "LRT_p"),
 value = c(n_obs, events, scalar(n_strata),
 scalar(stats::logLik(fit)), scalar(lrt_p)),
 stringsAsFactors = FALSE)
 attr(out, "checks") <- checks
 attr(out, "fit") <- fit
 out
}

#' Model comparison
#'
#' Computes AIC, BIC,
#' A *self-implemented* small-sample-corrected AICc and Akaike weights from AICc,
#' Then ranks models by weight.
#'
#' @details
#' Purpose: Compare candidate regression models by information criteria and
#' Select the relatively best one.
#' Method (original): AIC and BIC come from base R; AICc = AIC + 2k(k+1)/(n-k-1)
#' Where k is the number of parameters and n is the sample size (models with
#' N < k+1 are set to NA); then delta = AICc - min is computed, and the Akaike
#' Weight w = exp(-0.5*delta) / sum(exp(-0.5*delta)); a larger weight means a
#' Better model.
#' Arguments: \code{...} are several fitted model objects (e.g. lm/glm/coxph)
#' That support \code{AIC()}/\code{BIC()}/\code{nobs()}; \code{n} is the sample
#' Size (when omitted, it is attempted from \code{nobs()}).
#' Result: Returns a model-comparison table sorted by weight descending, with k
#' (number of parameters), n, AIC, BIC, AICc, delta_AICc, and weight.
#' Note: Models compared must be based on the same dataset and the same response
#' Variable, otherwise the criteria are not comparable. The weight reflects the
#' Relative probability of each model being the best, not an absolute
#' Goodness-of-fit.
#' @param ... Fitted model objects.
#' @param n Sample size; auto-detected from the model if NULL.
#' @return A ranked comparison table sorted by Akaike weight.
#' @param lrt If `TRUE`, add pairwise likelihood-ratio tests between nested models ordered by parameter count.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' m1 <- stats::lm(sbp_week12 ~ age + bmi, data = d)
#' m2 <- stats::lm(sbp_week12 ~ age + bmi + treatment, data = d)
#' mc_model_compare(m1, m2, lrt = TRUE)
#' @export
mc_model_compare <- function(..., n = NULL, lrt = FALSE) {
 models <- list(...)
 if (length(models) < 2) stop("`mc_model_compare` needs at least two models.")
 nm <- if (is.null(names(models))) {
 paste0("model", seq_along(models))
 } else {
 ifelse(nzchar(names(models)), names(models), paste0("model", seq_along(models)))
 }
 aic <- sapply(models, function(m) stats::AIC(m))
 bic <- sapply(models, function(m) stats::BIC(m))
 k <- sapply(models, function(m) length(stats::coef(m)))
 if (is.null(n)) {
 nob <- tryCatch(sapply(models, function(m) stats::nobs(m)), error = function(e) NULL)
 n <- if (is.null(nob)) stop("Could not auto-detect sample size; please pass n explicitly.") else nob[1]
 }
 aicc <- ifelse(n - k - 1 > 0, aic + 2 * k * (k + 1) / (n - k - 1), NA_real_)
 delta <- aicc - min(aicc, na.rm = TRUE)
 w <- if (any(is.na(aicc))) {
 rep(NA_real_, length(aicc))
 } else {
 e <- exp(-0.5 * delta)
 e / sum(e, na.rm = TRUE)
 }
 ## Comparability audit: information criteria are only meaningful for the
 ## same response fitted to the same observations.
 get_family <- function(m) tryCatch({
 if (!is.null(m$family)) m$family$family else class(m)[1]
 }, error = function(e) class(m)[1])
 get_response <- function(m) tryCatch(
 paste(deparse(stats::formula(m)[[2]]), collapse = ""),
 error = function(e) NA_character_)
 fam <- vapply(models, get_family, character(1))
 resp <- vapply(models, get_response, character(1))
 comparable <- length(unique(fam)) == 1L &&
 length(unique(resp[!is.na(resp)])) == 1L
 if (!comparable)
 warning("Models do not share the same family/response; AIC, BIC and AICc are not comparable.")
 lrt_p <- rep(NA_real_, length(models))
 if (isTRUE(lrt) && comparable) {
 ord <- order(k)
 for (i in seq_len(length(ord) - 1L)) {
 m1 <- models[[ord[i]]]; m2 <- models[[ord[i + 1L]]]
 q <- tryCatch({
 a <- if (inherits(m1, "glm") && inherits(m2, "glm"))
 stats::anova(m1, m2, test = "Chisq") else stats::anova(m1, m2)
 pp <- a[nrow(a), ]
 as.numeric(pp[grep("Pr", names(pp))[1]])
 }, error = function(e) NA_real_)
 lrt_p[ord[i + 1L]] <- q
 }
 }
 tab <- data.frame(model = nm, k = k, n = n, AIC = aic, BIC = bic,
 AICc = aicc, delta_AICc = delta, weight = w,
 family = fam, response = resp,
 comparable = comparable, lrt_p = lrt_p,
 stringsAsFactors = FALSE)
 tab[order(-tab$weight), ]
}
