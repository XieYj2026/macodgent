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
#' @examples
#' # multivariable logistic regression for low birth weight
#' birthwt <- MASS::birthwt
#' mc_logistic(birthwt, low ~ age + lwt + smoke + factor(race))
#' @export
mc_logistic <- function(data, formula, conf = 0.95) {
 fit <- eval(substitute(stats::glm(formula, data = DATA,
 family = stats::binomial),
 list(DATA = data)))
 td <- broom::tidy(fit, conf.int = TRUE, conf.level = conf, exponentiate = TRUE)
 td
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
#' @examples
#' # : univariate screening then multivariable modelling ("screen then model")
#' birthwt <- MASS::birthwt
#' bw <- birthwt
#' bw$race <- factor(bw$race)
#' bw$smoke <- factor(bw$smoke)
#' res <- mc_logistic_screen(bw, y = "low",
#' x_vars = c("age","lwt","race","smoke","ht","ui"))
#' cat("==== univariate ====\n"); print(res$univariate, row.names = FALSE)
#' cat("\n==== selected ====\n"); print(res$selected)
#' cat("\n==== multivariable ====\n"); print(as.data.frame(res$multivariable), row.names = FALSE)
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
#' @examples
#' # aIC-based stepwise variable selection
#' birthwt <- MASS::birthwt
#' bw <- birthwt
#' bw$race <- factor(bw$race)
#' bw$smoke <- factor(bw$smoke)
#' mc_stepwise(bw, low ~ age + lwt + race + smoke + ht + ui)$coefficients
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
#' @examples
#' \donttest{
#' # ordinal logistic  -  ordinal outcome of number of prenatal visits (ftv)
#' birthwt <- MASS::birthwt
#' bw <- birthwt
#' bw$ftv_ord <- ordered(cut(bw$ftv, c(-1, 0, 2, Inf),
#' labels = c("0", "1-2", "3+")))
#' bw$smoke <- factor(bw$smoke)
#' mc_ordinal_logistic(bw, ftv_ord ~ age + lwt + smoke)
#' }
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
#' @examples
#' # demonstrating conditional logistic regression with the built-in matched case-control dataset infert
#' # case: infertility status (1 = case); stratum: matched set; spontaneous/induced: miscarriage history
#' infert <- datasets::infert
#' mc_clogistic(infert, case ~ spontaneous + induced + survival::strata(stratum))
#' @export
mc_clogistic <- function(data, formula) {
 if (!requireNamespace("survival", quietly = TRUE)) {
 stop("`mc_clogistic` requires the survival package; please install.packages('survival').")
 }
 # survival::clogit internally calls coxph() via eval(call, parent.frame()),
 # so it fails unless survival is on the search path. Evaluate it from within a
 # child environment of the survival namespace so coxph resolves correctly.
 e <- new.env(parent = asNamespace("survival"))
 e$formula <- formula
 e$data <- data
 fit <- eval(call("clogit", formula = formula, data = data), envir = e)
 beta <- as.numeric(fit$coefficients)
 nam <- names(fit$coefficients)
 se <- sqrt(diag(as.matrix(fit$var)))
 z <- beta / se
 p <- 2 * stats::pnorm(-abs(z))
 or <- exp(beta)
 ci <- exp(beta + c(-1, 1) * stats::qnorm(0.975) * se)
 data.frame(variable = nam, estimate = beta, std_error = se, z = z,
 p.value = p, OR = or, OR_lower = ci[1], OR_upper = ci[2],
 stringsAsFactors = FALSE)
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
#' @param... Fitted model objects.
#' @param n Sample size; auto-detected from the model if NULL.
#' @return A ranked comparison table sorted by Akaike weight.
#' @examples
#' # comparing multiple linear models predicting birth weight (aIC/bIC/aICc/weights)
#' birthwt <- MASS::birthwt
#' m1 <- lm(bwt ~ lwt, data = birthwt)
#' m2 <- lm(bwt ~ lwt + age, data = birthwt)
#' m3 <- lm(bwt ~ lwt + age + smoke + factor(race), data = birthwt)
#' mc_model_compare(m1, m2, m3)
#' @export
mc_model_compare <- function(..., n = NULL) {
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
 tab <- data.frame(model = nm, k = k, n = n, AIC = aic, BIC = bic,
 AICc = aicc, delta_AICc = delta, weight = w,
 stringsAsFactors = FALSE)
 tab[order(-tab$weight), ]
}
