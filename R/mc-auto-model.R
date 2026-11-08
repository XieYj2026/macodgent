# =============================================================================
# macodgent - unified model dispatcher
# One entry point for every modelling family in the package. The returned
# bundle has the same shape as mc_auto_analyze()$model, so it can be printed,
# reported and exported by the same code paths.
# =============================================================================

# Normalise a family-specific coefficient table into the shared schema.
.mc_norm_coef <- function(term, estimate, effect = "Beta", conf.low = NA_real_,
                          conf.high = NA_real_, p.value = NA_real_,
                          std.error = NA_real_, p_adjust = "BH",
                          alpha = 0.05) {
  n <- length(term)
  rep_n <- function(x) rep_len(as.numeric(x), n)
  out <- data.frame(
    term = as.character(term), effect = rep_len(as.character(effect), n),
    estimate = rep_n(estimate), std.error = rep_n(std.error),
    conf.low = rep_n(conf.low), conf.high = rep_n(conf.high),
    p.value = rep_n(p.value), stringsAsFactors = FALSE)
  out$p.adjusted <- as.numeric(stats::p.adjust(out$p.value, method = p_adjust))
  out$significant <- !is.na(out$p.value) & out$p.value < alpha
  out$n <- NA_integer_
  rownames(out) <- NULL
  out
}

# Normalise a family-specific checks table.
.mc_norm_checks <- function(x, family) {
  need <- c("check", "value", "threshold", "status", "note")
  if (is.data.frame(x) && all(c("check", "value", "status") %in% names(x))) {
    for (k in need) if (!k %in% names(x)) x[[k]] <- ""
    return(x[, need, drop = FALSE])
  }
  data.frame(check = paste0(family, " model fitted"), value = 1,
             threshold = "fit succeeds", status = "PASS",
             note = "no family-specific check panel available",
             stringsAsFactors = FALSE)
}

.mc_norm_perf <- function(...) {
  kv <- list(...)
  if (!length(kv)) return(data.frame(metric = "n", value = NA_real_))
  data.frame(metric = names(kv),
             value = vapply(kv, function(v) {
               z <- suppressWarnings(as.numeric(v)[1])
               if (length(z)) z else NA_real_
             }, numeric(1)), stringsAsFactors = FALSE)
}

#' Fit any macodgent model family with one call
#'
#' Dispatches to the package's modelling functions and returns a bundle with
#' the same shape as `mc_auto_analyze()$model`: `coefficients`, `performance`,
#' `checks`, `diagnostics`, `fit` and the raw family-specific result. Use it
#' directly for a single model, or through `mc_auto_analyze(family = ...)` to
#' combine it with the automated audit, descriptives, comparisons and report.
#'
#' @details
#' Supported families: `linear`, `logistic`, `cox`, `ordinal`, `quantile`,
#' `gam`, `rcs`, `pca`, `elastic_net`, `nnet`, `competing_risk` and
#' `multistate`. Families backed by packages in `Suggests` fail with a clear
#' message when the package is unavailable. `multistate` needs long-format
#' data and a transition matrix supplied via `trans = ...`.
#'
#' Predictors are inferred as all numeric columns except the outcome, time,
#' status and grouping columns; pass `predictors` to override. Rows with
#' missing values in the analysed columns are dropped (complete-case analysis).
#'
#' @param data A `data.frame`.
#' @param outcome Name of the outcome column.
#' @param time,status Column names for survival / competing-risk families.
#' @param family Model family; see Details.
#' @param predictors Optional character vector of predictor column names.
#' @param group Optional grouping column, used by the competing-risk family.
#' @param conf Confidence level; default `0.95`.
#' @param alpha Significance level; default `0.05`.
#' @param p_adjust Multiplicity method for coefficient p values; default `"BH"`.
#' @param seed Random seed for resampling-based families; default `1L`.
#' @param ... Extra arguments, e.g. `trans` for `family = "multistate"`.
#' @return A list of class `mc_auto_model` with `family`, `formula`,
#'   `predictors`, `n`, `events`, `coefficients`, `performance`, `checks`,
#'   `diagnostics`, `fit` and `raw`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' m <- mc_auto_model(d, outcome = "sbp_week12", family = "linear")
#' m$performance
#' m$checks
#' print(m)
#' @export
mc_auto_model <- function(data, outcome = NULL, time = NULL, status = NULL,
                          family = c("linear", "logistic", "cox", "ordinal",
                                     "quantile", "gam", "rcs", "pca",
                                     "elastic_net", "nnet",
                                     "competing_risk", "multistate"),
                          predictors = NULL, group = NULL, conf = 0.95,
                          alpha = 0.05, p_adjust = "BH", seed = 1L, ...) {
  family <- match.arg(family)
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  extra <- list(...)
  if (family %in% c("cox", "competing_risk") &&
      (is.null(time) || is.null(status)))
    stop("family '", family, "' requires both `time` and `status`.")
  if (!family %in% c("cox", "competing_risk", "multistate") && is.null(outcome))
    stop("family '", family, "' requires `outcome`.")
  if (is.null(predictors)) {
    num <- names(data)[vapply(data, is.numeric, logical(1))]
    predictors <- setdiff(num, c(outcome, time, status, group))
  }
  bad <- setdiff(c(outcome, time, status, group, predictors), names(data))
  if (length(bad)) stop("Unknown column(s): ", paste(bad, collapse = ", "))
  if (!length(predictors)) stop("No numeric predictors available.")

  need <- unique(c(outcome, time, status, group, predictors))
  d <- data[, need, drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) < 10L) stop("At least 10 complete rows are required.")

  bt <- function(x) paste0("`", gsub("`", "", x), "`")
  rhs <- paste(vapply(predictors, bt, character(1)), collapse = " + ")
  lhs <- if (!is.null(outcome)) bt(outcome) else ""
  fml <- switch(family,
    cox = stats::as.formula(paste0("survival::Surv(", bt(time), ", ",
                                   bt(status), ") ~ ", rhs)),
    multistate = stats::as.formula(paste0(
      "survival::Surv(", bt(extra$start %||% "Tstart"), ", ",
      bt(extra$stop %||% "Tstop"), ", ", bt(extra$event %||% "status"),
      ") ~ ", rhs)),
    stats::as.formula(paste(lhs, "~", rhs)))
  fml_txt <- paste(deparse(fml), collapse = " ")

  co <- NULL; checks <- NULL; perf <- NULL; fit <- NULL
  raw <- NULL; diagnostics <- NULL; events <- NA_integer_

  if (family %in% c("linear", "logistic", "cox")) {
    fit <- if (family == "cox") .mc_fit(fml, d, "cox") else
      if (family == "linear") stats::lm(fml, data = d) else
        stats::glm(fml, data = d, family = stats::binomial)
    co <- .mc_coef(fit, family, conf, p_adjust, alpha)
    checks <- .mc_norm_checks(.mc_model_checks(fit, family, alpha), family)
    perf <- .mc_norm_perf(n = nrow(d), AIC = stats::AIC(fit),
                          BIC = stats::BIC(fit),
                          logLik = as.numeric(stats::logLik(fit)))
    if (family == "cox") events <- sum(d[[status]] == 1)
    if (family == "logistic") events <- {
      y <- d[[outcome]]
      tb <- table(if (is.factor(y)) as.character(y) else y)
      if (length(tb) >= 2L) as.integer(min(tb)) else 0L
    }
    diagnostics <- list(fit = fit)
  } else if (family == "ordinal") {
    y <- d[[outcome]]
    if (!is.ordered(y)) d[[outcome]] <- ordered(y)
    raw <- mc_ordinal_logistic(d, fml)
    k <- !raw$is_intercept
    co <- .mc_norm_coef(raw$variable[k], raw$OR[k], "OR",
                        raw$OR_lower[k], raw$OR_upper[k], raw$p.value[k],
                        raw$std_error[k], p_adjust, alpha)
    checks <- .mc_norm_checks(attr(raw, "checks"), family)
    perf <- .mc_norm_perf(n = nrow(d), n_levels = length(levels(d[[outcome]])))
    diagnostics <- attr(raw, "checks")
  } else if (family == "quantile") {
    raw <- mc_quantile_reg(d, fml, taus = 0.5, conf = conf)
    co <- .mc_norm_coef(raw$term, raw$estimate, "Beta",
                        raw$conf.low, raw$conf.high, raw$p.value,
                        raw$std_error, p_adjust, alpha)
    checks <- .mc_norm_checks(NULL, family)
    perf <- .mc_norm_perf(n = nrow(d), tau = 0.5,
                          n_coefficients = nrow(raw))
    diagnostics <- attr(raw, "trend")
  } else if (family == "gam") {
    if (!requireNamespace("mgcv", quietly = TRUE))
      stop("family 'gam' requires the mgcv package.")
    s <- mgcv::s
    ## Smooth only genuinely continuous columns; low-cardinality scores
    ## cannot support a spline basis and make mgcv fail on the degrees of
    ## freedom, so they enter as plain linear terms.
    cont <- predictors[vapply(d[predictors], function(v)
      length(unique(v)) > 10L, logical(1))]
    lin <- setdiff(predictors, cont)
    gam_terms <- c(
      vapply(cont, function(p) paste0("s(", bt(p), ", k = 5)"), character(1)),
      vapply(lin, bt, character(1)))
    if (!length(gam_terms)) gam_terms <- "1"
    gam_fml <- stats::as.formula(paste(lhs, "~",
                                       paste(gam_terms, collapse = " + ")))
    raw <- mc_gam(d, gam_fml)
    fml_txt <- paste(deparse(gam_fml), collapse = " ")
    td <- raw$terms
    ## broom::tidy(gam) reports edf / ref.df for smooth terms instead of
    ## estimate / std.error, so fill the missing columns explicitly.
    pick <- function(k) if (k %in% names(td)) as.numeric(td[[k]]) else
      rep(NA_real_, nrow(td))
    co <- .mc_norm_coef(td$term, pick("estimate"), "Beta",
                        rep(NA_real_, nrow(td)), rep(NA_real_, nrow(td)),
                        pick("p.value"), pick("std.error"), p_adjust, alpha)
    checks <- .mc_norm_checks(raw$checks, family)
    perf <- raw$performance
    diagnostics <- list(smooths = raw$smooths)
  } else if (family == "rcs") {
    x <- predictors[1L]
    binary <- length(unique(d[[outcome]])) == 2L
    raw <- mc_rcs(d, x = x, y = outcome, nk = 3L,
                  family = if (binary) "binomial" else "gaussian",
                  conf = conf)
    td <- raw$coefficients
    co <- .mc_norm_coef(td$term, td$estimate,
                        if (binary) "OR" else "Beta",
                        td$conf.low, td$conf.high, td$p.value,
                        td$std.error, p_adjust, alpha)
    checks <- .mc_norm_checks(raw$checks, family)
    perf <- .mc_norm_perf(n = nrow(d), knots = length(raw$knots),
                          nonlinearity_p = raw$nonlinearity_p)
    diagnostics <- list(curve = raw$curve, knots = raw$knots)
    fml_txt <- paste0(lhs, " ~ rcs(", bt(x), ")")
  } else if (family == "pca") {
    raw <- mc_pca_regression(d, fml, seed = seed, conf = conf,
                             validation = "holdout")
    k <- raw$coef_original$term != "(Intercept)"
    cf <- raw$coef_original[k, ]
    co <- .mc_norm_coef(cf$term, cf$estimate, "Beta", cf$conf.low,
                        cf$conf.high, cf$p.value, cf$std.error,
                        p_adjust, alpha)
    checks <- .mc_norm_checks(NULL, family)
    perf <- raw$performance
    diagnostics <- list(external = raw$external, cv_table = raw$cv_table)
  } else if (family == "elastic_net") {
    raw <- mc_elastic_net(d, x = predictors, y = outcome, seed = seed)
    cf <- raw$coefficients
    co <- .mc_norm_coef(cf$variable, cf$coefficient, "Beta",
                        p.value = rep(NA_real_, nrow(cf)))
    checks <- .mc_norm_checks(NULL, family)
    perf <- raw$performance
    diagnostics <- list(cv_table = raw$cv_table)
    fml_txt <- paste0(lhs, " ~ penalised(", paste(predictors, collapse = " + "), ")")
  } else if (family == "nnet") {
    raw <- mc_nnet_reg(d, fml, seed = seed, validation = "holdout")
    ## A neural network has no inferential coefficient table; report the
    ## out-of-sample performance instead of pretending to have p values.
    co <- NULL
    checks <- .mc_norm_checks(NULL, family)
    perf <- raw$performance
    diagnostics <- list(validation = raw$validation,
                        predictions = raw$predictions)
  } else if (family == "competing_risk") {
    if (is.null(group)) group <- predictors[1L]
    raw <- mc_competing_risk(d, time = time, status = status,
                             group = group, conf = conf)
    cf <- raw$coefficients
    co <- .mc_norm_coef(cf$variable, cf$HR, "SHR", cf$HR_lower,
                        cf$HR_upper, cf$p.value, cf$se, p_adjust, alpha)
    checks <- .mc_norm_checks(NULL, family)
    perf <- raw$checks
    diagnostics <- list(cif = raw$cif, cause_counts = raw$cause_counts)
    events <- sum(d[[status]] != 0)
    fml_txt <- paste0("Fine-Gray(", time, ", ", status, ") ~ ", group)
  } else if (family == "multistate") {
    if (is.null(extra$trans))
      stop("family 'multistate' requires a transition matrix via `trans = ...` and long-format data.")
    raw <- mc_multistate(d, trans = extra$trans, formula = fml)
    td <- tryCatch(as.data.frame(broom::tidy(raw$model)),
                   error = function(e) NULL)
    if (!is.null(td))
      co <- .mc_norm_coef(td$term, td$estimate, "HR",
                          rep(NA_real_, nrow(td)), rep(NA_real_, nrow(td)),
                          td$p.value, td$std.error, p_adjust, alpha)
    checks <- .mc_norm_checks(NULL, family)
    perf <- raw$checks
    diagnostics <- list(probability_table = raw$probability_table)
  }

  if (is.null(co))
    checks <- rbind(checks,
                    .mc_check_row("Inferential coefficients", 0,
                                  "coefficient table available", FALSE,
                                  "this family is predictive, not inferential"))
  out <- list(family = family, formula = fml_txt,
              predictors = predictors, n = nrow(d), events = events,
              coefficients = co, performance = perf,
              checks = checks, diagnostics = diagnostics,
              fit = fit, raw = raw)
  class(out) <- c("mc_auto_model", "list")
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' @export
print.mc_auto_model <- function(x, ...) {
  cat("macodgent model family:", x$family, "\n")
  cat("n =", x$n, "| predictors =", length(x$predictors), "\n")
  cat(x$formula, "\n")
  if (!is.null(x$coefficients)) {
    cat("\nCoefficients\n")
    print(x$coefficients[, c("term", "effect", "estimate", "p.value")],
          row.names = FALSE, digits = 3)
  }
  cat("\nChecks\n")
  print(x$checks[, c("check", "value", "status")],
        row.names = FALSE, digits = 3)
  invisible(x)
}
