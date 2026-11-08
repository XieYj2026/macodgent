#' Run a complete automated analysis on one data frame
#'
#' One call produces an auditable bundle of results: a per-column audit with
#' automatic role assignment, descriptive statistics with distribution checks,
#' a stratified baseline table, group comparisons with effect sizes and
#' multiplicity-adjusted p values, a correlation analysis, univariable
#' screening, a multivariable model with an explicit diagnostic panel, written
#' findings and a log of every automatic decision. Set `report` to write a
#' self-contained HTML report at the same time.
#'
#' @details
#' Purpose: remove the repetitive glue work between loading a cohort and having
#' a defensible first-pass analysis, without hiding any of the choices made.
#'
#' Automatic decisions (all recorded in `decisions`):
#' \itemize{
#' \item Columns that are constant, mostly missing, identifier-like, of an
#' unsupported type, or high-cardinality categorical variables are excluded.
#' Numeric 0/1 columns other than the outcome are treated as two-level factors.
#' \item The model family is inferred from the outcome: continuous outcomes use
#' linear regression, two-level outcomes use logistic regression, and supplying
#' `time` plus `status` selects a Cox model. `status` is recoded automatically
#' when it is not already 0/1 (the largest code, or the last factor level, is
#' taken as the event).
#' \item Pearson correlation is used when all numeric columns look normal,
#' otherwise Spearman. Numeric predictors correlated at or above
#' `collinear_cut` are pruned greedily, keeping the variable with the stronger
#' univariable association first.
#' \item The number of predictors is capped by `events / epv_min` for binary
#' and survival outcomes and by `n / epv_min` for continuous outcomes, unless
#' `max_predictors` is supplied. When `screening = TRUE`, variables are ranked
#' by their univariable p value before the cap is applied. A fit that fails or
#' yields a rank-deficient design is retried with fewer predictors.
#' }
#'
#' Diagnostics are family specific: linear models get residual normality,
#' homoscedasticity, collinearity, influence and global fit checks; logistic
#' models reuse [mc_logistic_report()] for discrimination, calibration, EPV,
#' VIF and influence; Cox models get proportional hazards, concordance, EPV
#' and global fit checks.
#'
#' Note: automated test selection and variable pre-selection are exploratory.
#' They do not replace a prospectively specified analysis plan, clinical
#' confounder selection, or inspection of the raw distributions. Every
#' excluded column, recoded variable and dropped predictor is listed in
#' `decisions`; `warnings` collects the methodological caveats.
#'
#' @param data A `data.frame`.
#' @param outcome Name of the outcome column. Continuous values give a linear
#'   model, two-level values a logistic model. `NULL` skips regression.
#' @param group Name of a grouping column for the baseline table and group
#'   comparisons. If `NULL` and `outcome` is binary, the outcome is used.
#' @param time,status Column names for survival analysis. Supplying both
#'   overrides `outcome` and fits a Cox model.
#' @param vars Optional character vector of columns to analyse; default all.
#' @param max_missing Maximum proportion of missing values allowed in a
#'   predictor; default `0.5`.
#' @param collinear_cut Absolute correlation at which one of a pair of numeric
#'   predictors is dropped; default `0.8`.
#' @param p_adjust Multiplicity method used for all families of tests.
#' @param alpha Significance level; default `0.05`.
#' @param conf Confidence level; default `0.95`.
#' @param epv_min Minimum events (or observations) per predictor used to cap
#'   model size; default `10`.
#' @param max_predictors Optional hard cap on the number of predictors.
#' @param screening Rank predictors by univariable p value before applying the
#'   cap; default `TRUE`.
#' @param auc_boot Number of bootstrap resamples for the logistic AUC interval;
#'   default `200`. Use `0` for the analytic approximation.
#' @param report Optional path to an HTML report written on the same call.
#' @param seed Seed for bootstrap and simulated tests; default `1L`.
#' @param verbose Print each automatic decision as it is made; default `FALSE`.
#' @param family Model family. `"auto"` (default) infers linear, logistic or
#'   Cox from the outcome. Any other value is routed through
#'   [mc_auto_model()], which supports ordinal, quantile, GAM, RCS, PCR,
#'   elastic net, neural-network, competing-risk, multistate and the three
#'   classical families. The automated audit, descriptives, baseline table,
#'   group comparisons and correlation analysis still run in that case.
#' @return An object of class `mc_auto_analysis`: a list with `overview`,
#'   `columns`, `descriptives`, `baseline`, `comparison`, `correlation`,
#'   `screening`, `model`, `findings`, `decisions`, `warnings`, `data` and
#'   `call`. Tabular components are `data.frame` objects.
#' @seealso \code{\link{mc_auto_report}}, \code{\link{mc_compare_groups}},
#'   \code{\link{mc_logistic_report}}, \code{\link{mc_baseline_table}}.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' auto <- mc_auto_analyze(d, outcome = "bp_control", group = "treatment",
#'   vars = c("age", "bmi", "sbp_baseline", "ldl", "glucose", "smoke", "diabetes"),
#'   auc_boot = 0L)
#' auto$findings
#' auto$model$checks
#' auto$decisions
#' @export
mc_auto_analyze <- function(data, outcome = NULL, group = NULL,
                            time = NULL, status = NULL, vars = NULL,
                            max_missing = 0.5, collinear_cut = 0.8,
                            p_adjust = c("holm", "bonferroni", "BH", "BY", "none"),
                            alpha = 0.05, conf = 0.95, epv_min = 10,
                            max_predictors = NULL, screening = TRUE,
                            auc_boot = 200L, report = NULL, seed = 1L,
                            verbose = FALSE, family = "auto") {
  p_adjust <- match.arg(p_adjust)
  family_arg <- family
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (nrow(data) < 5L) stop("`data` needs at least 5 rows.")
  if (xor(is.null(time), is.null(status)))
    stop("Supply both `time` and `status`, or neither.")
  if (!is.null(vars)) {
    bad <- setdiff(vars, names(data))
    if (length(bad)) stop("Unknown column(s): ", paste(bad, collapse = ", "))
    # `vars` restricts candidate predictors only. Outcome, group, time and
    # status are structural columns and must survive the subsetting step.
    keep <- unique(c(vars, outcome, group, time, status))
    data <- data[, keep, drop = FALSE]
  }
  special <- c(outcome, group, time, status)
  bad <- setdiff(special, names(data))
  if (length(bad)) stop("Unknown column(s): ", paste(bad, collapse = ", "))

  dec <- new.env(parent = emptyenv())
  dec$rows <- list()
  add_dec <- function(step, decision, reason = "") {
    dec$rows[[length(dec$rows) + 1L]] <- data.frame(
      step = step, decision = decision, reason = reason,
      stringsAsFactors = FALSE)
    if (verbose) cat("[mc_auto] ", step, ": ", decision, "\n", sep = "")
  }
  warns <- character(0)
  add_warn <- function(w) warns <<- c(warns, w)

  n0 <- nrow(data)

  ## ---- model family -------------------------------------------------------
  family <- "none"
  if (!is.null(time)) {
    family <- "cox"
  } else if (!is.null(outcome)) {
    kx <- .mc_kind(data[[outcome]])
    nl <- length(unique(data[[outcome]][!is.na(data[[outcome]])]))
    if (kx == "numeric") family <- "linear"
    else if (kx == "binary" || (kx == "categorical" && nl == 2L))
      family <- "logistic"
    else {
      add_warn(sprintf(
        "Outcome '%s' has %d levels / unsupported type; regression skipped.",
        outcome, nl))
    }
  }
  add_dec("Model family", family,
          if (family == "cox") "time + status supplied" else
            if (family == "none") "no usable outcome supplied" else
              paste0("inferred from outcome '", outcome, "'"))

  ## ---- grouping -----------------------------------------------------------
  if (is.null(group) && family == "logistic") {
    group <- outcome
    add_dec("Grouping", paste0("group = outcome '", outcome, "'"),
            "a binary outcome is the conventional Table One stratifier")
  } else if (!is.null(group)) {
    add_dec("Grouping", paste0("group = '", group, "'"),
            "used for the baseline table and comparisons")
  }

  ## ---- explicit family override / extended model routing ------------------
  extended_model <- NULL
  if (!identical(family_arg, "auto")) {
    allowed <- c("linear", "logistic", "cox", "ordinal", "quantile", "gam",
                 "rcs", "pca", "elastic_net", "nnet", "competing_risk",
                 "multistate")
    fam_choice <- match.arg(family_arg, allowed)
    if (!fam_choice %in% c("linear", "logistic", "cox")) {
      extended_model <- mc_auto_model(data, outcome = outcome, time = time,
                                      status = status, family = fam_choice,
                                      group = group, conf = conf, alpha = alpha,
                                      p_adjust = p_adjust, seed = seed)
      family <- fam_choice
      add_dec("Model family (explicit)", fam_choice,
              "routed through mc_auto_model(); audit, descriptives and comparisons still run")
      add_dec("Extended model",
              sprintf("%s model fitted with %d predictor(s), n = %d",
                      fam_choice, length(extended_model$predictors),
                      extended_model$n),
              extended_model$formula)
    } else {
      if (fam_choice == "cox" && is.null(time))
        stop("family = 'cox' requires both `time` and `status`.")
      family <- fam_choice
      add_dec("Model family (explicit)", fam_choice,
              "user-specified family overrides outcome inference")
    }
  }

  ## ---- binary numeric columns become factors ------------------------------
  keep_numeric <- c(outcome, time, status)
  for (nm in names(data)) {
    if (nm %in% keep_numeric) next
    x <- data[[nm]]
    if (is.numeric(x) && .mc_kind(x) == "binary") {
      lv <- sort(unique(x[!is.na(x)]))
      data[[nm]] <- factor(as.character(x), levels = as.character(lv))
      add_dec("Type conversion", paste0(nm, " -> factor"),
              "numeric column with two distinct values treated as a binary factor")
    }
  }
  if (!is.null(group) && is.factor(data[[group]]))
    data[[group]] <- droplevels(data[[group]])

  ## ---- column audit -------------------------------------------------------
  audit <- do.call(rbind, lapply(names(data), function(nm) {
    x <- data[[nm]]
    kind <- .mc_kind(x)
    n_miss <- sum(is.na(x))
    u <- length(unique(x[!is.na(x)]))
    pct <- n_miss / nrow(data)
    role <- "predictor"
    reason <- "eligible"
    if (kind %in% c("date", "other")) {
      role <- "excluded"; reason <- "unsupported column type"
    } else if (u == 0L) {
      role <- "excluded"; reason <- "all values missing"
    } else if (u == 1L) {
      role <- "excluded"; reason <- "constant column"
    } else if (pct > max_missing) {
      role <- "excluded"
      reason <- sprintf("missing %.1f%% exceeds %.1f%% limit", 100 * pct,
                        100 * max_missing)
    } else if (u >= 10L && u > max(10L, min(20L, floor(n0 / 10))) &&
               kind == "categorical") {
      role <- "excluded"; reason <- "high-cardinality categorical column"
    } else if ((n0 - n_miss) >= 20L && u / (n0 - n_miss) > 0.95 &&
               (kind == "categorical" ||
                grepl("(^id$|_id$|id$|patient|subject|record|number|serial|\\u7f16\\u53f7)",
                      tolower(nm)))) {
      role <- "excluded"; reason <- "identifier-like (nearly unique values)"
    }
    if (!is.null(outcome) && nm == outcome) {
      role <- "outcome"; reason <- "supplied as outcome"
    } else if (!is.null(time) && nm == time) {
      role <- "time"; reason <- "supplied as survival time"
    } else if (!is.null(status) && nm == status) {
      role <- "status"; reason <- "supplied as survival status"
    } else if (!is.null(group) && nm == group) {
      role <- "group"; reason <- "supplied as grouping variable"
    }
    data.frame(variable = nm, class = paste(class(x), collapse = ","),
               type = kind, n_missing = n_miss,
               pct_missing = round(100 * pct, 1), n_unique = u,
               role = role, reason = reason, in_model = FALSE,
               stringsAsFactors = FALSE)
  }))
  rownames(audit) <- NULL
  excluded <- audit$variable[audit$role == "excluded"]
  if (length(excluded)) add_dec("Column audit",
                                paste0("excluded ", length(excluded), ": ",
                                       paste(utils::head(excluded, 8L),
                                             collapse = ", ")),
                                "see `columns$reason` for each column")

  ## ---- descriptives -------------------------------------------------------
  desc <- do.call(rbind, lapply(names(data), function(nm)
    .mc_desc_row(nm, data[[nm]], alpha)))
  rownames(desc) <- NULL

  normal_all <- length(desc$type[desc$type == "numeric"]) > 0L &&
    all(desc$distribution[desc$type == "numeric"] == "approximately normal")
  cor_method <- if (normal_all) "pearson" else "spearman"

  ## ---- baseline table and group comparison -------------------------------
  cmp_vars <- audit$variable[audit$role == "predictor"]
  baseline <- NULL
  comparison <- NULL
  if (!is.null(group) && length(cmp_vars)) {
    baseline <- suppressWarnings(tryCatch(
      mc_baseline_table(data, vars = cmp_vars, group = group,
                        p_value = TRUE, smd = TRUE, p_adjust = p_adjust),
      error = function(e) NULL))
    comparison <- suppressWarnings(tryCatch(
      mc_compare_groups(data, group = group, vars = cmp_vars,
                        p_adjust = p_adjust, conf = conf, alpha = alpha,
                        seed = seed),
      error = function(e) NULL))
    if (is.null(comparison))
      add_warn("Group comparison could not be computed; check the grouping column.")
  }

  ## ---- correlation --------------------------------------------------------
  num_vars <- audit$variable[audit$role %in% c("predictor", "group", "outcome") &
                               audit$type == "numeric"]
  correlation <- NULL
  if (length(num_vars) >= 2L) {
    correlation <- tryCatch(
      mc_correlate(data, vars = num_vars, method = cor_method,
                   conf = conf, p_adjust = p_adjust),
      error = function(e) NULL)
    add_dec("Correlation", paste0(cor_method, " on ", length(num_vars),
                                  " numeric columns"),
            if (normal_all) "all numeric columns look approximately normal" else
              "at least one numeric column is non-normal or untested")
  }
  strong <- if (!is.null(correlation))
    correlation[!is.na(correlation$r) & abs(correlation$r) >= collinear_cut, ,
                drop = FALSE] else NULL

  ## ---- modelling ----------------------------------------------------------
  screening_tab <- NULL
  model <- NULL
  if (!is.null(extended_model)) {
    model <- extended_model
  } else if (family == "none") {
    if (is.null(outcome) && is.null(time))
      add_dec("Regression", "skipped",
              "no outcome supplied; pass outcome = '...' or time/status to model")
  } else {
    predictors <- audit$variable[audit$role %in% c("predictor", "group") &
                                   audit$reason != "unsupported column type"]
    predictors <- setdiff(predictors, c(outcome, time, status))
    if (!length(predictors)) {
      add_warn("No eligible predictors remain after the column audit.")
      family <- "none"
    } else {
      need <- unique(c(outcome, time, status, predictors))
      md <- data[, need, drop = FALSE]
      md <- md[stats::complete.cases(md), , drop = FALSE]
      for (nm in names(md))
        if (is.factor(md[[nm]])) md[[nm]] <- droplevels(md[[nm]])

      ## survival status recoding
      if (family == "cox") {
        st <- md[[status]]
        if (is.factor(st) || is.character(st)) {
          lv <- levels(factor(st))
          md[[status]] <- as.integer(as.character(st) == lv[length(lv)])
          add_dec("Survival coding", paste0("event = level '",
                                            lv[length(lv)], "'"),
                  "factor status recoded to 0/1")
        } else if (!is.logical(st)) {
          u <- sort(unique(as.numeric(st)))
          if (!identical(u, c(0, 1))) {
            ev <- max(u)
            md[[status]] <- as.integer(as.numeric(st) == ev)
            add_dec("Survival coding", paste0("event = ", ev),
                    "status was not 0/1; largest code treated as the event")
          }
        } else {
          md[[status]] <- as.integer(st)
        }
      }
      lhs <- if (family == "cox")
        sprintf("survival::Surv(%s, %s)", .mc_bt(time), .mc_bt(status)) else
        .mc_bt(outcome)

      events <- if (family == "cox") sum(md[[status]] == 1) else
        if (family == "logistic") {
          y <- md[[outcome]]
          tb <- table(if (is.factor(y)) as.character(y) else y)
          if (length(tb) < 2L) 0L else min(tb)
        } else nrow(md)
      n_model <- if (family == "linear") nrow(md) else as.integer(events)

      ## Outcome-leakage guard: a predictor that alone separates the outcome
      ## almost perfectly is usually derived from it (e.g. birth weight and a
      ## low-birth-weight indicator). Keeping it would produce a meaningless
      ## model with an AUC of 1, so it is excluded before any model is fitted.
      y_leak <- if (family == "cox") as.numeric(md[[status]]) else {
        yy <- md[[outcome]]
        if (is.factor(yy) || is.character(yy))
          as.integer(factor(yy)) - 1L else as.numeric(yy)
      }
      leaky <- character(0)
      for (v in predictors) {
        xx <- md[[v]]
        if (!is.numeric(xx)) next
        a <- if (family == "linear") {
          suppressWarnings(abs(stats::cor(xx, y_leak,
                                          use = "complete.obs")))
        } else {
          au <- tryCatch(.mc_auc(as.numeric(xx), y_leak),
                         error = function(e) NA_real_)
          if (is.na(au)) NA_real_ else max(au, 1 - au)
        }
        if (is.finite(a) && a >= 0.98) leaky <- c(leaky, v)
      }
      if (length(leaky)) {
        predictors <- setdiff(predictors, leaky)
        add_dec("Outcome leakage",
                paste0("excluded ", paste(leaky, collapse = ", ")),
                "single-variable association >= 0.98; likely derived from the outcome")
        add_warn(sprintf(
          "Possible outcome leakage: %s predicts the outcome almost perfectly and was excluded; verify it is not a derived variable.",
          paste(leaky, collapse = ", ")))
      }

      ## univariable screening
      if (screening && nrow(md) >= 5L) {
        rows <- lapply(predictors, function(v) {
          f <- stats::as.formula(paste(lhs, "~", .mc_bt(v)))
          fit1 <- suppressWarnings(tryCatch(.mc_fit(f, md, family),
                                            error = function(e) NULL))
          if (is.null(fit1)) return(NULL)
          co <- .mc_coef(fit1, family, conf, "none", alpha)
          if (is.null(co) || !nrow(co)) return(NULL)
          data.frame(variable = v, effect = co$effect[1],
                     estimate = co$estimate[1], conf.low = co$conf.low[1],
                     conf.high = co$conf.high[1],
                     p.value = .mc_screen_p(md, family, fit1, co),
                     n = co$n[1], stringsAsFactors = FALSE)
        })
        rows <- Filter(Negate(is.null), rows)
        if (length(rows)) {
          screening_tab <- do.call(rbind, rows)
          rownames(screening_tab) <- NULL
          screening_tab$p.adjusted <- as.numeric(stats::p.adjust(
            screening_tab$p.value, method = p_adjust))
          screening_tab$significant <- !is.na(screening_tab$p.adjusted) &
            screening_tab$p.adjusted < alpha
          screening_tab <- screening_tab[order(screening_tab$p.value), ,
                                         drop = FALSE]
          rownames(screening_tab) <- NULL
        }
      }

      ## collinearity pruning (numeric predictors only)
      num_pred <- predictors[vapply(md[predictors], is.numeric, logical(1))]
      dropped_cc <- character(0)
      if (length(num_pred) >= 2L) {
        cm <- tryCatch(stats::cor(md[num_pred], use = "pairwise.complete.obs"),
                       error = function(e) NULL)
        if (!is.null(cm)) {
          rank_p <- if (!is.null(screening_tab))
            screening_tab$variable else num_pred
          ord <- c(intersect(rank_p, num_pred),
                   setdiff(num_pred, rank_p))
          keep <- character(0)
          for (v in ord) {
            if (length(keep) && any(abs(cm[v, keep, drop = FALSE]) >=
                                    collinear_cut, na.rm = TRUE)) {
              dropped_cc <- c(dropped_cc, v)
            } else keep <- c(keep, v)
          }
          if (length(dropped_cc)) {
            predictors <- setdiff(predictors, dropped_cc)
            add_dec("Collinearity", paste0("dropped ",
                                           paste(dropped_cc, collapse = ", ")),
                    sprintf("|r| >= %.2f with a higher-priority predictor",
                            collinear_cut))
            add_warn(sprintf(
              "%d numeric predictor(s) dropped for |r| >= %.2f; confounder selection must be confirmed clinically.",
              length(dropped_cc), collinear_cut))
          }
        }
      }

      ## EPV / sample-size cap
      k_cap <- max(1L, floor(n_model / epv_min))
      if (!is.null(max_predictors))
        k_cap <- min(k_cap, as.integer(max_predictors))
      ord <- if (!is.null(screening_tab)) screening_tab$variable else predictors
      predictors <- c(intersect(ord, predictors), setdiff(predictors, ord))
      if (length(predictors) > k_cap) {
        dropped_epv <- predictors[(k_cap + 1L):length(predictors)]
        predictors <- predictors[seq_len(k_cap)]
        add_dec("Model size", paste0("kept ", k_cap, " of ",
                                     k_cap + length(dropped_epv),
                                     " predictors"),
                sprintf("%s = %d, EPV limit %.0f; dropped %s",
                        if (family == "linear") "n" else "events",
                        n_model, epv_min,
                        paste(utils::head(dropped_epv, 8L), collapse = ", ")))
        add_warn(sprintf(
          "Model size capped at %d predictors (EPV limit %.0f). Pre-selection by univariable p value is exploratory.",
          k_cap, epv_min))
      }

      ## multivariable fit with graceful rank-deficiency fallback
      fit <- NULL
      fml <- NULL
      repeat {
        if (!length(predictors)) break
        fml <- stats::as.formula(paste(lhs, "~",
                                       paste(.mc_bt(predictors), collapse = " + ")))
        fit <- suppressWarnings(tryCatch(.mc_fit(fml, md, family),
                                         error = function(e) NULL))
        if (!is.null(fit)) {
          cf <- tryCatch(as.numeric(stats::coef(fit)), error = function(e) NULL)
          if (!is.null(cf) && length(cf) > 1L && all(is.finite(cf))) break
          fit <- NULL
        }
        gone <- predictors[length(predictors)]
        predictors <- predictors[-length(predictors)]
        add_dec("Model fitting", paste0("dropped ", gone),
                "fit failed or design matrix rank-deficient")
      }

      if (is.null(fit)) {
        add_warn("A multivariable model could not be fitted with the remaining predictors.")
        family <- "none"
      } else {
        audit$in_model[audit$variable %in% predictors] <- TRUE
        co <- .mc_coef(fit, family, conf, p_adjust, alpha)
        perf <- NULL
        checks <- NULL
        diagnostics <- list()
        # nobs.coxph() returns the number of events, not the number of rows.
        n_used <- if (family == "cox") as.integer(fit$n) else
          tryCatch(as.integer(stats::nobs(fit)),
                   error = function(e) nrow(md))
        n_coef <- length(stats::coef(fit))

        if (family == "linear") {
          s <- summary(fit)
          res <- as.numeric(stats::residuals(fit))
          fv <- as.numeric(stats::fitted(fit))
          sh_p <- if (length(res) >= 3L && length(res) <= 5000L)
            tryCatch(as.numeric(stats::shapiro.test(res)$p.value),
                     error = function(e) NA_real_) else NA_real_
          bp_p <- tryCatch({
            aux <- stats::lm(e2 ~ fv,
                             data = data.frame(e2 = res^2, fv = fv))
            stats::pchisq(n_used * summary(aux)$r.squared, df = 1,
                          lower.tail = FALSE)
          }, error = function(e) NA_real_)
          vif <- .mc_vif(fit)
          cd <- tryCatch(as.numeric(stats::cooks.distance(fit)),
                         error = function(e) rep(NA_real_, n_used))
          nflag <- sum(cd > 4 / length(cd), na.rm = TRUE)
          f_p <- tryCatch(as.numeric(stats::pf(
            s$fstatistic[1], s$fstatistic[2], s$fstatistic[3],
            lower.tail = FALSE)), error = function(e) NA_real_)
          perf <- data.frame(
            metric = c("n", "n_coefficients", "R2", "adjusted_R2",
                       "residual_sd", "AIC", "BIC", "F_p", "max_VIF",
                       "influential_rows"),
            value = c(n_used, n_coef, s$r.squared, s$adj.r.squared,
                      summary(fit)$sigma, stats::AIC(fit), stats::BIC(fit),
                      f_p, vif, nflag), stringsAsFactors = FALSE)
          checks <- rbind(
            .mc_check("Residual normality", sh_p, "Shapiro p > 0.05",
                      .mc_check_status(if (is.na(sh_p)) NA else sh_p > alpha),
                      "Shapiro-Wilk on residuals; over-sensitive when n is large"),
            .mc_check("Homoscedasticity", bp_p, "Breusch-Pagan p > 0.05",
                      .mc_check_status(if (is.na(bp_p)) NA else bp_p > alpha)),
            .mc_check("Multicollinearity", vif, "max VIF < 5",
                      .mc_check_status(if (is.na(vif)) NA else vif < 5)),
            .mc_check("Influential observations", nflag, "flagged = 0",
                      .mc_check_status(nflag == 0L), "Cook's D > 4/n"),
            .mc_check("Global model fit", f_p, "F test p < 0.05",
                      .mc_check_status(if (is.na(f_p)) NA else f_p < alpha)))
          diagnostics <- list(vif = vif, cooks_distance = cd,
                              residual_normality_p = sh_p,
                              breusch_pagan_p = bp_p)

        } else if (family == "logistic") {
          lr <- suppressWarnings(tryCatch(
            mc_logistic_report(md, fml, conf = conf,
                               auc_boot = auc_boot, seed = seed),
            error = function(e) NULL))
          if (!is.null(lr)) {
            perf <- lr$performance
            ck <- lr$checks[, c("check", "value", "threshold", "status"),
                            drop = FALSE]
            ck$note <- ""
            checks <- rbind(
              ck[, c("check", "value", "threshold", "status", "note")],
              .mc_check("Global model fit",
                        as.numeric(perf$value[perf$metric == "LikelihoodRatio_p"]),
                        "LRT p < 0.05",
                        .mc_check_status(
                          as.numeric(perf$value[perf$metric ==
                                                  "LikelihoodRatio_p"]) < alpha)))
            diagnostics <- list(calibration = lr$calibration, roc = lr$roc,
                                vif = lr$vif, influential = lr$influential,
                                conclusion = lr$conclusion, report = lr)
          } else {
            auc <- tryCatch({
              yr <- stats::model.response(stats::model.frame(fit))
              yv <- if (is.factor(yr) || is.character(yr))
                as.integer(factor(yr)) - 1L else as.numeric(yr)
              .mc_auc(as.numeric(stats::fitted(fit)), yv)
            }, error = function(e) NA_real_)
            perf <- data.frame(metric = c("n", "n_coefficients", "AUC",
                                          "AIC", "BIC", "max_VIF"),
                               value = c(n_used, n_coef, auc, stats::AIC(fit),
                                         stats::BIC(fit), .mc_vif(fit)),
                               stringsAsFactors = FALSE)
            checks <- .mc_check("Discrimination", auc, "AUC >= 0.70",
                                .mc_check_status(if (is.na(auc)) NA else
                                  auc >= 0.7))
          }

        } else {
          s <- summary(fit)
          conc <- tryCatch(as.numeric(s$concordance[1]),
                           error = function(e) NA_real_)
          lrt <- tryCatch(as.numeric(s$logtest[3]),
                          error = function(e) NA_real_)
          ph <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
          ph_p <- if (!is.null(ph) && "GLOBAL" %in% rownames(ph$table))
            as.numeric(ph$table["GLOBAL", "p"]) else NA_real_
          ph_tab <- if (!is.null(ph)) {
            tb <- as.data.frame(ph$table)
            tb$variable <- rownames(tb)
            tb[, c("variable", "chisq", "df", "p")]
          } else NULL
          ph_note <- "Schoenfeld residuals; per-covariate p values in diagnostics"
          if (!is.null(ph_tab)) {
            q <- ph_tab[ph_tab$variable != "GLOBAL" & !is.na(ph_tab$p), ,
                        drop = FALSE]
            if (nrow(q)) ph_note <- sprintf(
              "Schoenfeld residuals; smallest covariate PH p = %s (%s)",
              .mc_pfmt(q$p[which.min(q$p)]), q$variable[which.min(q$p)])
          }
          epv <- events / max(1L, n_coef)
          perf <- data.frame(
            metric = c("n", "events", "n_coefficients", "EPV",
                       "concordance", "LRT_p", "PH_global_p", "AIC", "BIC"),
            value = c(n_used, events, n_coef, epv, conc, lrt, ph_p,
                      stats::AIC(fit), stats::BIC(fit)),
            stringsAsFactors = FALSE)
          checks <- rbind(
            .mc_check("Proportional hazards", ph_p, "global PH p > 0.05",
                      .mc_check_status(if (is.na(ph_p)) NA else ph_p > alpha),
                      ph_note),
            .mc_check("Discrimination", conc, "concordance >= 0.70",
                      .mc_check_status(if (is.na(conc)) NA else conc >= 0.7)),
            .mc_check("Overfitting (EPV)", epv, "EPV >= 10",
                      .mc_check_status(epv >= epv_min)),
            .mc_check("Global model fit", lrt, "LRT p < 0.05",
                      .mc_check_status(if (is.na(lrt)) NA else lrt < alpha)))
          diagnostics <- list(ph_table = ph_tab, concordance = conc,
                              logrank_p = lrt, epv = epv)
        }

        model <- list(family = family, formula = paste(deparse(fml),
                                                       collapse = " "),
                      predictors = predictors, n = n_used, events = events,
                      coefficients = co, performance = perf,
                      checks = checks, diagnostics = diagnostics, fit = fit)
        add_dec("Model", sprintf("%s model fitted with %d predictor(s), n = %d",
                                 family, length(predictors), n_used),
                paste(deparse(fml), collapse = " "))
      }
    }
  }

  ## ---- written findings ---------------------------------------------------
  findings <- character(0)
  findings <- c(findings, sprintf(
    "Dataset: %d rows x %d columns; %d complete rows, %d duplicated rows, %d missing cells (%.1f%%).",
    n0, ncol(data), sum(stats::complete.cases(data)),
    sum(duplicated(data)), sum(is.na(data)),
    100 * sum(is.na(data)) / max(1L, n0 * ncol(data))))
  findings <- c(findings, sprintf(
    "Column audit: %d predictor(s) eligible, %d excluded (%s).",
    sum(audit$role == "predictor"), length(excluded),
    if (length(excluded)) paste(utils::head(excluded, 4L), collapse = ", ") else
      "none"))
  nd <- desc$type == "numeric"
  if (any(nd)) findings <- c(findings, sprintf(
    "Distribution: %d of %d numeric columns look approximately normal; correlation used %s.",
    sum(desc$distribution[nd] == "approximately normal"), sum(nd), cor_method))
  if (!is.null(comparison)) {
    ns <- sum(comparison$significant, na.rm = TRUE)
    top <- utils::head(comparison$variable[order(comparison$p.adjusted)], 3L)
    findings <- c(findings, sprintf(
      "Group comparison (%s): %d of %d tests significant after %s adjustment%s.",
      if (is.null(group)) "groups" else group, ns, nrow(comparison),
      p_adjust,
      if (ns) paste0("; smallest adjusted p: ", paste(top, collapse = ", ")) else
        ""))
  }
  if (!is.null(strong) && nrow(strong)) {
    strong <- strong[order(-abs(strong$r)), , drop = FALSE]
    findings <- c(findings, sprintf(
      "Correlation: %d pair(s) with |r| >= %.2f; strongest %s ~ %s (r = %.2f, adjusted p = %s).",
      nrow(strong), collinear_cut, strong$var1[1], strong$var2[1],
      strong$r[1], .mc_pfmt(strong$p.adjusted[1])))
  }
  if (!is.null(screening_tab)) {
    ns <- sum(screening_tab$significant, na.rm = TRUE)
    findings <- c(findings, sprintf(
      "Univariable screening: %d of %d predictors significant after %s adjustment.",
      ns, nrow(screening_tab), p_adjust))
  }
  if (!is.null(model)) {
    co <- model$coefficients
    sig <- if (!is.null(co)) co$term[co$significant] else character(0)
    findings <- c(findings, sprintf(
      "Multivariable %s model: %d coefficients, %d with p < %.2f%s.",
      model$family, if (is.null(co)) 0L else nrow(co), length(sig), alpha,
      if (length(sig)) paste0(" (", paste(utils::head(sig, 5L),
                                          collapse = ", "), ")") else ""))
    rev <- model$checks$check[model$checks$status == "REVIEW"]
    findings <- c(findings, sprintf(
      "Model checks: %d PASS, %d REVIEW%s.",
      sum(model$checks$status == "PASS"), length(rev),
      if (length(rev)) paste0(" (", paste(rev, collapse = ", "), ")") else ""))
  } else {
    findings <- c(findings,
                  "No regression model was fitted; descriptives and comparisons are still available.")
  }

  ## ---- overview -----------------------------------------------------------
  overview <- data.frame(
    n_rows = n0, n_cols = ncol(data),
    n_complete_rows = sum(stats::complete.cases(data)),
    n_duplicated_rows = sum(duplicated(data)),
    n_numeric = sum(audit$type == "numeric"),
    n_categorical = sum(audit$type %in% c("categorical", "binary")),
    n_missing_cells = sum(is.na(data)),
    pct_missing_cells = round(100 * sum(is.na(data)) /
                                max(1L, n0 * ncol(data)), 1),
    n_excluded = length(excluded),
    n_predictors = sum(audit$role == "predictor"),
    n_in_model = if (is.null(model)) 0L else length(model$predictors),
    analysis = if (is.null(model)) family else model$family,
    outcome = if (is.null(outcome)) "" else outcome,
    group = if (is.null(group)) "" else group,
    p_adjust = p_adjust, alpha = alpha, conf = conf,
    generated = as.character(Sys.time()), stringsAsFactors = FALSE)

  decisions <- if (length(dec$rows)) {
    d <- do.call(rbind, dec$rows)
    d$step_no <- seq_len(nrow(d))
    d[, c("step_no", "step", "decision", "reason")]
  } else {
    data.frame(step_no = integer(0), step = character(0),
               decision = character(0), reason = character(0),
               stringsAsFactors = FALSE)
  }
  if (!is.null(model))
    warns <- c(warns,
               "Automated variable pre-selection is exploratory; confirm confounder selection against the clinical analysis plan.")
  if (n0 < 30L)
    warns <- c(warns, "Small sample: asymptotic confidence intervals and p values are unreliable.")

  out <- list(overview = overview, columns = audit, descriptives = desc,
              baseline = baseline, comparison = comparison,
              correlation = correlation, screening = screening_tab,
              model = model, findings = findings, decisions = decisions,
              warnings = warns, data = data, call = match.call())
  class(out) <- "mc_auto_analysis"
  if (!is.null(report)) mc_auto_report(out, file = report, quiet = TRUE)
  out
}

#' Print a compact summary of an automated analysis
#'
#' @param x An object of class `mc_auto_analysis`.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.mc_auto_analysis <- function(x, ...) {
  ov <- x$overview
  cat("macodgent automated analysis\n")
  cat(paste(rep("=", 64), collapse = ""), "\n")
  cat(sprintf("Data: %s rows x %s columns | analysis: %s | outcome: %s\n",
              ov$n_rows, ov$n_cols, ov$analysis,
              if (nzchar(ov$outcome)) ov$outcome else "(none)"))
  cat(sprintf("Audit: %s predictors, %s excluded, %s in model | p adjustment: %s\n",
              ov$n_predictors, ov$n_excluded, ov$n_in_model, ov$p_adjust))
  if (!is.null(x$comparison)) cat(sprintf(
    "Group comparison: %d tests, %d significant after adjustment\n",
    nrow(x$comparison), sum(x$comparison$significant, na.rm = TRUE)))
  if (!is.null(x$correlation)) cat(sprintf(
    "Correlation: %d pairs, %d with |r| >= 0.8\n", nrow(x$correlation),
    sum(abs(x$correlation$r) >= 0.8, na.rm = TRUE)))
  if (!is.null(x$model)) {
    cat("\nModel checks\n")
    print(x$model$checks[, c("check", "value", "status")],
          row.names = FALSE, digits = 3)
  }
  cat("\nKey findings\n")
  cat(paste0("- ", x$findings, collapse = "\n"), "\n")
  if (length(x$warnings)) {
    cat("\nCaveats\n")
    cat(paste0("- ", x$warnings, collapse = "\n"), "\n")
  }
  cat("\nComponents: ",
      paste(setdiff(names(x), c("call", "data")), collapse = ", "), "\n",
      sep = "")
  invisible(x)
}
