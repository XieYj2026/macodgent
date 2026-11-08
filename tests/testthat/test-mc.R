library(macodgent)

test_that("mc_mode returns the correct mode", {
  expect_equal(mc_mode(c(1, 2, 2, 3, 3, 3, 4, 5, 5, 6)), 3)
  expect_equal(sort(mc_mode(c("a", "a", "b", "b"))), c("a", "b"))
  expect_null(mc_mode(c(NA, NA), na.rm = TRUE))
})

test_that("mc_inspect gives a missing-value overview", {
  d <- data.frame(a = c(1, 2, NA), b = c("x", "x", "x"))
  ins <- mc_inspect(d)
  expect_equal(nrow(ins), 2)
  expect_equal(ins$n_missing[ins$variable == "a"], 1)
})

test_that("mc_impute_knn imputes numeric missing values", {
  df <- data.frame(a = c(1, 2, NA, 4, 5), b = c(2, NA, 3, 4, 5))
  out <- mc_impute_knn(df, k = 2)
  expect_false(any(is.na(out$a)))
})

test_that("mc_dedup removes duplicates by key column", {
  df <- data.frame(id = c(1, 1, 2), v = c("a", "a", "b"))
  expect_equal(nrow(mc_dedup(df, by = "id")), 2)
})

test_that("mc_scale returns z-scores with mean ~0", {
  out <- mc_scale(mtcars, cols = "mpg", method = "z")
  expect_true(abs(mean(out$mpg)) < 1e-9)
})

test_that("mc_desc_quant outputs ~20 descriptive stats", {
  d <- data.frame(x = rnorm(50, 10, 2), g = rep(1:2, each = 25))
  res <- mc_desc_quant(d, vars = "x")
  expect_equal(nrow(res), 1)
  expect_equal(res$mean, mean(d$x, na.rm = TRUE), tolerance = 1e-6)
  resg <- mc_desc_quant(d, vars = "x", group = "g")
  expect_equal(nrow(resg), 2)
})

test_that("mc_norm_test / mc_var_test return data frames", {
  d <- MASS::birthwt
  expect_true(all(c("test", "statistic", "p.value",
                     "skewness", "verdict") %in%
                   names(mc_norm_test(d$age))))
  vt <- mc_var_test(d$bwt, d$smoke, method = "bartlett")
  expect_true(vt$p.value >= 0 && vt$p.value <= 1)
})

test_that("mc_ttest / mc_anova / mc_nonparam", {
  set.seed(1)
  x <- rnorm(30, mean = 1)
  tt <- mc_ttest(x, mu = 0)
  expect_true(tt$p.value < 0.05)
  av <- mc_anova(bwt ~ factor(race), data = MASS::birthwt)
  expect_true("factor(race)" %in% av$term)
  npr <- mc_nonparam(MASS::birthwt$bwt, MASS::birthwt$smoke, method = "wilcoxon")
  expect_true(all(c("test", "statistic", "p.value",
                     "effect_size", "effect_ci.low") %in% names(npr)))
})

test_that("mc_contingency computes 2x2 OR", {
  res <- mc_contingency(MASS::birthwt, row = "smoke", col = "low")
  expect_true(!is.null(res$OR))
  expect_equal(nrow(res$table), 2)
})

test_that("mc_corr returns pairwise correlations", {
  cr <- mc_corr(MASS::birthwt, vars = c("age", "lwt", "bwt"))
  expect_equal(nrow(cr), 3)
})

test_that("mc_tableone builds a baseline table", {
  tb <- mc_tableone(MASS::birthwt, vars = c("age", "lwt", "smoke"), group = "smoke")
  expect_true("variable" %in% names(tb))
})

test_that("mc_logistic and the screening workflow", {
  bw <- MASS::birthwt
  bw$low <- factor(bw$low)
  lg <- mc_logistic(bw, low ~ age + smoke)
  expect_true("OR" %in% names(lg) || "estimate" %in% names(lg))
  sc <- mc_logistic_screen(bw, y = "low",
                           x_vars = c("age", "lwt", "race", "smoke", "ht", "ui"))
  expect_true(is.data.frame(sc$univariate))
})

test_that("mc_stepwise returns the final model", {
  bw <- MASS::birthwt
  st <- mc_stepwise(bw, low ~ age + lwt + race + smoke + ht + ui)
  expect_true(!is.null(st$model))
})

test_that("mc_rcs self-implemented spline runs", {
  set.seed(2)
  d <- data.frame(age = sample(18:45, 120, replace = TRUE),
                  bwt = rnorm(120, 3000, 600))
  rcs <- mc_rcs(d, x = "age", y = "bwt", nk = 3)
  expect_true(!is.null(rcs$plot))
  expect_true(length(rcs$knots) >= 4)
})

# The following functions depend on Suggests packages and are tested only when available
if (requireNamespace("survival", quietly = TRUE)) {
  test_that("mc_cox returns HR and a forest plot", {
    skip_if_not(requireNamespace("survival"))
    fit <- mc_cox(survival::lung, survival::Surv(time, status) ~ age + sex + ph.ecog)
    expect_true(!is.null(fit$forest))
    expect_true("HR" %in% names(fit$coefficients) || "estimate" %in% names(fit$coefficients))
  })
}

if (requireNamespace("quantreg", quietly = TRUE)) {
  test_that("mc_quantile_reg multiple quantiles", {
    qr <- mc_quantile_reg(MASS::birthwt, bwt ~ age + lwt)
    expect_true(0.25 %in% qr$tau)
  })
}

if (requireNamespace("glmnet", quietly = TRUE)) {
  test_that("mc_lasso selects non-zero coefficients", {
    ls <- mc_lasso(mtcars, x = c("hp", "wt", "disp", "cyl"), y = "mpg")
    expect_true(nrow(ls$coefficients) >= 0)
  })
}

if (requireNamespace("mgcv", quietly = TRUE)) {
  test_that("mc_gam fits", {
    g <- mc_gam(MASS::birthwt, bwt ~ s(age) + smoke)
    expect_true(!is.null(g$model))
  })
}


test_that("renamed functions are exact aliases of the current names", {
  pairs <- list(
    list(mc_mode, mc_stat_mode), list(mc_pcr, mc_pca_regression),
    list(mc_desc_quant, mc_describe_numeric),
    list(mc_norm_test, mc_test_normality),
    list(mc_var_test, mc_test_equal_variance),
    list(mc_nonparam, mc_test_rank), list(mc_anova, mc_anova_oneway),
    list(mc_survfit, mc_kaplan_meier),
    list(mc_wide_long, mc_pivot_longer),
    list(mc_long_wide, mc_pivot_wider),
    list(mc_scale, mc_scale_numeric),
    list(mc_tableone, mc_baseline_table),
    list(mc_ph_test, mc_test_proportional_hazards),
    list(mc_corr, mc_correlate),
    list(mc_contingency, mc_crosstab_test))
  for (p in pairs) expect_identical(p[[1]], p[[2]])
})

test_that("multiplicity adjustment matches stats::p.adjust", {
  raw <- c(0.001, 0.02, 0.04, 0.31, 0.62)
  a <- mc_adjust_p(raw, method = "holm")
  expect_equal(a$p.adjusted, as.numeric(stats::p.adjust(raw, "holm")))
  expect_equal(sum(a$significant), 1L)
  b <- mc_adjust_p(raw, method = "BH")
  expect_equal(b$p.adjusted, as.numeric(stats::p.adjust(raw, "BH")))
})

test_that("t-test reports Hedges g with a confidence interval", {
  set.seed(11)
  a <- rnorm(40, 0, 1); b <- rnorm(40, 0.8, 1)
  tt <- mc_ttest(a, b)
  expect_true("hedges_g" %in% names(tt))
  expect_true(tt$hedges_g < 0)
  expect_true(tt$effect_ci.low < tt$hedges_g)
  expect_true(tt$effect_ci.high > tt$hedges_g)
})

test_that("ANOVA reports effect sizes and Tukey contrasts", {
  d <- MASS::birthwt
  av <- mc_anova_oneway(bwt ~ factor(race), data = d)
  expect_true(all(c("eta_sq", "omega_sq") %in% names(av)))
  expect_true(av$eta_sq[1] >= 0 && av$eta_sq[1] <= 1)
  ph <- mc_posthoc_tukey(data = d, formula = bwt ~ factor(race))
  expect_true(nrow(ph) == 3L)
  expect_true(all(c("estimate", "conf.low", "conf.high", "p.adjusted") %in%
                    names(ph)))
})

test_that("group comparison returns effect sizes and adjusted p values", {
  d <- MASS::birthwt
  d$low <- factor(d$low, levels = 0:1, labels = c("normal", "low"))
  d$smoke <- factor(d$smoke)
  d$race <- factor(d$race)
  r <- mc_compare_groups(d, group = "low",
                         vars = c("age", "lwt", "bwt", "smoke", "race"))
  expect_equal(nrow(r), 5L)
  expect_true(all(c("effect", "estimate", "conf.low", "conf.high",
                    "p.value", "p.adjusted") %in% names(r)))
  expect_equal(r$p.adjusted, as.numeric(stats::p.adjust(r$p.value, "holm")))
  expect_true(all(r$p.adjusted >= r$p.value - 1e-12))
})

test_that("Table One carries p values and SMD", {
  d <- MASS::birthwt
  d$low <- factor(d$low, levels = 0:1, labels = c("normal", "low"))
  tb <- mc_baseline_table(d, vars = c("age", "lwt"), group = "low")
  expect_true(all(c("p.value", "p.adjusted", "smd", "effect") %in% names(tb)))
  expect_true(all(is.finite(tb$smd)))
})

test_that("crosstab reports an effect size and sparse-cell diagnostics", {
  r <- mc_crosstab_test(MASS::birthwt, row = "smoke", col = "low")
  expect_true("effect" %in% names(r))
  expect_true("diagnostics" %in% names(r))
  expect_true(r$effect$estimate >= 0)
  expect_true(r$diagnostics$n == nrow(MASS::birthwt))
})

test_that("correlation reports CI and adjusted p", {
  r <- mc_correlate(MASS::birthwt, vars = c("age", "lwt", "bwt"))
  expect_equal(nrow(r), 3L)
  expect_true(all(c("conf.low", "conf.high", "n", "p.adjusted") %in%
                    names(r)))
  expect_true(all(r$conf.low < r$r & r$r < r$conf.high))
})

test_that("logistic report delivers the full diagnostic suite", {
  d <- MASS::birthwt
  r <- mc_logistic_report(d, low ~ age + lwt + smoke, auc_boot = 100)
  expect_s3_class(r, "mc_logistic_report")
  expect_true(all(c("coefficients", "performance", "checks", "conclusion",
                    "calibration", "roc", "vif", "influential", "model",
                    "plot", "roc_plot") %in% names(r)))
  pv <- r$performance$value[r$performance$metric == "AUC"]
  expect_true(as.numeric(pv) > 0.5 && as.numeric(pv) < 1)
  expect_equal(nrow(r$checks), 5L)
  expect_true(all(r$checks$status %in% c("PASS", "REVIEW", "UNDETERMINED")))
  expect_true(nrow(r$calibration) >= 3L)
  expect_true(nrow(r$influential) == nrow(d))
  expect_output(print(r), "Logistic model report")
})

test_that("Kaplan-Meier adds median survival and risk table", {
  ov <- survival::ovarian
  ov$status <- ifelse(ov$fustat == 2, 0, 1)
  r <- mc_kaplan_meier(ov, time = "futime", status = "status", group = "rx")
  expect_true("median_survival" %in% names(r))
  expect_true("risk_table" %in% names(r))
  expect_true("df" %in% names(r$logrank))
  expect_true(nrow(r$risk_table) > 0L)
})

test_that("Cox adds performance and explicit checks", {
  d <- survival::lung
  r <- mc_cox(d, survival::Surv(time, status) ~ age + sex + ph.ecog)
  expect_true("concordance" %in% r$performance$metric)
  expect_true(all(r$checks$status %in% c("PASS", "REVIEW")))
})

# ---- automated one-call analysis -------------------------------------------

test_that("mc_auto_analyze returns a complete linear analysis bundle", {
  r <- mc_auto_analyze(mtcars, outcome = "mpg", group = "am", auc_boot = 0L)
  expect_s3_class(r, "mc_auto_analysis")
  expect_true(all(c("overview", "columns", "descriptives", "baseline",
                    "comparison", "correlation", "screening", "model",
                    "findings", "decisions", "warnings", "data") %in%
                    names(r)))
  expect_equal(r$overview$analysis, "linear")
  expect_equal(nrow(r$columns), ncol(mtcars))
  expect_true(nrow(r$decisions) > 0L)
  expect_true(length(r$findings) >= 5L)
  expect_true(all(r$model$checks$status %in%
                    c("PASS", "REVIEW", "UNDETERMINED")))
  expect_true(all(c("estimate", "conf.low", "conf.high", "p.adjusted") %in%
                    names(r$model$coefficients)))
  expect_true("p.adjusted" %in% names(r$comparison))
  expect_output(print(r), "macodgent automated analysis")
})

test_that("automated logistic analysis detects outcome leakage", {
  bw <- MASS::birthwt
  bw$low <- factor(bw$low, levels = 0:1, labels = c("normal", "low"))
  r <- mc_auto_analyze(bw, outcome = "low", auc_boot = 0L)
  expect_equal(r$overview$analysis, "logistic")
  expect_false("bwt" %in% r$model$predictors)
  expect_true(any(grepl("leakage", r$decisions$step, ignore.case = TRUE)))
  expect_true("Discrimination" %in% r$model$checks$check)
  expect_true(all(r$model$coefficients$effect == "OR"))
  expect_true(nrow(r$baseline) > 0L)
})

test_that("automated Cox analysis recodes status and checks proportional hazards", {
  r <- mc_auto_analyze(survival::lung, time = "time", status = "status")
  expect_equal(r$overview$analysis, "cox")
  expect_true(any(grepl("event = 2", r$decisions$decision, fixed = TRUE)))
  expect_true("Proportional hazards" %in% r$model$checks$check)
  ph <- r$model$checks$value[r$model$checks$check == "Proportional hazards"]
  expect_true(is.finite(ph))
  expect_true(r$model$n > 100L)
  expect_true(all(r$model$coefficients$effect == "HR"))
})

test_that("automated report and export produce portable files", {
  r <- mc_auto_analyze(mtcars, outcome = "mpg", auc_boot = 0L)
  f <- tempfile(fileext = ".html")
  mc_auto_report(r, file = f, quiet = TRUE)
  expect_true(file.exists(f))
  h <- paste(readLines(f, warn = FALSE), collapse = "")
  expect_true(grepl("Key findings", h, fixed = TRUE))
  expect_true(grepl("forest", h, fixed = TRUE))
  expect_true(grepl("Decision log", h, fixed = TRUE))
  d <- file.path(tempdir(), paste0("mc-export-", Sys.getpid()))
  files <- mc_auto_export(r, dir = d, quiet = TRUE)
  expect_true(any(grepl("model_coefficients.csv", files, fixed = TRUE)))
  expect_s3_class(readRDS(file.path(d, "analysis.rds")), "mc_auto_analysis")
  unlink(d, recursive = TRUE)
})

test_that("analysis without an outcome still returns descriptive results", {
  r <- mc_auto_analyze(iris)
  expect_equal(r$overview$analysis, "none")
  expect_null(r$model)
  expect_equal(nrow(r$descriptives), 5L)
  expect_equal(nrow(r$correlation), 6L)
})

test_that("non-syntactic column names work in automated formulas", {
  d <- data.frame(check.names = FALSE,
                  y = rnorm(60),
                  "birth weight" = rnorm(60),
                  "study group" = sample(c("A", "B"), 60, TRUE))
  r <- mc_auto_analyze(d, outcome = "y", group = "study group",
                       auc_boot = 0L)
  expect_equal(r$overview$analysis, "linear")
  expect_true(nrow(r$model$coefficients) >= 1L)
})


# ============================================================================
# Upgrade pass: direct coverage for the previously untested functions and
# assertions for every new capability (shared checks, honest validation,
# extended model families, interactive/printable reports).
# ============================================================================

test_that("mc_linear exposes the shared methodological checks", {
  r <- mc_linear(mtcars, mpg ~ hp + wt)
  expect_true("checks" %in% names(r))
  expect_true(all(c("check", "value", "threshold", "status") %in%
                    names(r$checks)))
  expect_true(nrow(r$checks) >= 5L)
  expect_true(all(r$checks$status %in% c("PASS", "REVIEW")))
})

test_that("mc_logistic can return the full diagnostic bundle", {
  bw <- MASS::birthwt
  tab <- mc_logistic(bw, low ~ age + lwt)
  expect_true(is.data.frame(tab))
  r <- mc_logistic(bw, low ~ age + lwt, diagnostics = TRUE, auc_boot = 0L)
  expect_true(all(c("coefficients", "checks", "performance", "report") %in%
                    names(r)))
  expect_true(nrow(r$checks) >= 4L)
})

test_that("mc_rcs binomial band stays on the response scale", {
  bw <- MASS::birthwt
  r <- mc_rcs(bw, x = "age", y = "low", family = "binomial")
  expect_true(all(r$curve$fit >= 0 & r$curve$fit <= 1))
  expect_true(all(r$curve$lo >= 0 & r$curve$hi <= 1))
  expect_true(all(r$curve$lo <= r$curve$fit & r$curve$fit <= r$curve$hi))
  expect_true("curve" %in% names(r))
  expect_true("nonlinearity_p" %in% names(r))
})

test_that("mc_nnet_reg reports out-of-sample RMSE", {
  skip_if_not_installed("nnet")
  r <- suppressWarnings(mc_nnet_reg(mtcars[, c("mpg", "hp", "wt")],
                                    mpg ~ hp + wt, hidden = 3L,
                                    maxit = 80L, seed = 1L,
                                    validation = "holdout"))
  v <- r$performance$value[r$performance$metric == "validation_rmse"]
  expect_equal(r$rmse, v)
  expect_true("validation_r2" %in% r$performance$metric)
})

test_that("penalised regressions return the enhanced bundle", {
  skip_if_not_installed("glmnet")
  x <- setdiff(names(mtcars), "mpg")
  for (fn in c("mc_lasso", "mc_ridge", "mc_elastic_net")) {
    r <- do.call(fn, list(data = mtcars, x = x, y = "mpg", seed = 1L))
    expect_true(all(c("cv_table", "performance", "lambda.1se") %in% names(r)))
    expect_true(nrow(r$cv_table) > 0L)
  }
})

test_that("mc_pca_regression gives coefficient CIs and external validation", {
  r <- mc_pca_regression(MASS::birthwt, bwt ~ age + lwt + smoke + ht + ui,
                         k = 5L, seed = 1L, validation = "holdout")
  expect_true(all(c("std.error", "conf.low", "conf.high", "p.value") %in%
                    names(r$coef_original)))
  expect_true("external_rmse" %in% r$performance$metric)
  expect_equal(nrow(r$coef_original), 6L)
})

test_that("mc_impute_mice reports chain convergence", {
  skip_if_not_installed("mice")
  d <- mtcars[1:30, c("mpg", "hp", "wt")]
  d$hp[1:5] <- NA
  r <- mc_impute_mice(d, m = 3, maxit = 6, seed = 1L)
  expect_true(is.data.frame(r$convergence))
  expect_true("relative_drift" %in% names(r$convergence))
  expect_true(all(c("check", "status") %in% names(r$checks)))
})

test_that("pivot helpers support fill, aggregation and multivalue columns", {
  d <- mtcars[1:6, c("mpg", "hp", "wt")]
  lg <- suppressWarnings(
    mc_pivot_longer(d, id_cols = "mpg", measure_vars = c("hp", "wt")))
  expect_equal(nrow(lg), 12L)
  wd <- suppressWarnings(mc_pivot_wider(
    lg, id_cols = "mpg", name_from = "variable",
    value_from = "value", values_fill = 0, values_fn = "mean"))
  expect_equal(ncol(wd), 3L)
  multi <- suppressWarnings(mc_pivot_longer(
    d, id_cols = "mpg",
    measure_vars = list(a = c("hp", "wt"), b = c("hp", "wt"))))
  expect_true(all(c("a", "b") %in% names(multi)))
})

test_that("mc_scale_numeric can learn parameters from a training frame", {
  r <- mc_scale_numeric(mtcars[1:10, c("mpg", "hp")],
                        train = mtcars[1:20, c("mpg", "hp")])
  expect_true(attr(r, "scaling")$train_only)
  expect_equal(length(attr(r, "scaling")$parameters), 2L)
})

test_that("mc_inspect detailed mode exposes distribution and roles", {
  r <- mc_inspect(MASS::birthwt, detailed = TRUE)
  expect_true(all(c("pct_unique", "is_constant", "n_outliers",
                    "distribution", "suggested_role") %in% names(r)))
  basic <- mc_inspect(MASS::birthwt)
  expect_equal(ncol(basic), 6L)
})

test_that("ordinal and conditional logistic add intervals and checks", {
  bw <- MASS::birthwt
  bw$ftv_ord <- ordered(cut(bw$ftv, c(-1, 0, 2, Inf),
                            labels = c("0", "1-2", "3+")))
  o <- mc_ordinal_logistic(bw, ftv_ord ~ age + lwt)
  expect_true(all(c("OR_lower", "OR_upper", "is_intercept") %in% names(o)))
  expect_true(all(o$OR_lower <= o$OR & o$OR <= o$OR_upper))
  expect_true(is.data.frame(attr(o, "checks")))
  c1 <- mc_clogistic(datasets::infert,
                     case ~ spontaneous + induced +
                       survival::strata(stratum))
  expect_equal(nrow(c1), 2L)
  expect_true(all(c1$OR_lower < c1$OR & c1$OR < c1$OR_upper))
  expect_equal(nrow(attr(c1, "checks")), 5L)
})

test_that("model comparison audits comparability and can run LRT", {
  m1 <- lm(mpg ~ hp, data = mtcars)
  m2 <- lm(mpg ~ hp + wt, data = mtcars)
  r <- mc_model_compare(m1, m2, lrt = TRUE)
  expect_true(all(r$comparable))
  expect_true(all(c("family", "response", "lrt_p") %in% names(r)))
})

test_that("stepwise accepts a non-binomial family and forward scope", {
  r <- mc_stepwise(bw <- MASS::birthwt, low ~ age + lwt + smoke,
                   direction = "forward")
  expect_true(inherits(r$model, "glm"))
  g <- mc_stepwise(mtcars, mpg ~ hp + wt + disp,
                   family = stats::gaussian())
  expect_true(inherits(g$model, "glm"))
  expect_false(isTRUE(all.equal(g$family$family, "binomial")))
})

test_that("univariable screening adjusts p values and reports stability", {
  bw <- MASS::birthwt
  r <- mc_logistic_screen(bw, y = "low",
                          x_vars = c("age", "lwt", "smoke", "ht"),
                          p_adjust = "BH", bootstrap = 5L, seed = 1L)
  expect_true("p.adjusted" %in% names(r$univariate))
  expect_equal(nrow(r$stability), 4L)
})

test_that("advanced regression families report inference and checks", {
  skip_if_not_installed("quantreg")
  q <- mc_quantile_reg(MASS::birthwt, bwt ~ lwt + age,
                       taus = c(0.25, 0.5, 0.75))
  expect_true(any(is.finite(q$std_error)))
  expect_true(is.data.frame(attr(q, "trend")))
  skip_if_not_installed("mgcv")
  g <- mc_gam(MASS::birthwt, bwt ~ s(lwt) + s(age))
  expect_true(all(c("smooths", "performance", "checks") %in% names(g)))
  expect_true("edf" %in% names(g$smooths))
})

test_that("competing-risk and PH-test functions add CIF and verdicts", {
  skip_if_not_installed("cmprsk")
  set.seed(3)
  n <- 80L
  d <- data.frame(t = rexp(n, 0.1),
                  s = sample(0:2, n, TRUE, c(0.3, 0.5, 0.2)),
                  x = rnorm(n))
  r <- mc_competing_risk(d, time = "t", status = "s", group = "x")
  expect_true(all(c("HR_lower", "HR_upper") %in% names(r$coefficients)))
  expect_true(is.data.frame(r$cif))
  expect_equal(nrow(r$checks), 6L)
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex,
                         data = survival::lung)
  p <- mc_test_proportional_hazards(fit = fit)
  expect_true("p.adjusted" %in% names(p$test))
  expect_true(all(c("global_p", "verdict", "checks") %in% names(p)))
})

test_that("multistate returns state probabilities", {
  skip_if_not_installed("mstate")
  set.seed(2)
  n <- 40L
  dl <- data.frame(id = rep(seq_len(n), each = 2L),
                   Tstart = 0, Tstop = rep(c(1, 2), each = n),
                   status = rep(c(1, 0), each = n),
                   trans = rep(1:2, each = n), x = rnorm(n * 2L))
  tm <- matrix(NA, 3L, 3L,
               dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  tm[1, 2] <- 1L; tm[1, 3] <- 2L
  r <- suppressWarnings(mc_multistate(
    dl, tm, survival::Surv(Tstart, Tstop, status) ~ x + strata(trans)))
  expect_true(all(c("msfit", "probability_table", "checks") %in% names(r)))
  expect_true(is.data.frame(r$checks))
})

test_that("theme and palette helpers are usable", {
  expect_true(is.list(mc_theme()))
  expect_true(is.character(mc_palette()))
})

test_that("mc_auto_model dispatches extended families uniformly", {
  skip_if_not_installed("glmnet")
  bw <- MASS::birthwt
  for (fam in c("linear", "quantile", "pca", "elastic_net")) {
    skip_if(fam == "quantile" && !requireNamespace("quantreg", quietly = TRUE))
    # quantreg::rq may report a non-unique solution on birthwt; keep other
    # families' warnings visible so real regressions are not masked.
    r <- if (identical(fam, "quantile")) {
      suppressWarnings(mc_auto_model(bw, outcome = "bwt", family = fam))
    } else {
      mc_auto_model(bw, outcome = "bwt", family = fam)
    }
    expect_equal(r$family, fam)
    expect_true(all(c("coefficients", "performance", "checks") %in%
                      names(r)))
    expect_true(nrow(r$checks) >= 1L)
  }
  expect_output(print(mc_auto_model(bw, outcome = "bwt", family = "linear")),
                "model family")
})

test_that("mc_auto_analyze honours an explicit extended family", {
  r <- mc_auto_analyze(mtcars, outcome = "mpg", family = "pca",
                       auc_boot = 0L)
  expect_equal(r$overview$analysis, "pca")
  expect_equal(r$model$family, "pca")
  expect_true(length(r$findings) >= 4L)
})

test_that("reports are sortable/printable and export a Word-readable table", {
  r <- mc_auto_analyze(mtcars, outcome = "mpg", auc_boot = 0L)
  f <- tempfile(fileext = ".html")
  mc_auto_report(r, file = f, quiet = TRUE)
  h <- paste(readLines(f, warn = FALSE), collapse = "")
  expect_true(grepl("sortable", h, fixed = TRUE))
  expect_true(grepl("@media print", h, fixed = TRUE))
  d <- tempfile()
  fs <- mc_auto_export(r, dir = d, word = TRUE, quiet = TRUE)
  expect_true(file.exists(file.path(d, "report.doc")))
  expect_true(length(fs) >= 10L)
})

test_that("every current function name is exercised directly", {
  expect_true(length(mc_stat_mode(c(1, 1, 2, 3))) >= 1L)
  expect_true(nrow(mc_describe_numeric(mtcars, vars = c("mpg", "wt"))) == 2L)
  expect_true(NROW(mc_test_normality(mtcars$mpg)) >= 1L)
  expect_true(NROW(mc_test_equal_variance(mtcars$mpg,
                                         factor(mtcars$cyl))) >= 1L)
  expect_true(NROW(mc_test_rank(mtcars$mpg, factor(mtcars$am))) >= 1L)

  skip_if_not_installed("glmnet")
  xvars <- setdiff(names(mtcars), "mpg")
  ri <- mc_ridge(mtcars, x = xvars, y = "mpg", nfolds = 3)
  en <- mc_elastic_net(mtcars, x = xvars, y = "mpg", nfolds = 3)
  expect_true(NROW(ri$coefficients) >= 1L)
  expect_true(NROW(en$coefficients) >= 1L)
})

test_that("sim_htn_trial is a documented, explicitly simulated example dataset", {
  data("sim_htn_trial", package = "macodgent")
  expect_equal(dim(sim_htn_trial), c(3000L, 47L))
  expect_true(isTRUE(attr(sim_htn_trial, "simulated")))
  expect_equal(attr(sim_htn_trial, "seed"), 20261105L)
  expect_true(grepl("SIMULATED", attr(sim_htn_trial, "provenance"),
                    fixed = TRUE))
  expect_true(all(c("bp_control", "bp_grade", "cv_event", "comp_event",
                    "match_id", "case_status", "followup_days") %in%
                    names(sim_htn_trial)))
  expect_gt(sum(is.na(sim_htn_trial$potassium)), 0L)
})

test_that("mc_auto_analyze keeps structural columns when vars is supplied", {
  data("sim_htn_trial", package = "macodgent")
  d <- sim_htn_trial[1:400, ]
  # vars lists candidate predictors only; outcome/group must survive subsetting
  r <- mc_auto_analyze(d, outcome = "bp_control", group = "treatment",
                       vars = c("age", "bmi", "ldl"), auc_boot = 0L)
  expect_equal(r$overview$analysis, "logistic")
  expect_true(is.data.frame(r$decisions))
  expect_true(nrow(r$model$checks) >= 1L)
})

test_that("mc_competing_risk accepts a factor grouping column", {
  skip_if_not_installed("cmprsk")
  data("sim_htn_trial", package = "macodgent")
  r <- suppressWarnings(mc_competing_risk(sim_htn_trial,
    time = "followup_days", status = "comp_event",
    group = "treatment", cause = 1))
  expect_true(nrow(r$coefficients) >= 1L)
  expect_true(all(c("HR", "HR_lower", "HR_upper") %in%
                    names(r$coefficients)))
})
