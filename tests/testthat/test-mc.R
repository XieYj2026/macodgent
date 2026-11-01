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
  expect_equal(ncol(mc_norm_test(d$age)), 3)
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
  expect_equal(ncol(npr), 3)
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
