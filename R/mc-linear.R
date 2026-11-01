# =============================================================================
# 软件名称    : macodgent 医学数据科学研究工具包
# 著作权人    : Dr.Xie
# 版本号      : 1.0.0
# 开发语言    : R
# 开发完成日期: 2026-10-30
# 权利声明    : 本程序为著作权人独立开发，受《计算机软件保护条例》保护；
#              未经书面许可，不得复制、发行、出租或将其用于商业目的。
# -----------------------------------------------------------------------------
# 模块        : 线性回归（第 5 章 常见类型回归分析的补充）
# 本文件包含  : mc_linear（一元/多元线性回归，含自研诊断模块）
# =============================================================================

#' 线性回归（含自研诊断模块）/ Linear regression with self-implemented diagnostics
#'
#' 对应第 5 章「常见类型回归分析」（5.2 一元/多元线性回归、5.6 模型诊断绘图）。
#' 在基础 \code{stats::lm} 之上，额外提供三项**自研**诊断：① 基于重抽样的
#' 自助法（bootstrap）系数置信区间；② 基于辅助回归的方差膨胀因子（VIF）多重
#' 共线性诊断；③ 四张标准回归诊断图（残差 vs 拟合、正态 Q-Q、尺度-位置、
#' Cook 距离）。所有结果整理为整洁 \code{data.frame}，便于直接入报告。
#'
#' Corresponds to Chapter 5 ("Linear regression", 5.2 & 5.6). Wraps
#' \code{stats::lm} and adds three *self-implemented* diagnostics: a bootstrap
#' confidence interval for coefficients, variance-inflation-factor (VIF)
#' multicollinearity checks via auxiliary regressions, and the four standard
#' regression diagnostic plots.
#'
#' @details
#' 用途：拟合一元或多元线性回归，并自动产出可报告的系数表与诊断结果。
#' 方法：以 \code{stats::lm} 做最小二乘拟合；在此之上自研实现三项诊断——
#' （1）\strong{自助法 CI}：对观测行做有放回重抽样、反复重拟合，取系数经验
#' 分位数作为区间，不依赖任何正态假设（设置 \code{bootstrap > 0} 开启）；
#' （2）\strong{VIF}：对每个预测变量以其余预测变量做辅助线性回归，取
#' 1/(1-R²) 作为该变量的方差膨胀因子，VIF>5 通常提示多重共线性；
#' （3）\strong{四张诊断图}：残差 vs 拟合（看非线性/异方差）、正态 Q-Q（看残差
#' 正态性）、尺度-位置（看离散度是否随拟合值变化）、Cook 距离（看强影响点）。
#' 参数：\code{data} 为数据框；\code{formula} 为模型公式（如 \code{y ~ x1 + x2}）；
#' \code{bootstrap} 为重抽样次数（默认 0 即不计算自助法 CI）；\code{conf_level}
#' 为置信水平，默认 0.95。
#' 结果：返回列表——\code{coef}（系数表：估计/标准误/t 值/p 值/正态 CI，若开启
#' 自助法另含 boot 下限/上限）、\code{performance}（R²/调整 R²/F 统计量/AIC/BIC）、
#' \code{vif}（各预测变量 VIF）、\code{diagnostics}（诊断用数据框）、
#' \code{diag_plots}（四张 ggplot 对象列表）、\code{plot}（残差诊断主图）。
#' 注意：VIF 仅对数值型预测变量计算；自助法会显著增加运行时间，样本量很大时
#' 建议先用少量（如 200）试探。本函数仅依赖基础 R 与 ggplot2，无需额外安装包。
#' @param data 数据框（a `data.frame`）。
#' @param formula 回归公式，如 `y ~ x1 + x2`（a model `formula`）。
#' @param bootstrap 自助法重抽样次数；`0` 表示不计算自助法置信区间
#'   （number of bootstrap resamples; `0` = skip bootstrap CI）。
#' @param conf_level 置信水平，默认 `0.95`（confidence level; default `0.95`）。
#' @return 列表，含 `coef` / `performance` / `vif` / `diagnostics` /
#'   `diag_plots` / `plot`（a list with coefficient table, performance metrics,
#'   VIF, diagnostics, diagnostic plots and the main residual plot）。
#' @examples
#' # 第 5 章 5.2/5.6 节：以母亲年龄、孕期体重、吸烟预测出生体重
#' birthwt <- MASS::birthwt
#' res <- mc_linear(birthwt, bwt ~ age + lwt + smoke + factor(race),
#'                  bootstrap = 200)
#' res$coef                       # 系数表（含自助法置信区间）
#' res$performance                # R² / 调整 R² / AIC / BIC
#' res$vif                        # 方差膨胀因子（多重共线性诊断）
#' print(res$plot)                # 残差 vs 拟合诊断图
#' @export
mc_linear <- function(data, formula, bootstrap = 0L, conf_level = 0.95) {
  if (!inherits(data, "data.frame"))
    stop("`data` must be a data.frame.")
  if (!inherits(formula, "formula"))
    stop("`formula` must be a formula.")
  bootstrap <- as.integer(bootstrap)
  if (bootstrap < 0L) stop("`bootstrap` must be a non-negative integer.")

  fit <- stats::lm(formula, data = data)
  sm  <- summary(fit)
  cf  <- as.data.frame(sm$coefficients)
  names(cf) <- c("estimate", "std_error", "t_value", "p_value")
  cf$term <- rownames(cf)
  ci <- stats::confint(fit, level = conf_level)
  cf$ci_lower <- ci[, 1]
  cf$ci_upper <- ci[, 2]

  ## ---- 自研：自助法系数 CI ----
  boot_ci <- NULL
  if (bootstrap > 0L) {
    n <- nrow(data)
    bmat <- matrix(NA_real_, bootstrap, length(coef(fit)))
    colnames(bmat) <- names(coef(fit))
    for (b in seq_len(bootstrap)) {
      idx <- sample.int(n, n, replace = TRUE)
      bf  <- tryCatch(stats::lm(formula, data = data[idx, , drop = FALSE]),
                      error = function(e) NULL)
      if (!is.null(bf)) bmat[b, names(coef(bf))] <- coef(bf)
    }
    a <- (1 - conf_level) / 2
    boot_ci <- apply(bmat, 2, function(v) {
      v <- v[!is.na(v)]
      if (length(v) < 2) c(NA_real_, NA_real_)
      else stats::quantile(v, c(a, 1 - a), na.rm = TRUE)
    })
    boot_ci <- as.data.frame(t(boot_ci))
    names(boot_ci) <- c("boot_lower", "boot_upper")
    cf <- cbind(cf, boot_ci[cf$term, , drop = FALSE])
  }
  cf <- cf[, c("term", "estimate", "std_error", "t_value", "p_value",
               "ci_lower", "ci_upper",
               if (bootstrap > 0L) c("boot_lower", "boot_upper") else NULL)]

  ## ---- 自研：VIF（基于辅助回归 1/(1-R^2)）----
  mm <- stats::model.matrix(formula, data = data)
  preds <- colnames(mm)[colnames(mm) != "(Intercept)"]
  vif <- data.frame(term = preds, VIF = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_along(preds)) {
    others <- setdiff(preds, preds[i])
    if (length(others) == 0) { vif$VIF[i] <- 1; next }
    r2 <- tryCatch({
      ar <- stats::lm(mm[, preds[i]] ~ mm[, others])
      summary(ar)$r.squared
    }, error = function(e) NA_real_)
    vif$VIF[i] <- if (is.na(r2) || r2 >= 1) NA_real_ else 1 / (1 - r2)
  }

  ## ---- 性能量 ----
  performance <- data.frame(
    metric = c("r_squared", "adj_r_squared", "f_statistic", "aic", "bic"),
    value  = c(sm$r.squared, sm$adj.r.squared,
               sm$fstatistic["value"], stats::AIC(fit), stats::BIC(fit)),
    stringsAsFactors = FALSE)

  ## ---- 诊断数据 ----
  fitted <- as.numeric(fit$fitted.values)
  resid  <- as.numeric(fit$residuals)
  stdres <- as.numeric(stats::rstandard(fit))
  cooks  <- as.numeric(stats::cooks.distance(fit))
  diag_df <- data.frame(fitted = fitted, resid = resid,
                        std_resid = stdres, cooks = cooks,
                        stringsAsFactors = FALSE)

  ## ---- 四张诊断图（ggplot2）----
  p1 <- ggplot2::ggplot(diag_df, ggplot2::aes(fitted, resid)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::geom_smooth(se = FALSE, color = "blue", size = 0.8) +
    ggplot2::labs(x = "Fitted values", y = "Residuals") + mc_theme()
  p2 <- ggplot2::ggplot(diag_df, ggplot2::aes(sample = std_resid)) +
    ggplot2::stat_qq(size = 1.2, alpha = 0.6) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red",
                         linetype = "dashed") +
    ggplot2::labs(x = "Theoretical quantiles", y = "Standardized residuals") +
    mc_theme()
  p3 <- ggplot2::ggplot(diag_df, ggplot2::aes(fitted, sqrt(abs(std_resid)))) +
    ggplot2::geom_point(size = 1.2, alpha = 0.6) +
    ggplot2::geom_smooth(se = FALSE, color = "blue", size = 0.8) +
    ggplot2::labs(x = "Fitted values", y = "Sqrt(|standardized residuals|)") +
    mc_theme()
  p4 <- ggplot2::ggplot(diag_df, ggplot2::aes(seq_along(cooks), cooks)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    ggplot2::labs(x = "Observation index", y = "Cook's distance") + mc_theme()

  list(coef = cf, performance = performance, vif = vif,
       diagnostics = diag_df, diag_plots = list(p1, p2, p3, p4),
       plot = p1, fit = fit)
}
