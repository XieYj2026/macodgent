# =============================================================================
# 软件名称    : macodgent 医学数据科学研究工具包
# 著作权人    : Dr.Xie
# 版本号      : 1.0.0
# 开发语言    : R
# 开发完成日期: 2026-10-30
# 权利声明    : 本程序为著作权人独立开发，受《计算机软件保护条例》保护；
#              未经书面许可，不得复制、发行、出租或将其用于商业目的。
# -----------------------------------------------------------------------------
# 模块        : 定性描述与关联分析模块（desc-qual）
# 本文件包含  : mc_contingency / mc_corr / mc_tableone
#              mc_contingency 覆盖卡方、Fisher 精确与 CMH 分层三种独立性检验；
#              mc_corr 的偏相关通过对控制变量回归取残差后计算，为原创实现。
# =============================================================================

#' 列联表与独立性检验（卡方 / Fisher / CMH）/ Contingency table & independence test
#'
#' 对应第 4 章 4.3 节「独立性检验，卡方检验，费希尔精确概率检验和
#' Cochran-Mantel-Haenszel 检验」。可构建两个分类变量的列联表并做三种检验：
#' ① `"chisq"`（默认）卡方独立性检验，对 2×2 表额外用 Fisher 法给出 OR 与 95% CI；
#' ② `"fisher"` 对任意 r×c 表做 Fisher 精确检验；
#' ③ `"cmh"` 在给定分层变量 `strata` 下做 Cochran-Mantel-Haenszel 分层检验，
#' 控制混杂后得到共同的（合并）比值比。
#'
#' Corresponds to Chapter 4, Section 4.3 (chi-square, Fisher's exact, and the
#' Cochran-Mantel-Haenszel test). Builds a contingency table and runs one of
#' three independence tests selected by `method`.
#'
#' @details
#' 用途：构建分类变量的列联表并做独立性/关联检验，是分析分类变量关联的第一步，且覆盖 4.3 节要求的三种方法。
#' 方法：
#' \itemize{
#'   \item \code{"chisq"}（默认）：对 2×2 及以上表做卡方独立性检验；\code{correct=TRUE} 时对 2×2 表加 Yates 连续性校正，防止小样本高估显著性。当表为 2×2 时，额外用 Fisher 精确法给出 OR 及其 95% CI。
#'   \item \code{"fisher"}：对任意 r×c 表用 \code{stats::fisher.test} 做精确检验（适合小样本或期望频数过小的情形）。
#'   \item \code{"cmh"}：需提供 \code{strata} 分层变量，构造三维列联表后用 \code{stats::mantelhaen.test} 做 Cochran-Mantel-Haenszel 检验，在控制分层因素后判断行-列关联是否仍然存在，并返回共同比值比（common OR）。
#' }
#' 参数：.data 为数据框；row/col 为行、列变量名；correct 控制卡方连续性校正；method 选择检验类型；strata 为 CMH 所需的层变量名。
#' 结果：返回列表，含 table（列联表，三维时分层展开）、test（检验数据框：test/statistic/df/p.value），并视方法附带 OR / commonOR。
#' 注意：卡方检验要求期望频数不能过小（一般要求足够比例的格子期望 >5），否则改用 \code{method="fisher"}。CMH 适用于存在混杂分层（如多中心研究）时估计「调整后的」关联。OR>1 提示正关联；置信区间跨 1 表示关联不显著。
#' @param .data 数据框（a `data.frame`）。
#' @param row 行变量名（分类）（row variable name, categorical）。
#' @param col 列变量名（分类）（column variable name, categorical）。
#' @param correct 卡方连续性校正，默认 `TRUE`（continuity correction for
#'   chi-square; default `TRUE`；仅对 `method="chisq"` 生效）。
#' @param method 检验类型：`"chisq"`（默认）、`"fisher"` 或 `"cmh"`
#'   （test type: chi-square / Fisher exact / Cochran-Mantel-Haenszel）。
#' @param strata CMH 分层变量名（仅 `method="cmh"` 时使用）（stratifying
#'   variable name for the CMH test）。
#' @return 列表，含 `table`（列联表）、`test`（检验数据框：test/statistic/df/
#'   p.value）；`method="chisq"` 且 2×2 表时附 `OR`（OR/CI.low/CI.high/fisher.p），
#'   `method="fisher"` 时附 `OR`（含 95% CI 与 p），`method="cmh"` 时附
#'   `commonOR`（共同比值比及其 CI）（a list with \code{table}, \code{test},
#'   and method-specific OR / commonOR）。
#' @examples
#' # 第 4 章 4.3 节：以 birthwt 演示三种独立性检验
#' birthwt <- MASS::birthwt
#'
#' ## ---- 卡方独立性检验（2x2 附 Fisher OR）----
#' res <- mc_contingency(birthwt, row = "smoke", col = "low")
#' cat("---- 列联表 / table ----\n"); print(res$table)
#' cat("\n---- 卡方检验 / chi-square ----\n"); print(res$test, row.names = FALSE)
#' cat("\n---- 比值比 OR / odds ratio ----\n"); print(res$OR, row.names = FALSE)
#'
#' ## ---- Fisher 精确检验（小样本/期望频数过小时）----
#' res2 <- mc_contingency(birthwt, row = "smoke", col = "low", method = "fisher")
#' print(res2$test, row.names = FALSE); print(res2$OR, row.names = FALSE)
#'
#' ## ---- CMH 分层检验：控制高血压(ht)后吸烟与低体重的共同 OR ----
#' res3 <- mc_contingency(birthwt, row = "smoke", col = "low",
#'                        method = "cmh", strata = "ht")
#' cat("\n---- CMH 检验 / Cochran-Mantel-Haenszel ----\n")
#' print(res3$test, row.names = FALSE)
#' print(res3$commonOR, row.names = FALSE)
#' @export
mc_contingency <- function(.data, row, col, correct = TRUE,
                           method = c("chisq", "fisher", "cmh"),
                           strata = NULL) {
  method <- match.arg(method)
  if (method == "cmh") {
    if (is.null(strata)) stop("`method = 'cmh'` requires a `strata` variable.")
    if (!requireNamespace("survival", quietly = TRUE))
      stop("`method = 'cmh'` requires the survival package.")
    xtab <- stats::xtabs(~ .data[[row]] + .data[[col]] + .data[[strata]])
    mh <- stats::mantelhaen.test(xtab)
    test <- data.frame(test = "Cochran-Mantel-Haenszel",
                       statistic = as.numeric(mh$statistic),
                       df = as.numeric(mh$parameter),
                       p.value = as.numeric(mh$p.value),
                       stringsAsFactors = FALSE)
    ## 共同比值比（MH 点估计），其 95% CI 用「精确条件 Logistic」
    ## （survival::coxph, method="exact"，等价于 clogit）求得，比 MH 渐进方差
    ## 更稳健；survival 已在 Imports 中。
    cor <- as.numeric(mh$estimate)
    dn <- dimnames(xtab)
    rows <- list()
    for (kk in seq_len(dim(xtab)[3])) {
      t2 <- xtab[, , kk]
      a <- t2[1, 1]; b <- t2[1, 2]; c <- t2[2, 1]; d <- t2[2, 2]
      rows[[length(rows) + 1]] <- data.frame(stratum = dn[[3]][kk], event = 1, exposure = 1, n = a)
      rows[[length(rows) + 1]] <- data.frame(stratum = dn[[3]][kk], event = 0, exposure = 1, n = b)
      rows[[length(rows) + 1]] <- data.frame(stratum = dn[[3]][kk], event = 1, exposure = 0, n = c)
      rows[[length(rows) + 1]] <- data.frame(stratum = dn[[3]][kk], event = 0, exposure = 0, n = d)
    }
    long <- do.call(rbind, rows)
    long <- long[rep(seq_len(nrow(long)), long$n), c("stratum", "event", "exposure")]
    long$stratum <- factor(long$stratum)
    ci <- tryCatch({
      cl <- survival::coxph(
        survival::Surv(rep(1, nrow(long)), long$event) ~ long$exposure +
          survival::strata(long$stratum), method = "exact")
      as.numeric(stats::confint(cl))
    }, error = function(e) c(NA_real_, NA_real_))
    res <- list(
      table = xtab,
      test = test,
      commonOR = data.frame(commonOR = cor,
                            CI.low = if (length(ci) == 2) exp(ci[1]) else NA_real_,
                            CI.high = if (length(ci) == 2) exp(ci[2]) else NA_real_,
                            stringsAsFactors = FALSE))
    return(res)
  }
  tab <- table(.data[[row]], .data[[col]])
  if (method == "fisher") {
    ft <- stats::fisher.test(tab)
    test <- data.frame(test = "Fisher's exact",
                       statistic = NA_real_,
                       df = NA_real_,
                       p.value = as.numeric(ft$p.value),
                       stringsAsFactors = FALSE)
    a <- tab[1, 1]; b <- tab[1, 2]; c <- tab[2, 1]; d <- tab[2, 2]
    or_point <- if (nrow(tab) == 2 && ncol(tab) == 2) (a * d) / (b * c) else NA_real_
    res <- list(
      table = tab,
      test = test,
      OR = data.frame(OR = or_point,
                      CI.low = as.numeric(ft$conf.int[1]),
                      CI.high = as.numeric(ft$conf.int[2]),
                      fisher.p = as.numeric(ft$p.value),
                      stringsAsFactors = FALSE))
    return(res)
  }
  ## default: chi-square
  ct <- stats::chisq.test(tab, correct = correct)
  res <- list(
    table = tab,
    test = data.frame(test = "Chi-square",
                      statistic = as.numeric(ct$statistic),
                      df = as.numeric(ct$parameter),
                      p.value = as.numeric(ct$p.value),
                      stringsAsFactors = FALSE)
  )
  if (nrow(tab) == 2 && ncol(tab) == 2) {
    a <- tab[1, 1]; b <- tab[1, 2]; c <- tab[2, 1]; d <- tab[2, 2]
    or <- (a * d) / (b * c)
    ft <- stats::fisher.test(tab)
    res$OR <- data.frame(OR = or,
                         CI.low = as.numeric(ft$conf.int[1]),
                         CI.high = as.numeric(ft$conf.int[2]),
                         fisher.p = as.numeric(ft$p.value),
                         stringsAsFactors = FALSE)
  }
  res
}

#' 相关分析（含偏相关）/ Correlation analysis (with partial correlation)
#'
#' 对应第 4 章：Pearson / Spearman 相关；可通过 `control` 指定控制变量，
#' 对所有变量先就控制变量做回归取残差，从而得到偏相关。
#'
#' Corresponds to Chapter 4: Pearson / Spearman correlation; optionally compute
#' partial correlation by residualising all variables on the control variable(s).
#'
#' @details
#' 用途：计算变量间的相关系数，支持 Pearson（线性）与 Spearman（秩/单调），并可计算偏相关（控制混杂变量后）。
#' 方法：pearson 衡量线性关联，假设近似正态；spearman 基于秩，对非线性单调关系与异常值更稳健；control 指定一个或多个控制变量，计算偏相关（剔除这些变量影响后的净相关）。
#' 参数：.data 为数据框；vars 为参与相关的变量名向量；method 选择相关类型；control 为可选的控制变量名。
#' 结果：返回相关矩阵（含相关系数），偏相关时返回控制后的矩阵。
#' 注意：相关系数只衡量关联强度与方向，不等于因果；接近 0 不代表无关（可能是非线性）。Spearman 更适合 ordinal 数据或明显偏态。偏相关能帮助识别「虚假相关」（由共同混杂变量造成）。
#' @param .data 数据框（a `data.frame`）。
#' @param vars 变量名；`NULL` 表示全部数值列（variable names; `NULL` = all
#'   numeric columns）。
#' @param method `"pearson"`（默认）或 `"spearman"`（correlation method）。
#' @param control 可选控制变量名（用于偏相关）（optional control variable names
#'   for partial correlation）。
#' @return 长格式数据框，列为 `var1` / `var2` / `r` / `p.value`（a long
#'   `data.frame` with columns `var1`, `var2`, `r`, `p.value`）。
#' @examples
#' # 第 4 章 4.5 节：birthwt 中出生体重、孕期体重的相关（Spearman）
#' birthwt <- MASS::birthwt
#' cat("==== Spearman 相关 / Spearman ====\n")
#' # 注：birthwt 含相同观测值（ties），cor.test 会以渐近法替代精确 p 值并提示，
#' # 属正常统计说明；此处抑制该提示以简化示例输出
#' suppressWarnings(print(mc_corr(birthwt, vars = c("bwt","lwt"), method = "spearman"),
#'                          row.names = FALSE))
#' cat("\n==== 控制母亲年龄后的偏相关 / partial on age ====\n")
#' print(mc_corr(birthwt, vars = c("bwt","lwt"), control = "age"),
#'       row.names = FALSE)
#' @export
mc_corr <- function(.data, vars = NULL, method = c("pearson", "spearman"),
                    control = NULL) {
  method <- match.arg(method)
  if (is.null(vars)) vars <- names(.data)[sapply(.data, is.numeric)]
  mat <- .data[vars]
  if (!is.null(control)) {
    if (method != "pearson") warning("Partial correlation is based on residuals; method forced to pearson.")
    method <- "pearson"
    for (v in vars) {
      f <- stats::as.formula(paste(v, "~", paste(control, collapse = "+")))
      fit <- stats::lm(f, data = .data)
      mat[[v]] <- stats::residuals(fit)
    }
  }
  n <- length(vars)
  pairs <- utils::combn(n, 2)
  rows <- list()
  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1, k]; j <- pairs[2, k]
    t <- stats::cor.test(mat[[i]], mat[[j]], method = method)
    rows[[length(rows) + 1]] <- data.frame(
      var1 = vars[i], var2 = vars[j],
      r = as.numeric(t$estimate), p.value = as.numeric(t$p.value),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

#' 基线特征表（Table One）/ Baseline characteristics table
#'
#' 对应第 4 章「基线表」。数值变量以「均值 ± 标准差」呈现，分类变量以
#' 频数与百分比呈现，支持按分组变量分层。
#'
#' Corresponds to Chapter 4 ("Table One"). Numeric variables are summarised as
#' mean +/- SD; categorical variables as counts and percentages; can be
#' stratified by a grouping variable.
#'
#' @details
#' 用途：生成论文中常见的「Table 1」基线特征表，一行一个变量，按暴露/分组列展示，便于直接写进 manuscript。
#' 方法：数值变量汇总为「均值 ± 标准差」（近似正态时）；分类变量汇总为频数与百分比。
#' 参数：.data 为数据框；vars 为要纳入的变量名；group 为可选的分组列（如治疗组），此时按组分别汇总；digits 控制小数位数，默认 2。
#' 结果：返回整齐的数据框，可直接 knitr::kable() 或 print() 成表。
#' 注意：分组比较的 p 值通常不在此表内（可另行用 mc_ttest/mc_nonparam/mc_contingency 计算）。分类变量若水平过多会撑宽表格，建议先合并稀有类别。这是描述性表格，不替代统计推断。
#' @param .data 数据框（a `data.frame`）。
#' @param vars 需要纳入的变量名（variable names to include）。
#' @param group 可选的分组变量名（optional grouping variable name）。
#' @param digits 小数位数，默认 `2`（number of decimals; default `2`）。
#' @return 数据框，列为 `variable` / `level` / `stratum` / `value`（a
#'   `data.frame` with columns `variable`, `level`, `stratum`, `value`）。
#' @examples
#' # 第 4 章：birthwt 基线特征表（按低出生体重 low 分组）
#' birthwt <- MASS::birthwt
#' bw <- birthwt
#' bw$race  <- factor(bw$race, levels = 1:3, labels = c("白人","黑人","其他"))
#' bw$smoke <- factor(bw$smoke, levels = 0:1, labels = c("否","是"))
#' bw$low   <- factor(bw$low, levels = 0:1, labels = c("正常","低体重"))
#' mc_tableone(bw, vars = c("age","lwt","race","smoke","bwt"), group = "low")
#' @export
mc_tableone <- function(.data, vars, group = NULL, digits = 2) {
  rows <- list()
  for (v in vars) {
    x <- .data[[v]]
    if (is.numeric(x)) {
      rows[[length(rows) + 1]] <- data.frame(
        variable = v, level = "", stratum = "Overall",
        value = sprintf("%.*f +/- %.*f", digits, mean(x, na.rm = TRUE),
                        digits, stats::sd(x, na.rm = TRUE)),
        stringsAsFactors = FALSE)
      if (!is.null(group)) for (lv in levels(factor(.data[[group]]))) {
        xx <- x[.data[[group]] == lv]
        rows[[length(rows) + 1]] <- data.frame(
          variable = v, level = "", stratum = as.character(lv),
          value = sprintf("%.*f +/- %.*f", digits, mean(xx, na.rm = TRUE),
                          digits, stats::sd(xx, na.rm = TRUE)),
          stringsAsFactors = FALSE)
      }
    } else {
      x <- factor(x)
      tb <- table(x)
      for (lvl in names(tb)) {
        pct <- 100 * tb[lvl] / sum(tb)
        rows[[length(rows) + 1]] <- data.frame(
          variable = v, level = lvl, stratum = "Overall",
          value = sprintf("%d (%.1f%%)", tb[lvl], pct),
          stringsAsFactors = FALSE)
      }
      if (!is.null(group)) for (g in levels(factor(.data[[group]]))) {
        xx <- factor(x[.data[[group]] == g])
        tb <- table(xx)
        for (lvl in names(tb)) {
          pct <- 100 * tb[lvl] / sum(tb)
          rows[[length(rows) + 1]] <- data.frame(
            variable = v, level = lvl, stratum = as.character(g),
            value = sprintf("%d (%.1f%%)", tb[lvl], pct),
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  do.call(rbind, rows)
}
