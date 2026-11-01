# =============================================================================
# 软件名称    : macodgent 医学数据科学研究工具包
# 著作权人    : Dr.Xie
# 版本号      : 1.0.0
# 开发语言    : R
# 开发完成日期: 2026-10-30
# 权利声明    : 本程序为著作权人独立开发，受《计算机软件保护条例》保护；
#              未经书面许可，不得复制、发行、出租或将其用于商业目的。
# -----------------------------------------------------------------------------
# 模块        : 定量描述与假设检验模块（desc-quant）
# 本文件包含  : .desc_one（内部）/ mc_desc_quant / mc_norm_test / mc_var_test /
#              mc_ttest / mc_anova / mc_nonparam
#              其中 .desc_one 自行实现偏度、峰度与 Shapiro-Wilk 正态性统计，
#              mc_var_test 的 Levene 法亦为基于 lm 的原创实现。
# =============================================================================

# Internal: compute a set of descriptive statistics for one variable
# (Fisher-Pearson skewness / kurtosis)
.desc_one <- function(x, conf = 0.95) {
  n_missing <- sum(is.na(x))
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3) {
    return(setNames(rep(NA_real_, 18),
      c("n", "n_missing", "mean", "median", "sd", "var", "min", "max",
        "range", "IQR", "SE.mean", "CI.mean.low", "CI.mean.high", "CV",
        "skewness", "kurtosis", "shapiro.W", "shapiro.p")))
  }
  m <- mean(x); s <- sd(x)
  se <- s / sqrt(n)
  qt <- stats::qt(1 - (1 - conf) / 2, df = n - 1)
  z <- (x - m) / s
  g1 <- (n / ((n - 1) * (n - 2))) * sum(z^3)
  g2 <- ((n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3))) * sum(z^4) -
    3 * (n - 1)^2 / ((n - 2) * (n - 3))
  sh <- tryCatch(stats::shapiro.test(x), error = function(e) NULL)
  c(
    n = n, n_missing = n_missing, mean = m, median = stats::median(x),
    sd = s, var = stats::var(x), min = min(x), max = max(x),
    range = max(x) - min(x), IQR = stats::IQR(x),
    SE.mean = se, CI.mean.low = m - qt * se, CI.mean.high = m + qt * se,
    CV = if (m != 0) s / abs(m) else NA_real_,
    skewness = g1, kurtosis = g2,
    shapiro.W = if (is.null(sh)) NA_real_ else as.numeric(sh$statistic),
    shapiro.p = if (is.null(sh)) NA_real_ else as.numeric(sh$p.value)
  )
}

#' 一站式定量描述（旗舰函数）/ One-shot quantitative description
#'
#' 对应第 3 章：一次调用产出约 20 个描述统计量（n、缺失数、均值、中位数、标准差、
#' 方差、最小值、最大值、极差、IQR、均值标准误、95% 置信区间、变异系数、偏度、
#' 峰度、Shapiro-Wilk 正态性检验 W 与 p 值）。支持按分组变量分别计算。
#'
#' Corresponds to Chapter 3 of Part 2: produces about 20 descriptive statistics
#' in one call (n, n missing, mean, median, sd, variance, min, max, range, IQR,
#' SE, 95% CI, CV, skewness, kurtosis, Shapiro-Wilk normality). Supports grouped
#' computation.
#'
#' @details
#' 用途：一键生成定量变量的描述性统计，覆盖约 20 个常用指标，是论文「基线表」和建模前探索的主力函数。
#' 参数：.data 为数据框；vars 指定要描述的数值列（默认全部数值列）；group 可指定分组列，按组分别输出；conf 为均值的置信水平，默认 0.95。
#' 结果：返回数据框，含样本量、均值、标准差、中位数、四分位距、最小值、最大值、偏度、峰度、正态性检验统计量等，以及均值置信区间。
#' 注意：偏度/峰度用于判断分布形态（偏度远偏离 0、峰度远大于 3 提示非正态）；若计划做 t 检验/方差分析，应先看正态性与方差齐性。分组时各组样本量可能不同，注意比较的公平性。
#' @param .data 数据框（a `data.frame`）。
#' @param vars 数值变量名；`NULL` 表示全部数值列（numeric variable names;
#'   `NULL` = all numeric columns）。
#' @param group 可选的分组变量名（optional grouping variable name）。
#' @param conf 置信水平，默认 `0.95`（confidence level; default `0.95`）。
#' @return 宽格式数据框，每行一个变量（分组时每组一行）。列包含 `variable`、
#'   （分组时的 `group`）、以及上述 18 个统计量（a wide `data.frame`, one row
#'   per variable, one row per group level when grouped; 18 statistic columns）。
#' @examples
#' # 第 3 章：低出生体重研究（birthwt，MASS 包）
#' birthwt <- MASS::birthwt
#' cat("==== 总体描述：出生体重/年龄/孕期体重 ====\n")
#' print(mc_desc_quant(birthwt, vars = c("bwt","age","lwt")), row.names = FALSE)
#' cat("\n==== 按母亲是否吸烟分组 ====\n")
#' print(mc_desc_quant(birthwt, vars = c("bwt","age"), group = "smoke"), row.names = FALSE)
#' @export
mc_desc_quant <- function(.data, vars = NULL, group = NULL, conf = 0.95) {
  if (!is.data.frame(.data)) stop("`.data` must be a data.frame.")
  if (is.null(vars)) vars <- names(.data)[sapply(.data, is.numeric)]
  stopifnot(all(vars %in% names(.data)))

  rows <- list()
  if (is.null(group)) {
    for (v in vars) rows[[length(rows) + 1]] <- c(variable = v, .desc_one(.data[[v]], conf))
  } else {
    gvec <- .data[[group]]
    for (lv in unique(gvec)) {
      sub <- .data[gvec == lv, vars, drop = FALSE]
      for (v in vars) {
        one <- .desc_one(sub[[v]], conf)
        row <- c(setNames(as.character(lv), group), variable = v, one)
        rows[[length(rows) + 1]] <- row
      }
    }
  }
  mat <- do.call(rbind, rows)
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  for (cn in names(df)) if (!cn %in% c("variable", group)) {
    df[[cn]] <- as.numeric(df[[cn]])
  }
  if (!is.null(group)) df <- df[, c(group, setdiff(names(df), group))]
  df
}

#' 正态性检验 / Normality test
#'
#' 对应第 3 章：默认 Shapiro-Wilk 检验（n <= 5000 稳健），或可选
#' Kolmogorov-Smirnov 检验。返回检验名称、统计量与 p 值。
#'
#' Corresponds to Chapter 3: Shapiro-Wilk (default, robust for n <= 5000) or
#' Kolmogorov-Smirnov. Returns the test name, statistic and p-value.
#'
#' @details
#' 用途：检验一个定量变量是否来自正态分布，是选择参数/非参数方法的前置判断。
#' 方法：shapiro（默认）为 Shapiro-Wilk 检验，对小到中等样本（通常 n<5000）敏感且功效好；ks 为 Kolmogorov-Smirnov 检验，默认与标准正态分布比较，对位置/尺度敏感。
#' 参数：x 为待检验的数值向量；method 选择上述两种方法。
#' 结果：返回检验统计量、p 值等。
#' 注意：大样本下正态性检验极易「显著」（p 很小）即便偏离程度实际无害，因此不要仅凭 p 值下结论，应结合直方图、Q-Q 图与偏度/峰度综合判断。Shapiro 检验在 n>5000 时可能不可靠或极慢。
#' @param x 数值向量（a numeric vector）。
#' @param method `"shapiro"`（默认）或 `"ks"`（test method）。
#' @param ... 传给 \code{stats::ks.test()} 的额外参数（extra arguments passed
#'   to \code{stats::ks.test()}）。
#' @return 数据框，列为 `test` / `statistic` / `p.value`（a `data.frame` with
#'   columns `test`, `statistic`, `p.value`）。
#' @examples
#' # 第 3 章：检验 birthwt 中定量变量是否符合正态
#' birthwt <- MASS::birthwt
#' mc_norm_test(birthwt$bwt)         # 出生体重
#' mc_norm_test(birthwt$age)         # 母亲年龄
#' @export
mc_norm_test <- function(x, method = c("shapiro", "ks"), ...) {
  method <- match.arg(method)
  x <- x[!is.na(x)]
  if (method == "shapiro") {
    t <- stats::shapiro.test(x)
    data.frame(test = "Shapiro-Wilk", statistic = as.numeric(t$statistic),
               p.value = as.numeric(t$p.value), stringsAsFactors = FALSE)
  } else {
    t <- stats::ks.test(x, "pnorm", mean(x), stats::sd(x), ...)
    data.frame(test = "Kolmogorov-Smirnov", statistic = as.numeric(t$statistic),
               p.value = as.numeric(t$p.value), stringsAsFactors = FALSE)
  }
}

#' 方差齐性检验 / Variance homogeneity test
#'
#' 对应第 3 章「三种方差齐性检验」：Bartlett、Fligner-Killeen、以及
#' 基于均值离差的 Levene 法。返回检验名称、统计量与 p 值。
#'
#' Corresponds to Chapter 3 ("three methods for variance homogeneity"):
#' Bartlett, Fligner-Killeen, and Levene (mean-deviation method). Returns the
#' test name, statistic and p-value.
#'
#' @details
#' 用途：检验多组数据的总体方差是否相等（方差齐性），是选择方差分析/t 检验前的必要前提（它们通常假设方差齐）。
#' 方法：bartlett（默认）对正态数据最敏感；fligner（Fligner-Killeen）对偏离正态更稳健；levene（Levene）基于离差，对厚尾分布稳健，最常用。
#' 参数：x 为数值向量；g 为分组因子/向量（长度须与 x 相同）；method 选择上述三种之一。
#' 结果：返回检验统计量、自由度与 p 值。
#' 注意：若方差不齐（p 很小），应改用 Welch 校正（如 mc_ttest 默认对两样本做 Welch）或非参数方法，而非标准 t 检验/ANOVA。Bartlett 对正态性假设较强，非正态数据优先 Levene 或 Fligner。
#' @param x 数值向量（a numeric vector）。
#' @param g 分组因子（a grouping factor）。
#' @param method `"bartlett"`（默认）/ `"fligner"` / `"levene"`（test method）。
#' @return 数据框，列为 `test` / `statistic` / `p.value`（a `data.frame` with
#'   columns `test`, `statistic`, `p.value`）。
#' @examples
#' # 第 3 章：按种族(race)检验不同组出生体重的方差齐性
#' birthwt <- MASS::birthwt
#' mc_var_test(birthwt$bwt, birthwt$race, method = "bartlett")
#' mc_var_test(birthwt$bwt, birthwt$race, method = "levene")
#' @export
mc_var_test <- function(x, g, method = c("bartlett", "fligner", "levene")) {
  method <- match.arg(method)
  g <- as.factor(g)
  if (method == "bartlett") {
    t <- stats::bartlett.test(x, g)
    nm <- "Bartlett"
  } else if (method == "fligner") {
    t <- stats::fligner.test(x, g)
    nm <- "Fligner-Killeen"
  } else {
    means <- tapply(x, g, mean, na.rm = TRUE)
    d <- abs(x - means[as.character(g)])
    a <- stats::anova(stats::lm(d ~ g))
    return(data.frame(test = "Levene",
                      statistic = as.numeric(a["g", "F value"]),
                      p.value = as.numeric(a["g", "Pr(>F)"]),
                      stringsAsFactors = FALSE))
  }
  data.frame(test = nm, statistic = as.numeric(t$statistic),
             p.value = as.numeric(t$p.value), stringsAsFactors = FALSE)
}

#' t 检验（统一入口）/ t-test (unified interface)
#'
#' 对应第 3 章：单样本、独立两样本、以及配对 t 检验的统一入口。通过
#' 是否提供 `y` 与 `paired` 自动分流。
#'
#' Corresponds to Chapter 3: one-sample, independent two-sample, and paired
#' t-tests through a single entry point.
#'
#' @details
#' 用途：统一入口做 t 检验，覆盖三种常见情形：单样本（与已知常数比较）、独立两样本、配对样本。
#' 方法：两独立样本默认采用 Welch 校正的 t 检验（不要求方差齐），更符合实际；可通过参数切换为 student t。
#' 参数：x 为第一组或单样本向量；y 为第二组向量（独立/配对时用，单样本留空）；paired 为 TRUE 时做配对 t 检验；mu 为单样本检验的参照值，默认 0；conf/alternative 控制置信水平与备择方向。
#' 结果：返回均值差（或样本均值）、标准误、t 统计量、自由度、p 值与置信区间。
#' 注意：配对检验要求 x、y 一一对应（如同一人前后测量）；若两组不独立却用了独立检验会低估误差。样本明显偏态或小样本时，考虑 mc_nonparam() 的 Wilcoxon 检验。
#' @param x 数值向量（a numeric vector）。
#' @param y 第二个向量（独立/配对时用）；单样本时为 `NULL`（a second vector
#'   for independent/paired; `NULL` for one-sample）。
#' @param paired 是否配对，默认 `FALSE`（logical; paired? default `FALSE`）。
#' @param mu 单样本检验的假设均值，默认 `0`（hypothesised mean for one-sample;
#'   default `0`）。
#' @param conf 置信水平，默认 `0.95`（confidence level; default `0.95`）。
#' @param alternative `"two.sided"`（默认）/ `"less"` / `"greater"`（alternative
#'   hypothesis）。
#' @return 数据框，列为 `method` / `estimate` / `conf.low` / `conf.high` /
#'   `statistic` / `p.value`（a `data.frame` with those six columns）。
#' @examples
#' # 第 3 章
#' birthwt <- MASS::birthwt
#' # 单样本：平均出生体重是否等于 3000 g？
#' mc_ttest(birthwt$bwt, mu = 3000)
#' # 独立两样本（Welch 校正）：吸烟 vs 不吸烟母亲的出生体重
#' mc_ttest(birthwt$bwt[birthwt$smoke == 1], birthwt$bwt[birthwt$smoke == 0])
#' @export
mc_ttest <- function(x, y = NULL, paired = FALSE, mu = 0, conf = 0.95,
                     alternative = "two.sided") {
  if (is.null(y) && !paired) {
    tt <- stats::t.test(x, mu = mu, conf.level = conf, alternative = alternative)
  } else {
    tt <- stats::t.test(x, y, paired = paired, conf.level = conf, alternative = alternative)
  }
  est <- tt$estimate
  if (length(est) == 2) est <- est[1] - est[2]
  data.frame(method = tt$method, estimate = as.numeric(est),
             conf.low = as.numeric(tt$conf.int[1]),
             conf.high = as.numeric(tt$conf.int[2]),
             statistic = as.numeric(tt$statistic),
             p.value = as.numeric(tt$p.value), stringsAsFactors = FALSE)
}

#' 单因素方差分析 / One-way ANOVA
#'
#' 对应第 3 章「方差分析（ANOVA）」。基于 \code{stats::aov}，返回各项的
#' 自由度、平方和、均方、F 统计量与 p 值。
#'
#' Corresponds to Chapter 3 ("ANOVA"). Built on \code{stats::aov}, returns the
#' degrees of freedom, sum of squares, mean square, F statistic and p-value of
#' each term.
#'
#' @details
#' 用途：单因素方差分析（one-way ANOVA），比较三个及以上组的均值是否总体相等。
#' 方法：基于线性模型分解变异，比较组间与组内均方，得到 F 统计量。
#' 参数：formula 形如 y ~ g（因变量 ~ 分组因子）；data 为包含这两列的数据框。
#' 结果：返回各项的自由度、平方和、均方、F 值与 p 值。
#' 注意：ANOVA 只告诉你「至少两组不同」，具体哪两组不同需事后多重比较（如 Tukey HSD）。其前提假设为各组正态且方差齐——先用 mc_norm_test()、mc_var_test() 核验。若方差不齐或分布偏态，考虑 Kruskal-Wallis（见 mc_nonparam）。
#' @param formula 公式，如 \code{y ~ group}（a formula like \code{y ~ group}）。
#' @param data 包含变量的数据框（a `data.frame` containing the variables）。
#' @return 数据框，列为 `term` / `df` / `sumsq` / `meansq` / `statistic` /
#'   `p.value`（a `data.frame` with those six columns）。
#' @examples
#' # 第 3 章：不同种族(race)母亲的出生体重是否有差异（单因素 ANOVA）
#' birthwt <- MASS::birthwt
#' mc_anova(bwt ~ factor(race), data = birthwt)
#' @export
mc_anova <- function(formula, data) {
  a <- stats::aov(formula, data = data)
  s <- summary(a)[[1]]
  data.frame(term = rownames(s), df = s[, "Df"], sumsq = s[, "Sum Sq"],
             meansq = s[, "Mean Sq"], statistic = s[, "F value"],
             p.value = s[, "Pr(>F)"], stringsAsFactors = FALSE)
}

#' 非参数组间差异检验 / Non-parametric group difference test
#'
#' 对应第 3 章「组间非参数检验」：两组用 Wilcoxon 秩和检验，多组用
#' Kruskal-Wallis 检验。返回检验名称、统计量与 p 值。
#'
#' Corresponds to Chapter 3 ("non-parametric tests between groups"): Wilcoxon
#' rank-sum (two groups) or Kruskal-Wallis (multiple groups). Returns the test
#' name, statistic and p-value.
#'
#' @details
#' 用途：当数据不满足正态/方差齐假设，或仅为等级资料时，用非参数方法比较组间差异。
#' 方法：wilcoxon 为 Wilcoxon 秩和检验（Mann-Whitney U），用于两组独立样本；kruskal 为 Kruskal-Wallis 检验，用于三组及以上。它们比较的是分布位置（秩）而非均值。
#' 参数：x 为数值向量；g 为分组向量（两组用 wilcoxon，多组用 kruskal）；method 选择方法。
#' 结果：返回检验统计量、p 值（Kruskal 还含自由度）。
#' 注意：非参数检验功效通常低于对应的参数检验（数据确属正态时），但更稳健；它检验「分布是否相同/位置是否偏移」，结论表述应为「中位数/分布差异」而非「均值差异」。两组配对请用 Wilcoxon 符号秩检验（此函数未覆盖，可用 stats::wilcox.test(paired=TRUE)）。
#' @param x 数值向量（a numeric vector）。
#' @param g 分组向量/因子（grouping vector/factor）。
#' @param method `"wilcoxon"`（默认）或 `"kruskal"`（test method）。
#' @return 数据框，列为 `test` / `statistic` / `p.value`（a `data.frame` with
#'   columns `test`, `statistic`, `p.value`）。
#' @examples
#' # 第 3 章：出生体重按吸烟状态（两组 Wilcoxon）与种族（多组 Kruskal-Wallis）
#' birthwt <- MASS::birthwt
#' mc_nonparam(birthwt$bwt, birthwt$smoke, method = "wilcoxon")
#' mc_nonparam(birthwt$bwt, factor(birthwt$race), method = "kruskal")
#' @export
mc_nonparam <- function(x, g, method = c("wilcoxon", "kruskal")) {
  method <- match.arg(method)
  g <- as.factor(g)
  if (method == "wilcoxon") {
    t <- stats::wilcox.test(x ~ g)
    nm <- "Wilcoxon rank-sum"
  } else {
    t <- stats::kruskal.test(x ~ g)
    nm <- "Kruskal-Wallis"
  }
  data.frame(test = nm, statistic = as.numeric(t$statistic),
             p.value = as.numeric(t$p.value), stringsAsFactors = FALSE)
}
