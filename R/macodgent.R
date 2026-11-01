# =============================================================================
# 软件名称    : macodgent 医学数据科学研究工具包
# 著作权人    : Dr.Xie
# 版本号      : 1.0.0
# 开发语言    : R
# 开发完成日期: 2026-10-30
# 权利声明    : 本程序为著作权人独立开发，受《计算机软件保护条例》保护；
#              未经书面许可，不得复制、发行、出租或将其用于商业目的。
# -----------------------------------------------------------------------------
# 模块        : 包级说明文档（package-level documentation）
# 本文件包含  : macodgent-package 总体说明、模块与函数清单、设计原则
# =============================================================================

#' macodgent：医学数据科学例行性分析工具箱
#'
#' 一套从医学数据科学例行性分析技术提炼出的
#' 可复用分析函数集合，覆盖数据预处理、定量/定性描述统计、常规与进阶回归、
#' 生存分析、竞争风险与多状态模型。所有导出函数均以 `mc_` 为前缀，并返回
#' 整洁的 `data.frame`，可直接落入报告与下游流程。
#'
#' macodgent is a curated toolkit of reusable analysis functions distilled from
#' Part 2 ("Routine Analysis Techniques") of the routine medical-data-science tutorial.
#' It covers data preprocessing, quantitative and qualitative
#' descriptive statistics, conventional and advanced regression, survival
#' analysis, competing-risk and multi-state models. All exported functions share
#' the `mc_` prefix and return tidy `data.frame` objects that drop straight into
#' reports and downstream pipelines.
#'
#' @section 模块与函数 / Modules and functions:
#' \describe{
#'   \item{第 1 章 认识数据 / Meet the data}{`mc_mode()`, `mc_inspect()`,
#'     `mc_theme()`, `mc_palette()`}
#'   \item{第 2 章 预处理 / Preprocessing}{`mc_impute_knn()`,
#'     `mc_impute_mice()`, `mc_dedup()`, `mc_scale()`, `mc_wide_long()`,
#'     `mc_long_wide()`}
#'   \item{第 3 章 定量描述 / Quantitative}{`mc_desc_quant()`,
#'     `mc_norm_test()`, `mc_var_test()`, `mc_ttest()`, `mc_anova()`,
#'     `mc_nonparam()`}
#'   \item{第 4 章 定性描述 / Qualitative}{`mc_contingency()`, `mc_corr()`,
#'     `mc_tableone()`}
#'   \item{第 5 章 常规回归 / Conventional regression}{`mc_linear()`,
#'     `mc_logistic()`, `mc_logistic_screen()`, `mc_stepwise()`,
#'     `mc_ordinal_logistic()`, `mc_clogistic()`, `mc_model_compare()`}
#'   \item{第 6 章 生存分析 / Survival}{`mc_survfit()`, `mc_cox()`,
#'     `mc_ph_test()`}
#'   \item{第 7 章 进阶回归 / Advanced regression}{`mc_gam()`, `mc_lasso()`,
#'     `mc_ridge()`, `mc_elastic_net()`, `mc_rcs()`, `mc_quantile_reg()`,
#'     `mc_pcr()`, `mc_nnet_reg()`}
#'   \item{第 8 章 进阶生存 / Advanced survival}{`mc_competing_risk()`,
#'     `mc_multistate()`}
#' }
#'
#' @section 设计原则 / Design principles:
#' \enumerate{
#'   \item 所有导出函数以 `mc_` 前缀命名，避免与基础 R 及其它包冲突。
#'   \item 分析函数一律返回整洁 `data.frame`（而非仅打印输出），便于报告与管道。
#'   \item 重方法（GAM、Lasso、竞争风险、多状态等）放在 `Suggests`，通过
#'     `requireNamespace()` 优雅降级；核心仅依赖 `survival`、
#'     `ggplot2`、`broom`、`stats`、`utils`。
#'   \item All exported functions use the `mc_` prefix to avoid name clashes.
#'   \item Analysis functions return tidy `data.frame` objects, ready for
#'     reporting and piping.
#'   \item Heavy methods live in `Suggests` with graceful namespace guards.
#' }
#'
#' @docType package
#' @name macodgent-package
#' @aliases macodgent
#' @keywords internal
#' @seealso 入口函数 / Entry points: \code{\link{mc_inspect}},
#'   \code{\link{mc_desc_quant}}, \code{\link{mc_tableone}},
#'   \code{\link{mc_logistic}}, \code{\link{mc_cox}}.
NULL
