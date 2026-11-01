# =============================================================================
# 软件名称    : macodgent 医学数据科学研究工具包
# 著作权人    : Dr.Xie
# 版本号      : 1.0.0
# 开发语言    : R
# 开发完成日期: 2026-10-30
# 权利声明    : 本程序为著作权人独立开发，受《计算机软件保护条例》保护；
#              未经书面许可，不得复制、发行、出租或将其用于商业目的。
# -----------------------------------------------------------------------------
# 模块        : 全局变量声明（global variable declarations）
# 本文件包含  : 非标准计算（NSE）引用的变量声明，避免 R CMD check 误报
# =============================================================================

#' @importFrom stats sd median IQR setNames
NULL

# Declare variable names referenced via non-standard evaluation (NSE) to avoid
# spurious "undefined global variable" notes from R CMD check.
utils::globalVariables(c("sv", "g", "d", "z", "variable", "estimate",
                         "conf.low", "conf.high", "strata", "lo", "hi"))
