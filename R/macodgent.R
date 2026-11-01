# =============================================================================
# Software name: macodgent medical data science research toolkit
# Copyright holder: Dr.Xie
# Version: 1.0.0
# Language: R
# Date completed: 2026-10-30
# Rights: This program was independently developed by the copyright holder and is
# protected under the Computer Software Protection Regulations; without
# written permission it may not be copied, distributed, rented, or used
# for commercial purposes.
# -----------------------------------------------------------------------------
# Module: package-level documentation
# This file contains: macodgent-package overview, module & function list,
# design principles
# =============================================================================

#' Macodgent: routine-analysis toolkit for medical data science
#'
#' Macodgent is a curated toolkit of reusable analysis functions implementing the
#' Routine Analysis Techniques for medical data science. It covers data
#' preprocessing, quantitative and qualitative
#' Descriptive statistics, conventional and advanced regression, survival
#' Analysis, competing-risk and multi-state models. All exported functions share
#' The `mc_` prefix and return tidy `data.frame` objects that drop straight into
#' Reports and downstream pipelines.
#'
#' @section Modules and functions:
#' \describe{
#' \item{Meet the data}{`mc_mode()`, `mc_inspect()`,
#' `mc_theme()`, `mc_palette()`}
#' \item{Preprocessing}{`mc_impute_knn()`,
#' `mc_impute_mice()`, `mc_dedup()`, `mc_scale()`, `mc_wide_long()`,
#' `mc_long_wide()`}
#' \item{Quantitative}{`mc_desc_quant()`,
#' `mc_norm_test()`, `mc_var_test()`, `mc_ttest()`, `mc_anova()`,
#' `mc_nonparam()`}
#' \item{Qualitative}{`mc_contingency()`, `mc_corr()`,
#' `mc_tableone()`}
#' \item{Conventional regression}{`mc_linear()`,
#' `mc_logistic()`, `mc_logistic_screen()`, `mc_stepwise()`,
#' `mc_ordinal_logistic()`, `mc_clogistic()`, `mc_model_compare()`}
#' \item{Survival}{`mc_survfit()`, `mc_cox()`,
#' `mc_ph_test()`}
#' \item{Advanced regression}{`mc_gam()`, `mc_lasso()`,
#' `mc_ridge()`, `mc_elastic_net()`, `mc_rcs()`, `mc_quantile_reg()`,
#' `mc_pcr()`, `mc_nnet_reg()`}
#' \item{Advanced survival}{`mc_competing_risk()`,
#' `mc_multistate()`}
#' }
#'
#' @section Design principles:
#' \enumerate{
#' \item All exported functions use the `mc_` prefix to avoid name clashes.
#' \item Analysis functions return tidy `data.frame` objects, ready for
#' Reporting and piping.
#' \item Heavy methods live in `Suggests` with graceful namespace guards.
#' }
#'
#' @docType package
#' @name macodgent-package
#' @aliases macodgent
#' @keywords internal
#' @seealso Entry points: \code{\link{mc_inspect}},
#' \code{\link{mc_desc_quant}}, \code{\link{mc_tableone}},
#' \code{\link{mc_logistic}}, \code{\link{mc_cox}}.
NULL
