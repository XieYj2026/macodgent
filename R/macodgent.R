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
#' The `mc_` prefix. Descriptive statistics and hypothesis tests return tidy
#' `data.frame` objects; model-fitting functions return a list holding the fitted
#' model, a coefficient table and, where applicable, a plot.
#'
#' @section Modules and functions:
#' \describe{
#' \item{Meet the data}{`mc_stat_mode()`, `mc_inspect()`,
#' `mc_theme()`, `mc_palette()`}
#' \item{Preprocessing}{`mc_impute_knn()`,
#' `mc_impute_mice()`, `mc_dedup()`, `mc_scale_numeric()`, `mc_pivot_longer()`,
#' `mc_pivot_wider()`}
#' \item{Quantitative}{`mc_describe_numeric()`,
#' `mc_test_normality()`, `mc_test_equal_variance()`, `mc_ttest()`, `mc_anova_oneway()`,
#' `mc_test_rank()`}
#' \item{Qualitative}{`mc_crosstab_test()`, `mc_correlate()`,
#' `mc_baseline_table()`}
#' \item{Conventional regression}{`mc_linear()`,
#' `mc_logistic()`, `mc_logistic_screen()`, `mc_stepwise()`,
#' `mc_ordinal_logistic()`, `mc_clogistic()`, `mc_model_compare()`}
#' \item{Survival}{`mc_kaplan_meier()`, `mc_cox()`,
#' `mc_test_proportional_hazards()`}
#' \item{Advanced regression}{`mc_gam()`, `mc_lasso()`,
#' `mc_ridge()`, `mc_elastic_net()`, `mc_rcs()`, `mc_quantile_reg()`,
#' `mc_pca_regression()`, `mc_nnet_reg()`}
#' \item{Advanced survival}{`mc_competing_risk()`,
#' `mc_multistate()`}
#' \item{Inference and reporting}{`mc_compare_groups()`,
#' `mc_adjust_p()`, `mc_posthoc_tukey()`, `mc_logistic_report()`}
#' }
#'
#' @section Inference layer:
#' Every hypothesis test reports an effect size with a confidence interval,
#' and every family of tests reports multiplicity-adjusted p values. Model
#' functions add explicit methodological checks (discrimination, calibration,
#' events per variable, collinearity, influence, proportional hazards).
#'
#' @section Design principles:
#' \enumerate{
#' \item All exported functions use the `mc_` prefix to avoid name clashes.
#' \item Descriptive statistics and hypothesis tests return tidy `data.frame`
#' objects, ready for reporting and piping. Model-fitting functions return a
#' list with the fitted model, a coefficient table and an optional plot.
#' \item Heavy methods live in `Suggests` with graceful namespace guards.
#' }
#'
#' @docType package
#' @name macodgent-package
#' @aliases macodgent
#' @keywords internal
#' @seealso Entry points: \code{\link{mc_inspect}},
#' \code{\link{mc_describe_numeric}}, \code{\link{mc_baseline_table}},
#' \code{\link{mc_logistic}}, \code{\link{mc_cox}}.
NULL
