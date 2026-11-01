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
# Module: global variable declarations
# This file contains: variable declarations referenced via non-standard
# evaluation (NSE), to avoid spurious R CMD check notes
# =============================================================================

#' @importFrom stats sd median IQR setNames
NULL

# Declare variable names referenced via non-standard evaluation (NSE) to avoid
# spurious "undefined global variable" notes from R CMD check.
utils::globalVariables(c("sv", "g", "d", "z", "variable", "estimate",
 "conf.low", "conf.high", "strata", "lo", "hi"))
