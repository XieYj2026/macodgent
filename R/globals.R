#' @importFrom stats sd median IQR setNames
NULL

# Declare variable names referenced via non-standard evaluation (NSE) to avoid
# spurious "undefined global variable" notes from R CMD check.
utils::globalVariables(c("sv", "g", "d", "z", "variable", "estimate",
                         "conf.low", "conf.high", "strata", "lo", "hi"))
