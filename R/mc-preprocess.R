# =============================================================================
# Software name : macodgent medical data-science research toolkit
# Copyright owner : Dr.Xie
# Version : 1.0.0
# Language : R
# Completion date : 2026-10-30
# Copyright notice : This program is independently developed by the copyright
# owner and protected under the Software Copyright Protection
# Regulation; without written permission, it may not be
# copied, distributed, rented, or used for commercial purposes.
# -----------------------------------------------------------------------------
# Module : Preprocessing module
# This file contains: mc_impute_knn / mc_impute_mice / mc_dedup / mc_scale /
# mc_wide_long / mc_long_wide
# where the wide/long conversion functions mc_wide_long /
# mc_long_wide are native base R implementations with no
# external pivoting package dependency.
# =============================================================================

#' K-nearest-neighbour imputation
#'
#' Imputes only numeric columns
#' Using standardised Euclidean distance; missing values in non-numeric columns
#' Are left untouched with a warning. No extra dependency is required.
#'
#' @details
#' Purpose: Impute missing values in numeric columns using K-nearest neighbours,
#' Suitable when missingness is low and variables are correlated.
#' Method: Only numeric columns are processed. Columns are standardised
#' (mean-centred, divided by SD), then for each missing point the k nearest
#' Non-missing observations are found by Euclidean distance and their mean is
#' Used to impute; columns with zero SD fall back to column-mean imputation.
#' Missing values in non-numeric columns are left untouched with a warning.
#' Arguments: `.data` is a data.frame; `k` is the number of neighbours, default 5.
#' Result: Returns an imputed data.frame with the same structure as the input.
#' Note: KNN is sensitive to scale, so standardisation is done internally;
#' However, if columns differ greatly in range, it is still recommended to run
#' mc_scale() first. Categorical variables are not used in KNN and must be
#' Handled separately (e.g. mode imputation or category merging).
#' @param.data A \code{data.frame}.
#' @param k Number of neighbours; default `5`.
#' @return The imputed \code{data.frame} (same structure as the input).
#' @examples
#' # demonstrating kNN missing-value imputation with iris
#' data(iris)
#' set.seed(105)
#' iris_na <- iris
#' iris_na$Sepal.Length[sample(150, 15)] <- NA
#' iris_na$Sepal.Width[sample(150, 15)] <- NA
#' head(mc_impute_knn(iris_na, k = 3))
#' @export
mc_impute_knn <- function(.data, k = 5) {
 if (!is.data.frame(.data)) stop("`.data` must be a data.frame.")
 num_cols <- names(.data)[sapply(.data, is.numeric)]
 if (length(num_cols) == 0) {
 warning("No numeric columns found; KNN imputation skipped.")
 return(.data)
 }
 out <-.data
 for (col in num_cols) {
 x <- out[[col]]
 miss <- is.na(x)
 if (!any(miss)) next
 present <- !miss
 # Standardise for distance computation
 sdx <- sd(x[present])
 if (sdx == 0 || is.na(sdx)) {
 x[miss] <- mean(x[present])
 out[[col]] <- x
 next
 }
 z <- (x - mean(x[present])) / sdx
 ref_idx <- which(present)
 for (i in which(miss)) {
 d <- abs(z[i] - z[ref_idx])
 ord <- order(d)[seq_len(min(k, length(ref_idx)))]
 x[i] <- mean(x[ref_idx[ord]])
 }
 out[[col]] <- x
 }
 non_num_miss <- sum(sapply(out, function(c) sum(is.na(c)))) -
 sum(sapply(out[num_cols], function(c) sum(is.na(c))))
 if (non_num_miss > 0) {
 warning("Missing values in non-numeric columns were not handled by KNN; please process categorical variables separately.")
 }
 out
}

#' Multiple imputation (MICE + diagnostics)
#'
#'
#' Requires the \pkg{mice} package (Suggests). On top of the MICE engine, it
#' **adds an original per-variable imputation diagnostic**.
#'
#' @details
#' Purpose: Use multiple imputation by chained equations (MICE) to handle
#' Missingness, suitable when the missingness mechanism is complex and you want
#' To retain uncertainty; it also provides a quantitative per-column imputation
#' Diagnostic afterwards.
#' Dependency: Requires the mice package (Suggests); an informative error is
#' Raised if it is not installed.
#' Method: For each variable with missingness, related variables are used as
#' Predictors in an iterative model that generates multiple complete datasets,
#' Preserving the uncertainty due to missingness.
#' Arguments: `.data` is a data.frame; `m` is the number of imputations
#' (datasets generated), default 5; `method` of NULL means automatic selection;
#' `...` is passed to mice::mice().
#' Result: Returns a list containing `mids` (the raw multiple-imputation object
#' For further analysis), `long` (the long table from
#' Mice::complete(action = "long")), and `summary`  -  an **original diagnostic
#' Table**: per variable it gives `n_missing_original` (original missing count),
#' `n_imputed` (count of imputed observations), and `pct_imputed` (imputation
#' Proportion).
#' Note: Analysis after multiple imputation should not rely on a single dataset;
#' Instead use the `mids` object with downstream functions supporting mice
#' (e.g. with()/pool()) to pool results. The `summary` is only for quickly
#' Checking which variables were imputed and by how much, and does not represent
#' Imputation quality itself.
#' @param.data A \code{data.frame}.
#' @param m Number of imputations; default `5`.
#' @param method Method passed to \code{mice()}; \code{NULL} = automatic.
#' @param... Other arguments passed to \code{mice::mice()}.
#' @return A list with \code{mids}, \code{long}, and \code{summary}: `mids` is
#' The raw \code{mids} object, `long` is the imputed long table (from
#' \code{mice::complete(action = "long")}), and `summary` is the per-variable
#' Imputation diagnostic table.
#' @examples
#' \dontrun{
#' # multiple imputation by chained equations (mICE) with per-variable diagnostic
#' set.seed(106)
#' n <- 44
#' df <- data.frame(
#' age = c(rnorm(n - 4, 50), NA, NA, rnorm(2, 55)),
#' bili = c(NA, rnorm(n - 4, 1.2), NA, rnorm(2, 1.5)),
#' chol = rnorm(n, 200, 40), stringsAsFactors = FALSE)
#' res <- mc_impute_mice(df, m = 3)
#' head(res$long)
#' print(res$summary) # original: per-variable imputation diagnostic
#' }
#' @export
mc_impute_mice <- function(.data, m = 5, method = NULL,...) {
 if (!requireNamespace("mice", quietly = TRUE)) {
 stop("`mc_impute_mice` requires the mice package; please install.packages('mice').")
 }
 fit <- mice::mice(.data, m = m, method = method, printFlag = FALSE,...)
 long <- mice::complete(fit, action = "long", include = FALSE)
 # ---- original diagnostic: per-variable imputation statistics ----
 orig_na <- colSums(is.na(.data))
 imp_n <- if (length(fit$imp) > 0) {
 v <- sapply(fit$imp, nrow)
 out <- rep(0L, ncol(.data))
 names(out) <- names(.data)
 out[names(v)] <- v
 out
 } else rep(0L, ncol(.data))
 diag <- data.frame(
 variable = names(.data),
 n_missing_original = as.integer(orig_na),
 n_imputed = as.integer(imp_n),
 pct_imputed = round(100 * imp_n / pmax(1, nrow(.data)), 2),
 stringsAsFactors = FALSE
 )
 list(mids = fit, long = long, summary = diag)
}

#' Remove duplicates
#'
#' Deduplicates by all columns by
#' Default, or by specified key columns, keeping the last row of each group.
#'
#' @details
#' Purpose: Remove duplicate rows, commonly used to clean duplicate entries or
#' Repeated follow-up records.
#' Arguments: `.data` is a data.frame; `by` names the key columns used to detect
#' Duplicates (e.g. subject ID), NULL means judge by all columns; `keep` chooses
#' Whether to retain the "last" (default) or "first" row of each group.
#' Result: Returns the deduplicated data.frame.
#' Note: When `by` is specified, only those key columns are compared; other
#' Columns are ignored even if they differ, and duplicates are collapsed to a
#' Single row. Ensure the key columns uniquely identify a record, or information
#' May be lost. This matches Excel's "Remove Duplicates" semantics.
#' @param.data A \code{data.frame}.
#' @param by Key column names; \code{NULL} = all columns.
#' @param keep Whether to keep \code{"first"} or \code{"last"}; default \code{"last"}.
#' @return The deduplicated \code{data.frame}.
#' @examples
#' # removing duplicate follow-up records (deduplicate by subject iD)
#' df <- data.frame(subj = c(1, 1, 2, 2, 3), measure = c("a", "a", "b", "b", "c"), visit = 1:5)
#' mc_dedup(df, by = "subj") # keep the last row of each group
#' mc_dedup(df, by = "subj", keep = "first") # keep the first row of each group
#' @export
mc_dedup <- function(.data, by = NULL, keep = c("last", "first")) {
 keep <- match.arg(keep)
 if (!is.null(by)) {
.data <-.data[!duplicated(.data[by], fromLast = (keep == "last")), ]
 return(.data)
 }
.data[!duplicated(.data, fromLast = (keep == "last")), ]
}

#' Feature scaling
#'
#' Supports z-score (default),
#' Min-max, and robust (median/IQR based) scaling. Only numeric columns are
#' Replaced.
#'
#' @details
#' Purpose: Scale numeric columns to remove scale differences, speed up model
#' Convergence, or prepare data for algorithms such as KNN or regularised
#' Regression.
#' Method: Three scaling options  -  z (default; z-score: subtract mean, divide by
#' SD, yielding mean 0 and SD 1), minmax (linear transform to [0, 1]), and
#' Robust (subtract median, divide by IQR; robust to outliers).
#' Arguments: `.data` is a data.frame; `cols` names the columns to scale, NULL
#' Means all numeric columns; `method` selects one of the three above.
#' Result: Returns a data.frame where only numeric columns are replaced; other
#' Columns are unchanged.
#' Note: Scaling discards the original units, so model coefficients are harder to
#' Interpret in real units; if interpretability is needed, keep the original
#' Variables or record the scaling parameters for back-transformation. Missing
#' Values are excluded column-wise when computing scaling.
#' @param.data A \code{data.frame}.
#' @param cols Numeric columns to scale; \code{NULL} = all numeric columns.
#' @param method Scaling method: \code{"z"} (default), \code{"minmax"}, or \code{"robust"}.
#' @return The scaled \code{data.frame}; only numeric columns are replaced.
#' @examples
#' # scaling numeric columns of iris (removing scale differences)
#' mc_scale(iris[, 1:4], method = "z") # z-score standardisation
#' mc_scale(iris[, 1:4], method = "minmax") # normalise to [0, 1]
#' mc_scale(iris[, 1:4], method = "robust") # robust scaling
#' @export
mc_scale <- function(.data, cols = NULL, method = c("z", "minmax", "robust")) {
 method <- match.arg(method)
 if (is.null(cols)) cols <- names(.data)[sapply(.data, is.numeric)]
 out <-.data
 for (col in cols) {
 x <- out[[col]]
 if (!is.numeric(x)) next
 if (method == "z") {
 x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
 } else if (method == "minmax") {
 rng <- range(x, na.rm = TRUE)
 x <- (x - rng[1]) / (rng[2] - rng[1])
 } else {
 med <- median(x, na.rm = TRUE)
 iqr <- IQR(x, na.rm = TRUE)
 x <- (x - med) / iqr
 }
 out[[col]] <- x
 }
 out
}

#' Wide to long (native base R)
#'
#' Pivots the specified measure columns of a wide table into a long layout.
#' Implemented natively in base R with **no external pivoting dependency**.
#'
#' @details
#' Purpose: Convert a wide table (multiple measurements of one subject in
#' Different columns) to a long table (one row per measurement); a core step in
#' Data processing.
#' Method (original, native implementation): First verify that `id_cols` and
#' `measure_vars` exist in the data.frame; then for each measure column, replicate
#' And stack the identifier columns by the number of measure columns, and merge
#' Each measure column name as the "variable" and its values as the "value",
#' Producing a tidy long table. Uses only base R with no external dependency.
#' Arguments: `.data` is a data.frame; `id_cols` are the identifier columns to
#' Keep (e.g. ID, group); `measure_vars` are the measure columns to lengthen;
#' `name_to`/`value_to` are the output column names for "variable" and "value",
#' Defaulting to variable / value.
#' Result: Returns a long-format data.frame.
#' Note: After conversion, rows = original rows x number of measure columns; this
#' Is the inverse of mc_long_wide(). If measure columns carry time suffixes
#' (e.g. t1/t2/t3), the result is a tidy "time-value" pair of columns, convenient
#' For repeated-measures analysis. Identifier columns that are factors are
#' Preserved as factors; measure columns should share a common type, otherwise
#' The "value" column is coerced to a common type.
#' @param.data A \code{data.frame}.
#' @param id_cols Identifier columns to keep.
#' @param measure_vars Measured columns to pivot longer.
#' @param name_to Output name column; default \code{"variable"}.
#' @param value_to Output value column; default \code{"value"}.
#' @return The long \code{data.frame}.
#' @examples
#' # convert wide table of longitudinal blood-pressure measurements to long (t1/t2/t3 visits)
#' wide <- data.frame(id = 1:4,
#' t1 = rnorm(4, 120), t2 = rnorm(4, 122), t3 = rnorm(4, 118))
#' mc_wide_long(wide, id_cols = "id", measure_vars = c("t1", "t2", "t3"))
#' @export
mc_wide_long <- function(.data, id_cols, measure_vars,
 name_to = "variable", value_to = "value") {
 if (!is.data.frame(.data)) stop("`.data` must be a data.frame.")
 miss <- setdiff(c(id_cols, measure_vars), names(.data))
 if (length(miss) > 0) {
 stop("The following columns do not exist in.data: ", paste(miss, collapse = ", "))
 }
 n <- nrow(.data)
 k <- length(measure_vars)
 id_block <-.data[rep(seq_len(n), each = k), id_cols, drop = FALSE]
 name_block <- rep(measure_vars, times = n)
 # Interleave values by id: the i-th id takes row i of each measure column, then assemble into the value column
 value_block <- unlist(lapply(seq_len(n), function(i)
 sapply(measure_vars, function(mv).data[[mv]][i])))
 out <- cbind(
 id_block,
 setNames(data.frame(name_block, stringsAsFactors = FALSE), name_to),
 setNames(data.frame(value_block, stringsAsFactors = FALSE), value_to)
 )
 rownames(out) <- NULL
 out
}

#' Long to wide (native base R)
#'
#' Spreads a long table back to wide. Implemented natively in base R with
#' **no external pivoting dependency**.
#'
#' @details
#' Purpose: Spread a long table back to wide by the name column; the inverse of
#' mc_wide_long(), commonly used to restore a "time-value" long table to a
#' Subject x time-point wide matrix.
#' Method (original, native implementation): First take the unique combinations
#' Of identifier columns as the skeleton of output rows; then for each value of
#' The name column, fill the corresponding `value_from` back into the new column
#' According to the "identifier key -> value" relationship. Spreading is done via
#' Custom key concatenation and matching, relying only on base R.
#' Arguments: `.data` is the long table; `id_cols` are the identifier columns;
#' `name_from` is the source column supplying new column names; `value_from` is
#' The source column supplying fill values.
#' Result: Returns a wide-format data.frame.
#' Note: If duplicate (identifier + name) combinations exist (same ID and name
#' Appearing multiple times), the spread takes the first occurrence and issues a
#' Warning; ensure identifier + name is unique before converting. This is the
#' Inverse of mc_wide_long().
#' @param.data A \code{data.frame}.
#' @param id_cols Identifier columns.
#' @param name_from Column providing names.
#' @param value_from Column providing values.
#' @return The wide \code{data.frame}.
#' @examples
#' \dontrun{
#' # restore a "visit-blood pressure" long table to a subject x visit wide table
#' long_df <- data.frame(id = 1:2, variable = c("t1", "t2"),
#' value = c(120, 118))
#' mc_long_wide(long_df, id_cols = "id",
#' name_from = "variable", value_from = "value")
#' }
#' @export
mc_long_wide <- function(.data, id_cols, name_from, value_from) {
 if (!is.data.frame(.data)) stop("`.data` must be a data.frame.")
 miss <- setdiff(c(id_cols, name_from, value_from), names(.data))
 if (length(miss) > 0) {
 stop("The following columns do not exist in.data: ", paste(miss, collapse = ", "))
 }
 id_block <- unique(.data[id_cols])
 key_cols <- function(d) do.call(paste, c(d[id_cols], sep = "\r"))
 out_keys <- key_cols(id_block)
 if (any(duplicated(.data[c(id_cols, name_from)]))) {
 warning("Duplicate (identifier + name) combinations exist; spreading takes the first occurrence.")
 }
 name_vals <- unique(.data[[name_from]])
 cols <- lapply(name_vals, function(nv) {
 sub <-.data[.data[[name_from]] == nv, , drop = FALSE]
 sub_keys <- key_cols(sub)
 vals <- sub[[value_from]][match(out_keys, sub_keys)]
 setNames(data.frame(vals, stringsAsFactors = FALSE), as.character(nv))
 })
 out <- cbind(id_block, do.call(cbind, cols))
 rownames(out) <- NULL
 out
}
