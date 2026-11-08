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
# This file contains: mc_impute_knn / mc_impute_mice / mc_dedup / mc_scale_numeric /
# mc_pivot_longer / mc_pivot_wider
# where the wide/long conversion functions mc_pivot_longer /
# mc_pivot_wider are native base R implementations with no
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
#' mc_scale_numeric() first. Categorical variables are not used in KNN and must be
#' Handled separately (e.g. mode imputation or category merging).
#' @param .data A \code{data.frame}.
#' @param k Number of neighbours; default `5`.
#' @return The imputed \code{data.frame} (same structure as the input).
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' x <- d[1:200, c("age", "bmi", "potassium", "quality_score")]
#' colSums(is.na(x))
#' imp <- mc_impute_knn(x, k = 5)
#' colSums(is.na(imp))
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
#' @param .data A \code{data.frame}.
#' @param m Number of imputations; default `5`.
#' @param method Method passed to \code{mice()}; \code{NULL} = automatic.
#' @param ... Other arguments passed to \code{mice::mice()}.
#' @return A list with \code{mids}, \code{long}, and \code{summary}: `mids` is
#' The raw \code{mids} object, `long` is the imputed long table (from
#' \code{mice::complete(action = "long")}), and `summary` is the per-variable
#' Imputation diagnostic table.
#' @param diagnostics If `TRUE` (default), assess chain-mean convergence and add a PASS / REVIEW check.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' \donttest{
#' x <- d[1:150, c("age", "bmi", "potassium", "quality_score")]
#' res <- mc_impute_mice(x, m = 3)
#' head(res$long)
#' print(res$summary)
#' print(res$convergence)
#' }
#' @export
mc_impute_mice <- function(.data, m = 5, method = NULL, ...,
                          diagnostics = TRUE) {
 if (!requireNamespace("mice", quietly = TRUE)) {
 stop("`mc_impute_mice` requires the mice package; please install.packages('mice').")
 }
 fit <- mice::mice(.data, m = m, method = method, printFlag = FALSE, ...)
 long <- mice::complete(fit, action = "long", include = FALSE)
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
 convergence <- NULL
 checks <- NULL
 if (isTRUE(diagnostics)) {
 cm <- tryCatch(fit$chainMean, error = function(e) NULL)
 if (!is.null(cm) && length(dim(cm)) >= 2L) {
 ## mice stores chainMean as variable x iteration x chain.
 vnames <- dimnames(cm)[[1]]
 if (is.null(vnames)) vnames <- names(.data)[seq_len(dim(cm)[1])]
 rows <- lapply(seq_along(vnames), function(j) {
 a <- tryCatch(cm[j, , , drop = FALSE], error = function(e) NULL)
 if (is.null(a)) return(NULL)
 mu <- apply(a, 2L, function(v) mean(as.numeric(v), na.rm = TRUE))
 mu <- mu[is.finite(mu)]
 if (length(mu) < 4L) return(NULL)
 half <- floor(length(mu) / 2L)
 drift <- abs(mean(mu[(half + 1L):length(mu)]) - mean(mu[seq_len(half)]))
 scale <- max(stats::sd(mu), 1e-8)
 data.frame(variable = vnames[j], n_iterations = length(mu),
 first_half_mean = mean(mu[seq_len(half)]),
 second_half_mean = mean(mu[(half + 1L):length(mu)]),
 relative_drift = drift / scale, stringsAsFactors = FALSE)
 })
 rows <- Filter(Negate(is.null), rows)
 convergence <- if (length(rows)) do.call(rbind, rows) else NULL
 if (!is.null(convergence) && nrow(convergence)) {
 convergence$status <- ifelse(convergence$relative_drift < 0.5,
 "PASS", "REVIEW")
 checks <- .mc_check_row(
 "MICE chain stability",
 suppressWarnings(max(convergence$relative_drift, na.rm = TRUE)),
 "max relative drift < 0.5",
 suppressWarnings(max(convergence$relative_drift, na.rm = TRUE)) < 0.5,
 "second-half vs first-half chain means")
 }
 }
 }
 list(mids = fit, long = long, summary = diag,
 convergence = convergence, checks = checks)
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
#' @param .data A \code{data.frame}.
#' @param by Key column names; \code{NULL} = all columns.
#' @param keep Whether to keep \code{"first"} or \code{"last"}; default \code{"last"}.
#' @return The deduplicated \code{data.frame}.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' x <- d[1:100, c("id", "age", "bmi", "treatment")]
#' dup <- rbind(x, x[1:5, ])            # simulated duplicate records
#' mc_dedup(dup, by = "id")
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
#' SD, yielding mean 0 and SD 1), minmax (linear transform to the 0-1 range), and
#' Robust (subtract median, divide by IQR; robust to outliers).
#' Arguments: `.data` is a data.frame; `cols` names the columns to scale, NULL
#' Means all numeric columns; `method` selects one of the three above.
#' Result: Returns a data.frame where only numeric columns are replaced; other
#' Columns are unchanged.
#' Note: Scaling discards the original units, so model coefficients are harder to
#' Interpret in real units; if interpretability is needed, keep the original
#' Variables or record the scaling parameters for back-transformation. Missing
#' Values are excluded column-wise when computing scaling.
#' @param .data A \code{data.frame}.
#' @param cols Numeric columns to scale; \code{NULL} = all numeric columns.
#' @param method Scaling method: \code{"z"} (default), \code{"minmax"}, or \code{"robust"}.
#' @return The scaled \code{data.frame}; only numeric columns are replaced.
#' @param train Optional training `data.frame` used to compute scaling parameters; prevents leakage into a test set.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' x <- d[1:200, c("age", "bmi", "sbp_baseline")]
#' mc_scale_numeric(x, method = "z")
#' mc_scale_numeric(x, cols = c("age", "bmi"), method = "minmax")
#' @export
mc_scale_numeric <- function(.data, cols = NULL,
 method = c("z", "minmax", "robust"),
 train = NULL) {
 method <- match.arg(method)
 ref <- if (is.null(train)) .data else train
 if (!is.data.frame(ref)) stop("`train` must be a data.frame when supplied.")
 if (is.null(cols)) cols <- names(.data)[sapply(.data, is.numeric)]
 out <- .data
 params <- list()
 for (col in cols) {
 x <- out[[col]]
 r <- ref[[col]]
 if (!is.numeric(x) || !is.numeric(r)) next
 if (method == "z") {
 mu <- mean(r, na.rm = TRUE); sc <- stats::sd(r, na.rm = TRUE)
 x <- (x - mu) / sc
 params[[col]] <- c(center = mu, scale = sc)
 } else if (method == "minmax") {
 rng <- range(r, na.rm = TRUE)
 x <- (x - rng[1]) / (rng[2] - rng[1])
 params[[col]] <- c(min = rng[1], max = rng[2])
 } else {
 med <- stats::median(r, na.rm = TRUE)
 iqr <- stats::IQR(r, na.rm = TRUE)
 x <- (x - med) / iqr
 params[[col]] <- c(median = med, iqr = iqr)
 }
 out[[col]] <- x
 }
 attr(out, "scaling") <- list(method = method, train_only = !is.null(train),
 parameters = params)
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
#' Is the inverse of mc_pivot_wider(). If measure columns carry time suffixes
#' (e.g. t1/t2/t3), the result is a tidy "time-value" pair of columns, convenient
#' For repeated-measures analysis. Identifier columns that are factors are
#' Preserved as factors; measure columns should share a common type, otherwise
#' The "value" column is coerced to a common type.
#' @param .data A \code{data.frame}.
#' @param id_cols Identifier columns to keep.
#' @param measure_vars Measured columns to pivot longer.
#' @param name_to Output name column; default \code{"variable"}.
#' @param value_to Output value column; default \code{"value"}.
#' @return The long \code{data.frame}.
#' @param check_duplicates Warn when `id_cols` do not uniquely identify rows; default `TRUE`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' w <- d[1:50, c("id", "sbp_baseline", "sbp_week4", "sbp_week8", "sbp_week12")]
#' long <- mc_pivot_longer(w, id_cols = "id",
#'   measure_vars = c("sbp_baseline", "sbp_week4", "sbp_week8", "sbp_week12"),
#'   name_to = "visit", value_to = "sbp")
#' head(long)
#' @export
mc_pivot_longer <- function(.data, id_cols, measure_vars,
 name_to = "variable", value_to = "value",
 check_duplicates = TRUE) {
 if (!is.data.frame(.data)) stop("`.data` must be a data.frame.")
 multi <- is.list(measure_vars) && !is.character(measure_vars)
 mv <- if (multi) unlist(measure_vars, use.names = FALSE) else measure_vars
 miss <- setdiff(c(id_cols, mv), names(.data))
 if (length(miss) > 0) {
 stop("The following columns do not exist in .data: ", paste(miss, collapse = ", "))
 }
 if (isTRUE(check_duplicates) && length(id_cols) &&
 anyDuplicated(.data[id_cols]) > 0L)
 warning("Identifier columns contain duplicate rows; the long output will repeat them.")
 n <- nrow(.data)
 k <- length(if (multi) measure_vars[[1L]] else measure_vars)
 id_block <- .data[rep(seq_len(n), each = k), id_cols, drop = FALSE]
 name_block <- rep(if (multi) measure_vars[[1L]] else measure_vars, times = n)
 values <- function(vars) unlist(lapply(seq_len(n), function(i)
 sapply(vars, function(mv2) .data[[mv2]][i])), use.names = FALSE)
 if (multi) {
 vnames <- names(measure_vars)
 if (is.null(vnames) || any(!nzchar(vnames)))
 vnames <- paste0(value_to[1L], seq_along(measure_vars))
 blocks <- lapply(measure_vars, values)
 val_df <- as.data.frame(setNames(blocks, vnames),
 stringsAsFactors = FALSE, check.names = FALSE)
 } else {
 val_df <- setNames(data.frame(values(measure_vars),
 stringsAsFactors = FALSE), value_to[1L])
 }
 out <- cbind(
 id_block,
 setNames(data.frame(name_block, stringsAsFactors = FALSE), name_to),
 val_df
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
#' mc_pivot_longer(), commonly used to restore a "time-value" long table to a
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
#' Inverse of mc_pivot_longer().
#' @param .data A \code{data.frame}.
#' @param id_cols Identifier columns.
#' @param name_from Column providing names.
#' @param value_from Column providing values.
#' @return The wide \code{data.frame}.
#' @param values_fill Value used to fill missing cells after spreading; default `NULL` (leave as NA).
#' @param values_fn Aggregation for duplicated identifier/name pairs: a function or `"mean"`, `"sum"`, `"median"`, `"max"`, `"min"`, `"first"`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' w <- d[1:50, c("id", "sbp_baseline", "sbp_week4", "sbp_week8", "sbp_week12")]
#' long <- mc_pivot_longer(w, id_cols = "id",
#'   measure_vars = c("sbp_baseline", "sbp_week4", "sbp_week8", "sbp_week12"),
#'   name_to = "visit", value_to = "sbp")
#' back <- mc_pivot_wider(long, id_cols = "id",
#'   name_from = "visit", value_from = "sbp")
#' head(back)
#' @export
mc_pivot_wider <- function(.data, id_cols, name_from, value_from,
 values_fill = NULL, values_fn = NULL) {
 if (!is.data.frame(.data)) stop("`.data` must be a data.frame.")
 miss <- setdiff(c(id_cols, name_from, value_from), names(.data))
 if (length(miss) > 0) {
 stop("The following columns do not exist in .data: ", paste(miss, collapse = ", "))
 }
 id_block <- unique(.data[id_cols])
 key_cols <- function(d) do.call(paste, c(d[id_cols], sep = "\r"))
 out_keys <- key_cols(id_block)
 dup <- anyDuplicated(.data[c(id_cols, name_from)]) > 0L
 fn <- NULL
 if (!is.null(values_fn)) {
 fn <- if (is.function(values_fn)) values_fn else {
 switch(as.character(values_fn)[1L],
 mean = function(v) mean(v, na.rm = TRUE),
 sum = function(v) sum(v, na.rm = TRUE),
 median = function(v) stats::median(v, na.rm = TRUE),
 max = function(v) max(v, na.rm = TRUE),
 min = function(v) min(v, na.rm = TRUE),
 first = function(v) v[1L],
 stop("`values_fn` must be a function or one of mean/sum/median/max/min/first."))
 }
 }
 if (dup && is.null(fn))
 warning("Duplicate (identifier + name) combinations exist; spreading takes the first occurrence.")
 name_vals <- unique(.data[[name_from]])
 cols <- lapply(name_vals, function(nv) {
 sub <- .data[.data[[name_from]] == nv, , drop = FALSE]
 if (!is.null(fn)) {
 kk <- key_cols(sub)
 vals <- vapply(out_keys, function(k) {
 v <- sub[[value_from]][kk == k]
 if (!length(v)) return(NA_real_)
 as.numeric(fn(v))[1L]
 }, numeric(1))
 } else {
 sub_keys <- key_cols(sub)
 vals <- sub[[value_from]][match(out_keys, sub_keys)]
 }
 if (!is.null(values_fill)) vals[is.na(vals)] <- values_fill
 setNames(data.frame(vals, stringsAsFactors = FALSE), as.character(nv))
 })
 out <- cbind(id_block, do.call(cbind, cols))
 rownames(out) <- NULL
 out
}
