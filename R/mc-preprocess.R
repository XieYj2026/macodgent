#' K-nearest-neighbour imputation (self-implemented, no extra dependency)
#'
#' Corresponds to Chapter 2 ("KNN imputation"). Imputes only numeric columns
#' using standardised Euclidean distance; missing values in non-numeric columns
#' are left untouched with a warning.
#'
#' @param .data A data.frame.
#' @param k Number of neighbours. Default 5.
#' @return The imputed data.frame.
#' @examples
#' df <- data.frame(a = c(1, 2, NA, 4, 5), b = c(2, NA, 3, 4, 5))
#' mc_impute_knn(df, k = 2)
#' @export
mc_impute_knn <- function(.data, k = 5) {
  if (!is.data.frame(.data)) stop("`.data` must be a data.frame.")
  num_cols <- names(.data)[sapply(.data, is.numeric)]
  if (length(num_cols) == 0) {
    warning("No numeric columns found; KNN imputation skipped.")
    return(.data)
  }
  out <- .data
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

#' Multiple imputation (MICE wrapper)
#'
#' Corresponds to Chapter 2 ("multiple imputation by chained equations").
#' Requires the \pkg{mice} package (Suggests). Returns the original \code{mids}
#' object and the combined long table.
#'
#' @param .data A data.frame.
#' @param m Number of imputations. Default 5.
#' @param method Method passed to mice; NULL means automatic.
#' @param ... Other arguments passed to \code{mice::mice()}.
#' @return A list with \code{mids} (the original object) and \code{long}
#'   (the imputed long table).
#' @examples
#' \dontrun{
#' mc_impute_mice(MASS::birthwt, m = 3)
#' }
#' @export
mc_impute_mice <- function(.data, m = 5, method = NULL, ...) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("`mc_impute_mice` requires the mice package; please install.packages('mice').")
  }
  fit <- mice::mice(.data, m = m, method = method, printFlag = FALSE, ...)
  long <- mice::complete(fit, action = "long", include = FALSE)
  list(mids = fit, long = long)
}

#' Remove duplicates
#'
#' Corresponds to Chapter 2 ("deduplication"). Deduplicates by all columns by
#' default, or by specified key columns, keeping the last row of each group.
#'
#' @param .data A data.frame.
#' @param by Key column names; NULL means all columns.
#' @param keep "first" or "last". Default "last".
#' @return The deduplicated data.frame.
#' @examples
#' df <- data.frame(id = c(1, 1, 2), v = c("a", "a", "b"))
#' mc_dedup(df, by = "id")
#' @export
mc_dedup <- function(.data, by = NULL, keep = c("last", "first")) {
  keep <- match.arg(keep)
  if (!is.null(by)) {
    .data <- .data[!duplicated(.data[by], fromLast = (keep == "last")), ]
    return(.data)
  }
  .data[!duplicated(.data, fromLast = (keep == "last")), ]
}

#' Feature scaling
#'
#' Corresponds to Chapter 2 ("feature scaling"). Supports z-score, min-max,
#' and robust (median/IQR based) scaling.
#'
#' @param .data A data.frame.
#' @param cols Numeric columns to scale; NULL means all numeric columns.
#' @param method "z" (default), "minmax" or "robust".
#' @return The scaled data.frame (only numeric columns are replaced).
#' @examples
#' mc_scale(mtcars, cols = c("mpg", "hp"), method = "z")
#' @export
mc_scale <- function(.data, cols = NULL, method = c("z", "minmax", "robust")) {
  method <- match.arg(method)
  if (is.null(cols)) cols <- names(.data)[sapply(.data, is.numeric)]
  out <- .data
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

#' Wide to long (tidyr wrapper)
#'
#' Corresponds to Chapter 2 ("wide/long conversion").
#'
#' @param .data A data.frame.
#' @param id_cols Identifier columns to keep.
#' @param measure_vars Measured columns to pivot longer.
#' @param name_to Output name column. Default "variable".
#' @param value_to Output value column. Default "value".
#' @return The long data.frame.
#' @examples
#' wide <- data.frame(id = 1:2, t1 = c(10, 20), t2 = c(15, 25))
#' mc_wide_long(wide, id_cols = "id", measure_vars = c("t1", "t2"))
#' @export
mc_wide_long <- function(.data, id_cols, measure_vars,
                         name_to = "variable", value_to = "value") {
  tidyr::pivot_longer(.data, cols = tidyselect_helpers(measure_vars, .data),
                      names_to = name_to, values_to = value_to)
}

#' Long to wide (tidyr wrapper)
#'
#' @param .data A data.frame.
#' @param id_cols Identifier columns.
#' @param name_from Column providing names.
#' @param value_from Column providing values.
#' @return The wide data.frame.
#' @examples
#' \dontrun{
#' mc_long_wide(long_df, id_cols = "id", name_from = "variable", value_from = "value")
#' }
#' @export
mc_long_wide <- function(.data, id_cols, name_from, value_from) {
  tidyr::pivot_wider(.data, id_cols = all_of_helpers(id_cols),
                     names_from = name_from, values_from = value_from)
}

# Internal: turn a character vector of column names into a tidyselect expression
tidyselect_helpers <- function(cols, .data) cols
all_of_helpers <- function(cols) cols
