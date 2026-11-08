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
# Module: utilities module (utils)
# This file contains: mc_stat_mode / mc_inspect / mc_theme / mc_palette
# (also includes the internal original algorithm.nice_breaks: "nice numbers"
# axis tick generator)
# =============================================================================

#' Mode (multiple modes supported)
#'
#' `mc_stat_mode()` is a teaching implementation of the mode, which base R does
#' Not provide. It returns *every* value tied for the highest frequency, so
#' Data with several modes is handled correctly instead of returning only one.
#'
#' @details
#' Purpose: find the value(s) with the highest frequency (the mode). Base R has
#' No built-in mode function, so this is a teaching-oriented reimplementation.
#' Arguments: x accepts a numeric or factor vector; na.rm controls whether
#' Missing values are removed before counting (default TRUE, so NAs are excluded
#' From the frequency tally).
#' Result: returns every value tied for the highest frequency, so multimodal
#' Data is handled naturally  -  a bimodal distribution returns both modes; the
#' Return type matches the input (numeric input yields a numeric vector).
#' Note: if every value occurs exactly once, all values are treated as modes and
#' Returned as-is, which is usually meaningless; run mc_inspect() first to
#' Understand the distribution. Returns NULL when all observations are missing.
#' @param x A numeric or factor vector.
#' @param na.rm Remove missing values first? Default `TRUE` (logical).
#' @return A vector of all mode(s), coerced back to the input type; `NULL` when
#' There are no non-missing values.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' mc_stat_mode(d$site)                 # most frequent simulated site
#' mc_stat_mode(as.numeric(d$race))
#' @export
mc_stat_mode <- function(x, na.rm = TRUE) {
 if (na.rm) x <- x[!is.na(x)]
 if (length(x) == 0) return(NULL)
 freq_table <- table(x)
 max_freq <- max(freq_table)
 modes <- names(freq_table)[freq_table == max_freq]
 # Restore the original type where possible (numeric case -> numeric)
 if (is.numeric(x)) modes <- as.numeric(modes)
 modes
}

#' Data frame overview
#'
#' One row per column: the column `class`, number of observations `n`,
#' Count of missing values `n_missing`, missing percentage `pct_missing`, and
#' Number of unique values `n_unique`. Handy for a quick data-quality check
#' Before modelling.
#'
#' @details
#' Purpose: a quick pre-modelling data-quality check; for each column it reports
#' The type, number of observations, number of missing values, missing
#' Percentage and number of unique values.
#' Arguments:.data is a data.frame or tibble; no other arguments.
#' Result: returns a data.frame with one row per input column. Columns: variable
#' (column name), class (type), n (valid observations), n_missing (missing
#' Count), pct_missing (missing percentage), n_unique (number of unique values).
#' Note: n is the count of valid (non-missing) observations, so n + n_missing
#' Equals the column's total length, and the missing percentage is computed
#' Against that total length. This is the first "health check" of the
#' Preprocessing pipeline  -  run it before any cleaning or modelling.
#' @param .data A `data.frame` or `tibble`.
#' @return A `data.frame` with those six columns, one row per input column.
#' @param detailed If `TRUE`, add unique-percentage, constant flag, outlier count, distribution verdict and a suggested variable role.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' mc_inspect(d[, c("id", "age", "bmi", "smoke", "crp")])
#' mc_inspect(d[, c("age", "bmi", "smoke", "crp")], detailed = TRUE)
#' @export
mc_inspect <- function(.data, detailed = FALSE) {
 if (!is.data.frame(.data)) stop("`.data` must be a data.frame or tibble.")
 rows <- lapply(names(.data), function(nm) {
 col <- .data[[nm]]
 n_miss <- sum(is.na(col))
 data.frame(
 variable = nm,
 class = paste(class(col), collapse = ","),
 n = length(col) - n_miss,
 n_missing = n_miss,
 pct_missing = round(100 * n_miss / length(col), 2),
 n_unique = length(unique(col[!is.na(col)])),
 stringsAsFactors = FALSE
 )
 })
 out <- do.call(rbind, rows)
 if (!isTRUE(detailed)) return(out)
 n <- nrow(.data)
 extra <- lapply(names(.data), function(nm) {
 col <- .data[[nm]]
 nun <- length(unique(col[!is.na(col)]))
 pct_unique <- round(100 * nun / n, 2)
 is_const <- nun <= 1L
 n_out <- NA_integer_
 distribution <- "categorical"
 if (is.numeric(col) && nun > 10L) {
 q <- stats::quantile(col, c(0.25, 0.75), na.rm = TRUE)
 iqr <- q[2] - q[1]
 n_out <- sum(col < q[1] - 1.5 * iqr | col > q[2] + 1.5 * iqr,
 na.rm = TRUE)
 sk <- tryCatch(mean((col - mean(col, na.rm = TRUE))^3, na.rm = TRUE) /
 stats::sd(col, na.rm = TRUE)^3,
 error = function(e) NA_real_)
 distribution <- if (!is.finite(sk)) "unknown" else
 if (sk > 0.5) "right-skewed" else if (sk < -0.5) "left-skewed" else
 "approximately symmetric"
 }
 role <- if (is_const) "constant (exclude)" else
 if (pct_unique >= 95 && nun > 10L &&
 grepl("(^id$|_id$|id$|patient|subject|record)", tolower(nm)))
 "identifier (exclude)" else
 if (pct_unique >= 95 && is.numeric(col)) "continuous predictor" else
 if (nun <= 10L) "categorical predictor" else "predictor"
 data.frame(pct_unique = pct_unique, is_constant = is_const,
 n_outliers = n_out, distribution = distribution,
 suggested_role = role, stringsAsFactors = FALSE)
 })
 out <- cbind(out, do.call(rbind, extra))
 rownames(out) <- NULL
 out
}

#' MA Codgent light ggplot2 theme
#'
#' A light, restrained ggplot2 theme with ample whitespace, ready for direct use
#' In reports. It is the default theme returned inside the `plot`/`forest`
#' Components of several analysis functions (e.g. `mc_kaplan_meier()`, `mc_cox()`,
#' `mc_rcs()`).
#'
#' @details
#' Purpose: a light, restrained ggplot2 theme with ample whitespace for report
#' Figures, matching the plot components returned by several analysis functions
#' (e.g. mc_kaplan_meier(), mc_cox(), mc_rcs()).
#' Arguments: base_size controls the base font size; default 12.
#' Result: returns a ggplot2 theme object that can be appended to any ggplot
#' With +.
#' Note: this is a pure theme (no geoms), so it must be combined with layers
#' Such as geom_point() or geom_line(); pair it with mc_palette() to change the
#' Colour scheme.
#' @param base_size Base font size; default `12`.
#' @return A `ggplot2` theme object, usable with `+`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' g <- ggplot2::ggplot(d[1:200, ], ggplot2::aes(sbp_baseline, sbp_week12)) +
#'   ggplot2::geom_point(alpha = 0.3) + mc_theme()
#' print(g)
#' @export
mc_theme <- function(base_size = 12) {
 ggplot2::theme_minimal(base_size = base_size) +
 ggplot2::theme(
 plot.title = ggplot2::element_text(size = base_size + 2, face = "bold"),
 panel.grid.minor = ggplot2::element_blank(),
 panel.grid.major = ggplot2::element_line(color = "#ECECEC"),
 axis.line = ggplot2::element_line(color = "#CCCCCC"),
 strip.background = ggplot2::element_rect(fill = "#F5F7FA", color = NA),
 legend.position = "bottom"
 )
}

#' Brand palette & palette generator
#'
#' Supplies a set of stable brand anchor colours (hex) and, instead of a static
#' Vector, **linearly interpolates** between anchors to generate any number of
#' Transition colours (a continuous gradient or a discrete ramp). Returns a named
#' Character vector usable with `ggplot2::scale_color_manual()` /
#' `scale_fill_manual()`.
#'
#' @details
#' Purpose: returns the brand anchor colours to keep plotting colours
#' Consistent; when the number of categories exceeds 6, it interpolates between
#' Anchors to generate the required number of colours, avoiding manual colour
#' Picking.
#' Method (original): internally it parses "#RRGGBB" into a 0-255 RGB triplet
#' (.hex_to_rgb), performs linear interpolation between anchor colours by
#' Position (.interp_palette), then converts the result back to hex
#' (.rgb_to_hex); the whole process uses no palette package.
#' Arguments: n is the number of colours to generate; NULL (default) returns the
#' Six fixed brand colours. type selects "brand" (default, six discrete brand
#' Colours), "gradient" (continuous gradient between brand colours, default 100
#' Colours), or "discrete" (a discrete ramp interpolated over the brand
#' Colours).
#' Result: by default returns a named character vector with the six brand
#' Colours: primary (blue), secondary (purple), accent (green), warn (red),
#' Info (cyan), muted (grey). When n is given, returns a colour vector of length
#' N (gradient is unnamed by default; discrete/brand keep anchor names).
#' Note: colours are generated purely algorithmically with no dependence on
#' Data; interpolation is in linear RGB space  -  switch to HCL space later if a
#' More perceptually uniform ramp is needed. Often paired with mc_theme() for
#' Style-consistent report figures.
#' @param n The number of colours to generate; `NULL` = the six fixed brand
#' Colours.
#' @param type Palette type: `"brand"`, `"gradient"` or `"discrete"`.
#' @return A named character vector of the six brand colours by default, or a
#' Vector of length `n`.
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' barplot(table(d$site), col = mc_palette(5), main = "Simulated site counts")
#' mc_palette(6, type = "gradient")
#' @export
mc_palette <- function(n = NULL, type = c("brand", "gradient", "discrete")) {
 brand <- c(primary = "#1A439C",
 secondary = "#511B78",
 accent = "#1C7331",
 warn = "#C0392B",
 info = "#0E7C86",
 muted = "#8A94A6")
 type <- match.arg(type)
 if (type == "brand" && is.null(n)) return(brand)
 if (type == "gradient") {
 if (is.null(n)) n <- 100L
 return(.interp_palette(c(brand["info"], brand["primary"], brand["secondary"]), n))
 }
 if (is.null(n)) n <- length(brand)
.interp_palette(brand, n)
}

# Internal (original): parse "#RRGGBB" into a 0-255 RGB triplet.
.hex_to_rgb <- function(hex) {
 hex <- sub("^#", "", hex)
 r <- strtoi(substr(hex, 1, 2), 16L)
 g <- strtoi(substr(hex, 3, 4), 16L)
 b <- strtoi(substr(hex, 5, 6), 16L)
 as.integer(c(r, g, b))
}

# Internal (original): convert a 0-255 RGB triplet back to "#RRGGBB".
.rgb_to_hex <- function(rgb) {
 rgb <- round(pmax(0, pmin(255, rgb)))
 sprintf("#%02X%02X%02X", rgb[1], rgb[2], rgb[3])
}

# Internal (original): linear interpolation between anchor hex colours,
# producing `n` transition colours (Heckbert-style segment stepping).
.interp_palette <- function(anchors, n) {
 if (n <= length(anchors)) {
 if (n == length(anchors)) return(anchors)
 return(anchors[round(seq(1, length(anchors), length.out = n))])
 }
 rgbs <- do.call(rbind, lapply(anchors,.hex_to_rgb))
 seg <- length(anchors) - 1L
 pos <- seq(0, 1, length.out = n)
 out <- character(n)
 for (i in seq_len(n)) {
 p <- pos[i] * seg
 k <- pmin(floor(p) + 1L, seg)
 t <- p - (k - 1L)
 a <- rgbs[k, ]; b <- rgbs[k + 1L, ]
 out[i] <-.rgb_to_hex(a + (b - a) * t)
 }
 out
}

# Internal (original): "nice numbers" axis tick generator (Heckbert 1990).
# Given a numeric vector and a desired tick count, returns a set of rounded,
# evenly spaced break values that fully cover the data range.
.nice_breaks <- function(x, n = 5) {
 x <- x[is.finite(x)]
 if (length(x) == 0) return(numeric(0))
 lo <- min(x); hi <- max(x)
 if (lo == hi) return(c(lo))
 rng <- hi - lo
 raw <- rng / max(1, (n - 1))
 mag <- 10 ^ floor(log10(raw))
 norm <- raw / mag
 nice <- if (norm < 1.5) 1 else if (norm < 3) 2 else if (norm < 7) 5 else 10
 step <- nice * mag
 start <- floor(lo / step) * step
 end <- ceiling(hi / step) * step
 seq(start, end, by = step)
}
