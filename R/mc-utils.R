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
# This file contains: mc_mode / mc_inspect / mc_theme / mc_palette
# (also includes the internal original algorithm.nice_breaks: "nice numbers"
# axis tick generator)
# =============================================================================

#' Mode (multiple modes supported)
#'
#' `mc_mode()` is a teaching implementation of the mode, which base R does
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
#' # low birthweight study dataset (birthwt, MASS)
#' Birthwt <- MASS::birthwt
#' Mc_mode(factor(birthwt$race)) # most common race category (1 = white)
#' Mc_mode(factor(birthwt$smoke)) # mode of smoking status (0 = non-smoker)
#' Mc_mode(c(NA, NA, 1, 2), na.rm = FALSE) # returns "NA" when NA present
#' @export
mc_mode <- function(x, na.rm = TRUE) {
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
#' @param.data A `data.frame` or `tibble`.
#' @return A `data.frame` with those six columns, one row per input column.
#' @examples
#' # data-quality check on the birthwt data frame (try injecting
#' # missing values to practice)
#' Birthwt <- MASS::birthwt
#' Mc_inspect(birthwt[, c("age", "lwt", "bwt", "smoke", "low")])
#' @export
mc_inspect <- function(.data) {
 if (!is.data.frame(.data)) stop("`.data` must be a data.frame or tibble.")
 rows <- lapply(names(.data), function(nm) {
 col <-.data[[nm]]
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
 do.call(rbind, rows)
}

#' MA Codgent light ggplot2 theme
#'
#' A light, restrained ggplot2 theme with ample whitespace, ready for direct use
#' In reports. It is the default theme returned inside the `plot`/`forest`
#' Components of several analysis functions (e.g. `mc_survfit()`, `mc_cox()`,
#' `mc_rcs()`).
#'
#' @details
#' Purpose: a light, restrained ggplot2 theme with ample whitespace for report
#' Figures, matching the plot components returned by several analysis functions
#' (e.g. mc_survfit(), mc_cox(), mc_rcs()).
#' Arguments: base_size controls the base font size; default 12.
#' Result: returns a ggplot2 theme object that can be appended to any ggplot
#' With +.
#' Note: this is a pure theme (no geoms), so it must be combined with layers
#' Such as geom_point() or geom_line(); pair it with mc_palette() to change the
#' Colour scheme.
#' @param base_size Base font size; default `12`.
#' @return A `ggplot2` theme object, usable with `+`.
#' @examples
#' Library(ggplot2)
#' # maternal pregnancy weight vs. birth weight scatter
#' # (coloured by smoking status)
#' Birthwt <- MASS::birthwt
#' P <- ggplot(birthwt, aes(lwt, bwt, color = factor(smoke))) +
#' Geom_point(size = 2, alpha = 0.6) + mc_theme(base_size = 13)
#' Print(p)
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
#' Pal <- mc_palette() # default: six brand colours
#' Print(pal)
#' Grad <- mc_palette(n = 20, type = "gradient") # 20-colour continuous gradient
#' Print(head(grad))
#' Library(ggplot2)
#' Df <- data.frame(name = factor(names(pal), names(pal)), col = pal)
#' P <- ggplot(df, aes(name, fill = name)) + geom_bar() +
#' Scale_fill_manual(values = pal) + theme_void()
#' Print(p)
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
