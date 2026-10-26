#' Compute the mode (possibly multiple)
#'
#' Teaching implementation of a mode function, since R has no built-in.
#' Returns all values with the highest frequency (allowing multiple modes).
#'
#' @param x A numeric or factor vector.
#' @param na.rm Logical; remove missing values first? Default TRUE.
#' @return A vector of all modes (NULL if no non-missing values).
#' @examples
#' mc_mode(c(1, 2, 2, 3, 3, 3, 4, 5, 5, 6))
#' mc_mode(c("a", "a", "b", "b", "c"))   # multiple modes
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
#' One row per column: class, number of observations, number of missing
#' values, missing percentage, and number of unique values.
#'
#' @param .data A data.frame or tibble.
#' @return A data.frame with columns variable / class / n / n_missing /
#'   pct_missing / n_unique.
#' @examples
#' mc_inspect(MASS::birthwt)
#' @export
mc_inspect <- function(.data) {
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
  do.call(rbind, rows)
}

#' MA Codgent light ggplot2 theme
#'
#' Light, restrained, with ample whitespace; ready for direct use in reports.
#'
#' @param base_size Base font size. Default 12.
#' @return A ggplot2 theme object.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + mc_theme()
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

#' MA Codgent brand palette
#'
#' A stable set of categorical colours (hex codes) for consistent plotting.
#'
#' @return A named character vector.
#' @examples
#' mc_palette()
#' @export
mc_palette <- function() {
  c(
    primary   = "#1A439C",
    secondary = "#511B78",
    accent    = "#1C7331",
    warn      = "#C0392B",
    info      = "#0E7C86",
    muted     = "#8A94A6"
  )
}
