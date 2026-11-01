# =============================================================================
# 软件名称    : macodgent 医学数据科学研究工具包
# 著作权人    : Dr.Xie
# 版本号      : 1.0.0
# 开发语言    : R
# 开发完成日期: 2026-10-30
# 权利声明    : 本程序为著作权人独立开发，受《计算机软件保护条例》保护；
#              未经书面许可，不得复制、发行、出租或将其用于商业目的。
# -----------------------------------------------------------------------------
# 模块        : 工具模块（utils）
# 本文件包含  : mc_mode / mc_inspect / mc_theme / mc_palette
#              （另含内部原创算法 .nice_breaks：坐标轴“美观刻度”生成）
# =============================================================================

#' 众数（支持多众数）/ Mode (multiple modes supported)
#'
#' 教学版众数函数：R 基础未提供。返回出现频率最高的全部取值，
#' 因此支持「多众数」场景（多个取值并列最高频）。
#'
#' `mc_mode()` is a teaching implementation of the mode, which base R does
#' not provide. It returns *every* value tied for the highest frequency, so
#' data with several modes is handled correctly instead of returning only one.
#'
#' @details
#' 用途：找出数据中出现频率最高的取值（众数）。R 基础语言没有内置众数函数，这里是面向教学的自实现版本。
#' 参数：x 接受数值或因子向量；na.rm 决定是否在计数前剔除缺失值（默认 TRUE，缺失值不计入频率）。
#' 结果：返回全部并列最高频的取值，因此天然支持「多众数」——双峰分布会同时返回两个众数；返回类型与输入保持一致（数值输入返回数值型向量）。
#' 注意：若每个取值都只出现一次，则所有值都会被判为「众数」并原样返回，这种情况通常没有实际意义；使用前可先用 mc_inspect() 了解分布形态。当所有观测都是缺失值时返回 NULL。
#' @param x 数值或因子向量（numeric or factor vector）。
#' @param na.rm 是否先剔除缺失值，默认 `TRUE`（logical; remove `NA`s before
#'   counting? Default `TRUE`）。
#' @return 众数向量（保持输入类型；数值输入返回数值型）；若无非缺失值则返回
#'   `NULL`（a vector of all mode(s), coerced back to the input type; `NULL` when
#'   there are no non-missing values）。
#' @examples
#' # 第 3 章：低出生体重研究数据集 birthwt（MASS 包）
#' birthwt <- MASS::birthwt
#' mc_mode(factor(birthwt$race))     # 种族变量的最常见类别（1 = 白人）
#' mc_mode(factor(birthwt$smoke))    # 吸烟状态众数（0 = 不吸烟）
#' mc_mode(c(NA, NA, 1, 2), na.rm = FALSE)   # 含 NA 时返回 "NA"
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

#' 数据框概览 / Data frame overview
#'
#' 逐列给出：类别（class）、观测数（n）、缺失数（n_missing）、缺失百分比
#' （pct_missing）与唯一值个数（n_unique）。常用于建模前快速把脉数据质量。
#'
#' One row per column: the column `class`, number of observations `n`,
#' count of missing values `n_missing`, missing percentage `pct_missing`, and
#' number of unique values `n_unique`. Handy for a quick data-quality check
#' before modelling.
#'
#' @details
#' 用途：建模前快速把脉数据质量，逐列给出类型、观测数、缺失数、缺失百分比与唯一值个数。
#' 参数：.data 为数据框或 tibble，无需其它参数。
#' 结果：返回一个数据框，每行对应输入的一列，列为 variable（列名）、class（类别）、n（有效观测数）、n_missing（缺失数）、pct_missing（缺失百分比）、n_unique（唯一值个数）。
#' 注意：n 是剔除缺失后的有效观测数，因此 n + n_missing 等于该列总长度，缺失百分比按列总长度计算。这是预处理流水线第一道「体检」步骤，建议在任何清洗/建模之前先跑一遍。
#' @param .data 数据框或 tibble（a `data.frame` or `tibble`）。
#' @return 数据框，列为 `variable` / `class` / `n` / `n_missing` /
#'   `pct_missing` / `n_unique`（a `data.frame` with those six columns, one row
#'   per input column）。
#' @examples
#' # 第 1 章：用 birthwt 体检数据框把脉数据质量（可先人为制造缺失演练）
#' birthwt <- MASS::birthwt
#' mc_inspect(birthwt[, c("age", "lwt", "bwt", "smoke", "low")])
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

#' MA Codgent 浅色 ggplot2 主题 / MA Codgent light ggplot2 theme
#'
#' 浅色、克制、留白充足的 ggplot2 主题，可直接用于报告出图，与各分析函数的
#' 返回值（如 `mc_survfit()`、`mc_cox()`、`mc_rcs()` 的 `plot` 组件）无缝衔接。
#'
#' A light, restrained ggplot2 theme with ample whitespace, ready for direct use
#' in reports. It is the default theme returned inside the `plot`/`forest`
#' components of several analysis functions (e.g. `mc_survfit()`, `mc_cox()`,
#' `mc_rcs()`).
#'
#' @details
#' 用途：提供一套浅色、克制、留白充足的 ggplot2 主题，用于报告出图，与多个分析函数返回值中的 plot 组件（如 mc_survfit()、mc_cox()、mc_rcs()）配色一致。
#' 参数：base_size 控制基础字号，默认 12。
#' 结果：返回一个 ggplot2 主题对象，可用 + 追加到任意 ggplot 图形之后。
#' 注意：这是纯主题（theme），不含几何对象（geom），必须配合 geom_point()、geom_line() 等图层使用；想换配色请配合 mc_palette() 一起使用。
#' @param base_size 基础字号，默认 `12`（base font size; default `12`）。
#' @return 一个 ggplot2 主题对象（a `ggplot2` theme object），可与 `+` 连用。
#' @examples
#' library(ggplot2)
#' # 第 3 章：母亲孕期体重与出生体重的散点（按是否吸烟着色）
#' birthwt <- MASS::birthwt
#' p <- ggplot(birthwt, aes(lwt, bwt, color = factor(smoke))) +
#'   geom_point(size = 2, alpha = 0.6) + mc_theme(base_size = 13)
#' print(p)
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

#' 品牌配色与调色板生成 / Brand palette & palette generator
#'
#' 提供一组稳定的品牌基准色（十六进制），并可**在基准色之间做线性插值**生成
#' 任意数量的过渡色（连续渐变或离散色板）。返回带名字的字符向量，可直接传给
#' ggplot2 的 `scale_color_manual()` / `scale_fill_manual()`。
#'
#' Supplies a set of stable brand anchor colours (hex) and, instead of a static
#' vector, **linearly interpolates** between anchors to generate any number of
#' transition colours (a continuous gradient or a discrete ramp). Returns a named
#' character vector usable with `ggplot2::scale_color_manual()` /
#' `scale_fill_manual()`.
#'
#' @details
#' 用途：返回品牌基准色并保证出图配色一致；当类别数超过 6 个时，可在基准色之间插值生成所需数量的颜色，避免手动配色。
#' 方法（原创）：内部先自行把 "#RRGGBB" 解析为 0–255 的 RGB 三元组（.hex_to_rgb），在锚点颜色间按位置做线性插值（.interp_palette），再把结果还原为十六进制（.rgb_to_hex）；整个过程不调用任何调色板包。
#' 参数：n 为需要生成的颜色个数，NULL（默认）时返回六个固定品牌色；type 选择 "brand"（默认，六个离散品牌色）、"gradient"（品牌色间的连续渐变，默认 100 色）、"discrete"（在品牌色基础上插值出的离散色板）。
#' 结果：默认返回一个带名字的字符向量，含 primary（主蓝）、secondary（紫）、accent（绿）、warn（红）、info（青）、muted（灰）六个品牌色；指定 n 时返回长度 n 的颜色向量（gradient 默认无名，discrete/brand 保留锚点名）。
#' 注意：这是纯算法生成的配色，不依赖任何数据；插值为线性 RGB 空间，如需感知更均匀的色板可后续改用 HCL 空间。常与 mc_theme() 搭配得到风格统一的报告图。
#' @param n 需要生成的颜色个数；`NULL`（默认）时返回六个固定品牌色
#'   （the number of colours to generate; `NULL` = the six fixed brand colours）。
#' @param type 调色板类型：`"brand"`（默认，品牌离散色）/ `"gradient"`
#'   （品牌色间的连续渐变）/ `"discrete"`（在品牌色基础上插值出的离散色）
#'   （palette type: `"brand"`, `"gradient"` or `"discrete"`）。
#' @return 命名字符向量（默认六个品牌色）或长度为 `n` 的颜色向量（a named
#'   character vector of the six brand colours by default, or a vector of length
#'   `n`）。
#' @examples
#' pal <- mc_palette()                 # 默认：六个品牌色
#' print(pal)
#' grad <- mc_palette(n = 20, type = "gradient")  # 20 色连续渐变
#' print(head(grad))
#' library(ggplot2)
#' df <- data.frame(name = factor(names(pal), names(pal)), col = pal)
#' p <- ggplot(df, aes(name, fill = name)) + geom_bar() +
#'   scale_fill_manual(values = pal) + theme_void()
#' print(p)
#' @export
mc_palette <- function(n = NULL, type = c("brand", "gradient", "discrete")) {
  brand <- c(primary   = "#1A439C",
             secondary = "#511B78",
             accent    = "#1C7331",
             warn      = "#C0392B",
             info      = "#0E7C86",
             muted     = "#8A94A6")
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
  rgbs <- do.call(rbind, lapply(anchors, .hex_to_rgb))
  seg <- length(anchors) - 1L
  pos <- seq(0, 1, length.out = n)
  out <- character(n)
  for (i in seq_len(n)) {
    p <- pos[i] * seg
    k <- pmin(floor(p) + 1L, seg)
    t <- p - (k - 1L)
    a <- rgbs[k, ]; b <- rgbs[k + 1L, ]
    out[i] <- .rgb_to_hex(a + (b - a) * t)
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
