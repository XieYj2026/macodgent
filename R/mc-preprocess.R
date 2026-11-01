# =============================================================================
# 软件名称    : macodgent 医学数据科学研究工具包
# 著作权人    : Dr.Xie
# 版本号      : 1.0.0
# 开发语言    : R
# 开发完成日期: 2026-10-30
# 权利声明    : 本程序为著作权人独立开发，受《计算机软件保护条例》保护；
#              未经书面许可，不得复制、发行、出租或将其用于商业目的。
# -----------------------------------------------------------------------------
# 模块        : 预处理模块（preprocess）
# 本文件包含  : mc_impute_knn / mc_impute_mice / mc_dedup / mc_scale /
#              mc_wide_long / mc_long_wide
#              其中长宽转换 mc_wide_long / mc_long_wide 为原生 base R 实现，
#              不依赖外部“长宽转换”包。
# =============================================================================

#' K 近邻填补（自实现，无额外依赖）/ K-nearest-neighbour imputation
#'
#' 对应第 2 章「KNN 填补」。仅对数值列做填补：先按列标准化，再用欧氏距离
#' 找最近的 k 个非缺失观测，取其均值回填。非数值列的缺失值原样保留并给出警告。
#'
#' Corresponds to Chapter 2 ("KNN imputation"). Imputes only numeric columns
#' using standardised Euclidean distance; missing values in non-numeric columns
#' are left untouched with a warning. No extra dependency is required.
#'
#' @details
#' 用途：用 K 近邻法填补数值列的缺失值，适合缺失比例不高、且变量间存在相关性的场景。
#' 方法：只对数值列处理。先按列标准化（减均值除标准差），再用欧氏距离为每个缺失点寻找最近的 k 个非缺失观测，取它们均值回填；标准差为 0 的列退化为用列均值填补。非数值列的缺失值原样保留并给出警告。
#' 参数：.data 为数据框；k 为近邻个数，默认 5。
#' 结果：返回与输入结构一致的填补后数据框。
#' 注意：KNN 对量纲敏感，故内部已做标准化；但不同列取值范围差异极大时仍建议先用 mc_scale()。分类变量不参与 KNN，需另行处理（如众数填补或合并类别）。
#' @param .data 数据框（a `data.frame`）。
#' @param k 近邻个数，默认 `5`（number of neighbours; default `5`）。
#' @return 填补后的数据框（the imputed `data.frame`；结构与输入一致）。
#' @examples
#' # 第 2 章：以 iris 演示 KNN 缺失值填补
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

#' 多重填补（MICE 引擎 + 原创回填诊断）/ Multiple imputation (MICE + diagnostics)
#'
#' 对应第 2 章「链式方程多重填补（MICE）」。依赖 `mice` 包（位于 Suggests）。
#' 在 MICE 引擎之外，**额外计算逐变量的缺失回填诊断**，便于快速核查填补质量。
#'
#' Corresponds to Chapter 2 ("multiple imputation by chained equations").
#' Requires the \pkg{mice} package (Suggests). On top of the MICE engine, it
#' **adds an original per-variable imputation diagnostic**.
#'
#' @details
#' 用途：用链式方程多重填补（MICE）处理缺失，适合缺失机制较复杂、希望保留不确定性的场景，并在事后给出每列缺失回填的量化诊断。
#' 依赖：需要 mice 包（位于 Suggests），未安装时会报错提示。
#' 方法：对每个缺失变量，以其它与它相关的变量为预测因子循环建模，迭代生成多个完整数据集，从而保留缺失带来的不确定性。
#' 参数：.data 为数据框；m 为填补次数（生成 m 个完整数据集），默认 5；method 传 NULL 表示自动选择；... 透传给 mice::mice()。
#' 结果：返回列表，含 mids（原始多重填补对象，可做进一步分析）、long（由 mice::complete(action="long") 得到的长表），以及 summary——**原创诊断表**：逐变量给出 n_missing_original（原始缺失数）、n_imputed（被回填的观测数）、pct_imputed（回填比例）。
#' 注意：多重填补后的分析不应只看某一个数据集，而应使用 mids 对象配合支持 mice 的下游函数（如 with()/pool()）合并结果；summary 仅用于快速核查哪些变量被填补、填补了多少，不代表填补质量本身。
#' @param .data 数据框（a `data.frame`）。
#' @param m 填补次数（imputations），默认 `5`（number of imputed datasets;
#'   default `5`）。
#' @param method 传给 `mice()` 的方法；`NULL` 表示自动选择（method passed to
#'   \code{mice::mice()}; `NULL` = automatic）。
#' @param ... 其它传给 \code{mice::mice()} 的参数（other arguments passed to
#'   \code{mice::mice()}）。
#' @return 列表，含 `mids`（原始 \code{mids} 对象）、`long`（填补后的长表，
#'   由 \code{mice::complete(action = "long")} 得到），以及 `summary`（逐变量
#'   缺失回填诊断表）（a list with \code{mids}, \code{long}, and \code{summary}）。
#' @examples
#' \dontrun{
#' # 第 2 章：链式方程多重填补（MICE），并给出逐变量回填诊断
#' set.seed(106)
#' n <- 44
#' df <- data.frame(
#'   age  = c(rnorm(n - 4, 50), NA, NA, rnorm(2, 55)),
#'   bili = c(NA, rnorm(n - 4, 1.2), NA, rnorm(2, 1.5)),
#'   chol = rnorm(n, 200, 40), stringsAsFactors = FALSE)
#' res <- mc_impute_mice(df, m = 3)
#' head(res$long)
#' print(res$summary)   # 原创：逐变量缺失回填诊断
#' }
#' @export
mc_impute_mice <- function(.data, m = 5, method = NULL, ...) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("`mc_impute_mice` requires the mice package; please install.packages('mice').")
  }
  fit <- mice::mice(.data, m = m, method = method, printFlag = FALSE, ...)
  long <- mice::complete(fit, action = "long", include = FALSE)
  # ---- 原创诊断：逐变量缺失回填统计 ----
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

#' 去重 / Remove duplicates
#'
#' 对应第 2 章「去重」。默认按全部列去重；也可指定键列（\code{by}），
#' 保留每组的最后一行（或第一行）。
#'
#' Corresponds to Chapter 2 ("deduplication"). Deduplicates by all columns by
#' default, or by specified key columns, keeping the last row of each group.
#'
#' @details
#' 用途：去除重复行，常用于清洗重复录入或重复随访记录。
#' 参数：.data 为数据框；by 指定用于判断重复的键列（如受试者 ID），NULL 表示按所有列判断；keep 选择保留每组的 "last"（默认）或 "first" 行。
#' 结果：返回去重后的数据框。
#' 注意：当指定 by 时，只比较这些键列是否相同，其余列即使不同也会被判为重复并只保留一行——因此请确保键列能唯一标识一条记录，否则可能误删信息。这与 Excel 的「删除重复项」语义一致。
#' @param .data 数据框（a `data.frame`）。
#' @param by 键列名；`NULL` 表示按所有列（key column names; `NULL` = all
#'   columns）。
#' @param keep 保留 `"first"` 还是 `"last"`，默认 `"last"`（which row to keep;
#'   default `"last"`）。
#' @return 去重后的数据框（the deduplicated `data.frame`）。
#' @examples
#' # 第 2 章：去除重复录入的随访记录（按受试者 ID 去重）
#' df <- data.frame(subj = c(1, 1, 2, 2, 3), measure = c("a", "a", "b", "b", "c"), visit = 1:5)
#' mc_dedup(df, by = "subj")                  # 保留每组最后一行
#' mc_dedup(df, by = "subj", keep = "first")  # 保留每组第一行
#' @export
mc_dedup <- function(.data, by = NULL, keep = c("last", "first")) {
  keep <- match.arg(keep)
  if (!is.null(by)) {
    .data <- .data[!duplicated(.data[by], fromLast = (keep == "last")), ]
    return(.data)
  }
  .data[!duplicated(.data, fromLast = (keep == "last")), ]
}

#' 特征缩放 / Feature scaling
#'
#' 对应第 2 章「特征缩放」。支持三种方式：z-score 标准化（默认）、
#' min-max 归一化、以及基于中位数/IQR 的稳健缩放。仅替换数值列。
#'
#' Corresponds to Chapter 2 ("feature scaling"). Supports z-score (default),
#' min-max, and robust (median/IQR based) scaling. Only numeric columns are
#' replaced.
#'
#' @details
#' 用途：对数值列做特征缩放，常用于消除量纲差异、加速模型收敛，或为某些算法（如 KNN、正则化回归）准备数据。
#' 方法：支持三种缩放——z（默认，z-score：减均值除标准差，结果均值 0、标准差 1）、minmax（线性归一到 [0,1]）、robust（减中位数除 IQR，对异常值稳健）。
#' 参数：.data 为数据框；cols 指定要缩放的列，NULL 表示全部数值列；method 选择上述三种之一。
#' 结果：返回仅数值列被替换、其余列不变的数据框。
#' 注意：缩放会抹掉原始量纲，模型系数因此难以直接解释实际单位；若需解释性，建议保留原始变量或记录缩放参数以便还原。缺失值在缩放时按列剔除计算。
#' @param .data 数据框（a `data.frame`）。
#' @param cols 需要缩放的数值列；`NULL` 表示全部数值列（numeric columns to
#'   scale; `NULL` = all numeric columns）。
#' @param method 缩放方式：`"z"`（默认）/ `"minmax"` / `"robust"`（scaling
#'   method: `"z"`, `"minmax"` or `"robust"`）。
#' @return 缩放后的数据框（仅数值列被替换）（the scaled `data.frame`; only
#'   numeric columns are replaced）。
#' @examples
#' # 第 2 章：对 iris 数值列做特征缩放（消除量纲差异）
#' mc_scale(iris[, 1:4], method = "z")        # z 标准化
#' mc_scale(iris[, 1:4], method = "minmax")   # 归一到 [0,1]
#' mc_scale(iris[, 1:4], method = "robust")   # 稳健缩放
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

#' 宽表转长表（原生 base R 实现）/ Wide to long (native base R)
#'
#' 把指定测量列从宽表透视成长表，是数据处理的核心步骤。本函数**不依赖任何
#' 外部“长宽转换”包**，使用 base R 自行完成透视。
#'
#' Pivots the specified measure columns of a wide table into a long layout.
#' Implemented natively in base R with **no external pivoting dependency**.
#'
#' @details
#' 用途：把宽表（同一受试者的多次测量放在不同列）转换为长表（每次测量占一行），是数据处理的核心步骤。
#' 方法（原创，原生实现）：先校验 id_cols 与 measure_vars 都存在于数据框；随后对每个测量列，将标识列按测量列数复制堆叠，并把测量列名作为「变量名」、对应取值作为「取值」逐列合并，输出整齐的长表。全程仅用 base R，无外部依赖。
#' 参数：.data 为数据框；id_cols 为需保留的标识列（如 ID、分组）；measure_vars 为要拉长的测量列；name_to/value_to 分别为输出中「变量名」与「取值」列名，默认 variable / value。
#' 结果：返回长格式数据框。
#' 注意：转换后行数 = 原行数 × 测量列数；这是 mc_long_wide() 的逆操作。若测量列带有时间后缀（如 t1/t2/t3），转换后得到整齐的「时间-取值」两列，便于做重复测量分析。标识列为因子时会被保留为因子；测量列应类型一致，否则「取值」列会被强制转换为公共类型。
#' @param .data 数据框（a `data.frame`）。
#' @param id_cols 需要保留的标识列（identifier columns to keep）。
#' @param measure_vars 需要拉长的测量列（measured columns to pivot longer）。
#' @param name_to 输出「名称」列名，默认 `"variable"`（output name column;
#'   default `"variable"`）。
#' @param value_to 输出「取值」列名，默认 `"value"`（output value column;
#'   default `"value"`）。
#' @return 长格式数据框（the long `data.frame`）。
#' @examples
#' # 第 2 章：纵向血压测量的宽表转为长表（t1/t2/t3 三次访视）
#' wide <- data.frame(id = 1:4,
#'   t1 = rnorm(4, 120), t2 = rnorm(4, 122), t3 = rnorm(4, 118))
#' mc_wide_long(wide, id_cols = "id", measure_vars = c("t1", "t2", "t3"))
#' @export
mc_wide_long <- function(.data, id_cols, measure_vars,
                         name_to = "variable", value_to = "value") {
  if (!is.data.frame(.data)) stop("`.data` must be a data.frame.")
  miss <- setdiff(c(id_cols, measure_vars), names(.data))
  if (length(miss) > 0) {
    stop("以下列不存在于 .data 中：", paste(miss, collapse = ", "))
  }
  n <- nrow(.data)
  k <- length(measure_vars)
  id_block <- .data[rep(seq_len(n), each = k), id_cols, drop = FALSE]
  name_block <- rep(measure_vars, times = n)
  # 按 id 交错取值：第 i 个 id 连续取各测量列的第 i 行，再拼成 value 列
  value_block <- unlist(lapply(seq_len(n), function(i)
    sapply(measure_vars, function(mv) .data[[mv]][i])))
  out <- cbind(
    id_block,
    setNames(data.frame(name_block, stringsAsFactors = FALSE), name_to),
    setNames(data.frame(value_block, stringsAsFactors = FALSE), value_to)
  )
  rownames(out) <- NULL
  out
}

#' 长表转宽表（原生 base R 实现）/ Long to wide (native base R)
#'
#' 把长表按名称列展开回宽表，是 mc_wide_long() 的逆操作。本函数**不依赖任何
#' 外部“长宽转换”包**，使用 base R 自行完成展开。
#'
#' Spreads a long table back to wide. Implemented natively in base R with
#' **no external pivoting dependency**.
#'
#' @details
#' 用途：把长表按名称列展开回宽表，是 mc_wide_long() 的逆操作，常用于把「时间-取值」长表还原为受试者×时间点的宽矩阵。
#' 方法（原创，原生实现）：先取标识列的唯一组合作为输出行的骨架；再对名称列的每个取值，按「标识键 → 取值」的关系把对应 value_from 填回该新列。使用自定义的键拼接与匹配完成展开，全程仅依赖 base R。
#' 参数：.data 为长表；id_cols 为标识列；name_from 提供新列名的来源列；value_from 提供填充值的来源列。
#' 结果：返回宽格式数据框。
#' 注意：若标识+名称存在重复组合（同一 ID 同一名称出现多次），展开时按首次出现取值并给出警告，转换前应确保标识+名称唯一。这是 mc_wide_long() 的逆操作。
#' @param .data 数据框（a `data.frame`）。
#' @param id_cols 标识列（identifier columns）。
#' @param name_from 提供列名的来源列（column providing names）。
#' @param value_from 提供取值的来源列（column providing values）。
#' @return 宽格式数据框（the wide `data.frame`）。
#' @examples
#' \dontrun{
#' # 第 2 章：把「访视-血压」长表还原为受试者×访视的宽表
#' long_df <- data.frame(id = 1:2, variable = c("t1", "t2"),
#'                       value = c(120, 118))
#' mc_long_wide(long_df, id_cols = "id",
#'   name_from = "variable", value_from = "value")
#' }
#' @export
mc_long_wide <- function(.data, id_cols, name_from, value_from) {
  if (!is.data.frame(.data)) stop("`.data` must be a data.frame.")
  miss <- setdiff(c(id_cols, name_from, value_from), names(.data))
  if (length(miss) > 0) {
    stop("以下列不存在于 .data 中：", paste(miss, collapse = ", "))
  }
  id_block <- unique(.data[id_cols])
  key_cols <- function(d) do.call(paste, c(d[id_cols], sep = "\r"))
  out_keys <- key_cols(id_block)
  if (any(duplicated(.data[c(id_cols, name_from)]))) {
    warning("存在重复的（标识 + 名称）组合，展开时按首次出现取值。")
  }
  name_vals <- unique(.data[[name_from]])
  cols <- lapply(name_vals, function(nv) {
    sub <- .data[.data[[name_from]] == nv, , drop = FALSE]
    sub_keys <- key_cols(sub)
    vals <- sub[[value_from]][match(out_keys, sub_keys)]
    setNames(data.frame(vals, stringsAsFactors = FALSE), as.character(nv))
  })
  out <- cbind(id_block, do.call(cbind, cols))
  rownames(out) <- NULL
  out
}
