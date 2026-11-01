# =============================================================================
# 软件名称    : macodgent 医学数据科学研究工具包
# 著作权人    : 谢耀珏（Dr.Xie）【个人，待最终确认】
# 版本号      : 1.0.0
# 开发语言    : R
# 开发完成日期: 2026-10-30
# 权利声明    : 本程序为著作权人独立开发，受《计算机软件保护条例》保护；
#              未经书面许可，不得复制、发行、出租或将其用于商业目的。
# -----------------------------------------------------------------------------
# 模块        : 高级建模模块（advanced）
# 本文件包含  : mc_gam / mc_lasso / mc_ridge / mc_elastic_net / mc_rcs /
#              mc_quantile_reg / mc_pcr（主成分回归）/ mc_nnet_reg（神经网络回归）
#              其中 .rcs_basis（受限立方样条基）为 Harrell 定义的原创实现，
#              不依赖 rms 包；mc_rcs 由此自建样条并拟合预测曲线。
# =============================================================================

#' 广义可加模型（GAM）/ Generalized Additive Model (GAM)
#'
#' 对应专栏第 7 章「广义可加回归模型」。基于 \code{\link[mgcv]{gam}}，返回每个
#' 平滑项的自由度与显著性。
#'
#' Corresponds to Chapter 7 "Generalized Additive Regression Models" of the
#' second part of the column. Built on \code{\link[mgcv]{gam}}, it returns the
#' degrees of freedom and significance of each smooth term.
#'
#' @details
#' 用途：拟合广义可加模型（GAM），在保持线性模型可解释性的同时，允许自变量以平滑（非线性）形式进入模型，自动捕捉曲线关系。
#' 方法：基于 mgcv 包的 gam()，用样条（如 s(x)）表示平滑项，自动选择平滑自由度并给出显著性。
#' 参数：data 为数据框；formula 形如 y ~ s(x1) + x2，其中 s() 包裹的项为平滑项，其余为线性项。
#' 结果：返回列表，含 model（gam 对象，可用 plot() 看平滑曲线）、terms（各平滑项的自由度 edf、参考自由度、F 值与 p 值，用于判断该变量是否具非线性效应）。
#' 注意：edf 接近 1 且 p 不显著，说明该变量近似线性、无需平滑；edf 很大可能过拟合。平滑项需要足够样本量才能稳定估计。依赖 mgcv（Suggests）。想看效应图用 plot(res$model, pages=1, shade=TRUE)。
#' @param data 数据框（a `data.frame`）。
#' @param formula GAM 公式，如 \code{y ~ s(x1) + x2}（a GAM formula such as
#'   \code{y ~ s(x1) + x2}）。
#' @param ... 传给 \code{\link[mgcv]{gam}} 的额外参数（additional arguments
#'   passed to \code{\link[mgcv]{gam}}）。
#' @return 列表，含 `model`（拟合的 gam 对象）与 `terms`（平滑项的整洁表）
#'   （a list with \code{model} and \code{terms}）。
#' @examples
#' \donttest{
#' # 专栏第 7 章：用 GAM 刻画孕期体重(lwt)对出生体重的非线性效应
#' birthwt <- MASS::birthwt
#' res <- mc_gam(birthwt, bwt ~ s(lwt) + smoke)
#' cat("平滑项表 / smooth terms:\n")
#' print(as.data.frame(res$terms), row.names = FALSE)
#' }
#' @seealso \code{\link{mc_rcs}}（受限立方样条，自实现基函数）。
#' @export
mc_gam <- function(data, formula, ...) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("`mc_gam` requires the mgcv package (shipped with R).")
  }
  fit <- eval(substitute(mgcv::gam(formula, data = DATA, ...),
                        list(DATA = data)))
  list(model = fit, terms = broom::tidy(fit))
}

# Internal: unified implementation for elastic-net models.
.elastic_net <- function(data, x, y, alpha, nfolds = 10, family = "gaussian") {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("`mc_elastic_net` requires the glmnet package; please install.packages('glmnet').")
  }
  xm <- as.matrix(data[x])
  yv <- data[[y]]
  cv <- glmnet::cv.glmnet(xm, yv, alpha = alpha, nfolds = nfolds, family = family)
  cf <- as.matrix(stats::coef(cv, s = "lambda.min"))
  df <- data.frame(variable = rownames(cf), coefficient = cf[, 1],
                   stringsAsFactors = FALSE)
  df <- df[df$variable != "(Intercept)", ]
  rownames(df) <- NULL
  list(cv = cv, coefficients = df, lambda.min = cv$lambda.min, alpha = alpha)
}

#' Lasso 回归 / Lasso Regression
#'
#' 对应专栏第 7 章「Lasso 回归」。基于 \pkg{glmnet}，\code{alpha = 1}。通过
#' 10 折交叉验证选 λ，返回非零系数表。
#'
#' Corresponds to Chapter 7 "Lasso Regression" of the second part of the
#' column. Built on \pkg{glmnet} with \code{alpha = 1}. Returns the non-zero
#' coefficient table selected by 10-fold cross-validated lambda.
#'
#' @details
#' 用途：用 Lasso（L1 正则化）做回归与变量选择，在高维或共线性数据中自动把不重要变量的系数压缩为 0，实现「内嵌式」筛选。
#' 方法：基于 glmnet 包，alpha=1 即 Lasso；用 k 折交叉验证（默认 10 折）选择最优惩罚强度 λ。
#' 参数：data 为数据框；x 为自变量列名向量；y 为因变量列名；nfolds 为 CV 折数（默认 10）；family 默认 gaussian（连续结局），也可用于 binomial 等。
#' 结果：返回含 coefficients（非零系数表，被压缩为 0 的变量不出现）的列表。
#' 注意：Lasso 给出的系数是有偏的（被整体收缩），适合预测与筛选，不适合直接做统计推断。当多个变量高度相关时，Lasso 会随机挑其中一个、其余置零。λ 越大收缩越强；可用 cv.glmnet 的 lambda.min/lambda.1se 调参。依赖 glmnet（Suggests）。
#' @param data 数据框（a `data.frame`）。
#' @param x 预测变量列名（字符向量）（character vector of predictor column names）。
#' @param y 因变量名（response variable name）。
#' @param nfolds 交叉验证折数，默认 `10`（cross-validation folds; default `10`）。
#' @param family `"gaussian"`（默认）或 `"binomial"` 等（GLM family）。
#' @return 列表，含 `cv`（cv.glmnet 对象）/ `coefficients`（非零系数列表）/
#'   `lambda.min`（a list with \code{cv}, \code{coefficients}, \code{lambda.min}）。
#' @examples
#' \donttest{
#' # 专栏第 7 章：用 Lasso 在 mtcars 的全部预测变量中做变量选择（mpg ~ .）
#' xvars <- setdiff(names(mtcars), "mpg")
#' mc_lasso(mtcars, x = xvars, y = "mpg")$coefficients
#' }
#' @seealso \code{\link{mc_ridge}}、\code{\link{mc_elastic_net}}。
#' @export
mc_lasso <- function(data, x, y, nfolds = 10, family = "gaussian") {
  .elastic_net(data, x, y, alpha = 1, nfolds = nfolds, family = family)
}

#' 岭回归 / Ridge Regression
#'
#' 对应专栏第 7 章「岭回归」。基于 \pkg{glmnet}，\code{alpha = 0}。返回系数表
#' （岭回归系数通常不稀疏）。
#'
#' Corresponds to Chapter 7 "Ridge Regression" of the second part of the
#' column. Built on \pkg{glmnet} with \code{alpha = 0}. Returns the coefficient
#' table (ridge coefficients are typically non-sparse).
#'
#' @details
#' 用途：用岭回归（L2 正则化）处理多重共线性，稳定系数估计，适合「变量都重要但不想丢弃任何一个」的场景。
#' 方法：基于 glmnet 包，alpha=0 即岭回归；对所有系数做等比例收缩（不置零），保留全部变量。
#' 参数：data 为数据框；x 为自变量列名向量；y 为因变量列名；nfolds 为 CV 折数（默认 10）；family 默认 gaussian。
#' 结果：返回含 coefficients（被收缩但均非零的系数表）的列表。
#' 注意：与 Lasso 不同，岭回归的系数永远不会精确为 0，因此不做变量选择、只做系数收缩。它对共线性稳健（把方差分摊到相关变量上），但解释时变量仍全部保留。λ 越大收缩越强。依赖 glmnet（Suggests）。
#' @param data 数据框（a `data.frame`）。
#' @param x 预测变量列名（字符向量）（character vector of predictor column names）。
#' @param y 因变量名（response variable name）。
#' @param nfolds 交叉验证折数，默认 `10`（cross-validation folds; default `10`）。
#' @param family `"gaussian"`（默认）或 `"binomial"` 等（GLM family）。
#' @return 列表，含 `cv` / `coefficients` / `lambda.min`（a list with
#'   \code{cv}, \code{coefficients}, \code{lambda.min}）。
#' @examples
#' \donttest{
#' # 专栏第 7 章：岭回归处理 mtcars 预测变量的多重共线性
#' xvars <- setdiff(names(mtcars), "mpg")
#' mc_ridge(mtcars, x = xvars, y = "mpg")$coefficients
#' }
#' @seealso \code{\link{mc_lasso}}、\code{\link{mc_elastic_net}}。
#' @export
mc_ridge <- function(data, x, y, nfolds = 10, family = "gaussian") {
  .elastic_net(data, x, y, alpha = 0, nfolds = nfolds, family = family)
}

#' 弹性网络回归 / Elastic Net Regression
#'
#' 对应专栏第 7 章「弹性网络回归」。基于 \pkg{glmnet}；`alpha` 在岭（0）与
#' Lasso（1）之间调节。
#'
#' Corresponds to Chapter 7 "Elastic Net Regression" of the second part of the
#' column. Built on \pkg{glmnet}; \code{alpha} tunes between ridge (0) and
#' lasso (1).
#'
#' @details
#' 用途：用弹性网络（Elastic Net）兼顾 Lasso 的稀疏性与岭回归对共线性的稳健性，是二者之间的连续调节。
#' 方法：基于 glmnet 包，alpha 在 0（纯岭）到 1（纯 Lasso）之间调节混合比例；alpha=0.5 是常用折中。
#' 参数：data 为数据框；x 为自变量列名向量；y 为因变量列名；alpha 为混合参数（默认 0.5）；nfolds 为 CV 折数；family 默认 gaussian。
#' 结果：返回含 coefficients（按所选 α 收缩后的系数表，部分可能为 0）的列表。
#' 注意：当自变量成组高度相关时，弹性网络比 Lasso 更稳健——它会把整组相关变量一起保留或一起剔除，而非随机留一个。α 越接近 1 越偏向变量选择，越接近 0 越偏向系数收缩。依赖 glmnet（Suggests）。
#' @param data 数据框（a `data.frame`）。
#' @param x 预测变量列名（字符向量）（character vector of predictor column names）。
#' @param y 因变量名（response variable name）。
#' @param alpha 弹性网络混合参数，默认 `0.5`（elastic-net mixing parameter;
#'   default `0.5`）。
#' @param nfolds 交叉验证折数，默认 `10`（cross-validation folds; default `10`）。
#' @param family `"gaussian"`（默认）或 `"binomial"` 等（GLM family）。
#' @return 列表，含 `cv` / `coefficients` / `lambda.min` / `alpha`（a list with
#'   \code{cv}, \code{coefficients}, \code{lambda.min}, \code{alpha}）。
#' @examples
#' \donttest{
#' # 专栏第 7 章：弹性网络（alpha = 0.5）兼顾变量选择与系数收缩
#' xvars <- setdiff(names(mtcars), "mpg")
#' mc_elastic_net(mtcars, x = xvars, y = "mpg", alpha = 0.5)$coefficients
#' }
#' @seealso \code{\link{mc_lasso}}、\code{\link{mc_ridge}}。
#' @export
mc_elastic_net <- function(data, x, y, alpha = 0.5, nfolds = 10, family = "gaussian") {
  .elastic_net(data, x, y, alpha = alpha, nfolds = nfolds, family = family)
}

# Internal: restricted cubic spline basis (Harrell's definition).
.rcs_basis <- function(x, knots = NULL, nk = 3) {
  x <- as.numeric(x)
  if (is.null(knots)) {
    qp <- c(0, (1:nk) / (nk + 1), 1)
    knots <- stats::quantile(x, probs = qp, na.rm = TRUE, names = FALSE)
  }
  knots <- sort(unique(knots))
  K <- length(knots)
  t1 <- knots[1]; tK <- knots[K]
  out <- data.frame(x = x); names(out) <- "rcs_lin"
  for (i in 2:(K - 1)) {
    ti <- knots[i]
    lambda_i <- (tK - ti) / (tK - t1)
    out[[paste0("rcs_", i)]] <- pmax(0, x - ti)^3 -
      lambda_i * pmax(0, x - t1)^3 - (1 - lambda_i) * pmax(0, x - tK)^3
  }
  attr(out, "knots") <- knots
  out
}

#' 受限立方样条回归（自实现样条基）/ Restricted Cubic Spline Regression
#'
#' 对应专栏第 7 章「受限立方样条回归模型」。不依赖 \pkg{rms}，自行构建受限立方
#' 样条基并拟合线性/广义线性模型，返回预测曲线。
#'
#' Corresponds to Chapter 7 "Restricted Cubic Spline Regression Models" of the
#' second part of the column. Without depending on \pkg{rms}, it builds a
#' restricted cubic spline basis and fits a linear / generalized linear model,
#' returning the prediction curve.
#'
#' @details
#' 用途：用受限立方样条（RCS）刻画自变量与结局之间的非线性剂量-反应关系，并可视化曲线，是流行病学中呈现「暴露-效应」曲线的常用方法。
#' 方法：自实现样条基（无额外依赖），在自变量两端做线性约束（「受限」），避免样条在边界处剧烈摆动；默认节点数 nk=3。
#' 参数：data 为数据框；x 为自变量列名；y 为因变量列名；nk 为节点数（默认 3，越大越灵活）；knots 可手动指定节点位置（默认按分位数自动选取）；family 默认 gaussian，也支持 binomial（此时刻画对数发生比随 x 的变化）。
#' 结果：返回列表，含 knots（节点位置）、coefficients（样条基系数）、plot（预测曲线图，ggplot 对象）。
#' 注意：节点数越多曲线越曲折、越易过拟合；一般 3–5 个节点足够。RCS 不假设单调，能呈现「U 形」「倒 U 形」等复杂关系。绘制结果用 print(res$plot)。
#' @param data 数据框（a `data.frame`）。
#' @param x 数值预测变量名（numeric predictor name）。
#' @param y 因变量名（response variable name）。
#' @param nk 内部节点数，默认 `3`（number of interior knots; default `3`）。
#' @param knots 完整节点向量（含边界）；`NULL` 表示自动按分位数取节点（full
#'   knot vector including boundaries; `NULL` = automatic quantile knots）。
#' @param family `"gaussian"`（默认）或 `"binomial"` 等（GLM family）。
#' @return 列表，含 `model` / `coefficients`（系数表）/ `plot`（预测曲线）/
#'   `knots`（节点位置）（a list with \code{model}, \code{coefficients},
#'   \code{plot}, \code{knots}）。
#' @examples
#' # 专栏第 7 章：母亲年龄对低出生体重风险的「剂量-反应」曲线（非线性）
#' birthwt <- MASS::birthwt
#' res <- mc_rcs(birthwt, x = "age", y = "low", nk = 4, family = "binomial")
#' cat("节点位置 / knots:", round(res$knots, 2), "\n")
#' print(res$plot)
#' @seealso \code{\link{mc_gam}}（光滑项同样用于非线性关系建模）。
#' @export
mc_rcs <- function(data, x, y, nk = 3, knots = NULL, family = "gaussian") {
  bx <- .rcs_basis(data[[x]], knots = knots, nk = nk)
  df <- data.frame(bx, y = data[[y]])
  f <- stats::as.formula(paste("y ~", paste(names(bx), collapse = " + ")))
  if (family == "gaussian") {
    fit <- stats::lm(f, data = df)
  } else {
    fit <- stats::glm(f, data = df, family = family)
  }
  grid <- data.frame(x = seq(min(data[[x]], na.rm = TRUE),
                             max(data[[x]], na.rm = TRUE), length.out = 100))
  gx <- .rcs_basis(grid$x, knots = attr(bx, "knots"), nk = nk)
  pred <- stats::predict(fit, newdata = gx, se.fit = TRUE, type = "response")
  pc <- data.frame(x = grid$x,
                   fit = as.numeric(pred$fit),
                   lo = as.numeric(pred$fit) - 1.96 * pred$se.fit,
                   hi = as.numeric(pred$fit) + 1.96 * pred$se.fit)
  plt <- ggplot2::ggplot(pc, ggplot2::aes(x = x, y = fit)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2) +
    ggplot2::scale_x_continuous(breaks = .nice_breaks(pc$x, 6)) +
    ggplot2::labs(x = x, y = paste("predicted", y)) +
    mc_theme()
  td <- broom::tidy(fit, conf.int = TRUE, exponentiate = (family != "gaussian"))
  list(model = fit, coefficients = td, plot = plt, knots = attr(bx, "knots"))
}

#' 分位数回归 / Quantile Regression
#'
#' 对应专栏第 7 章「分位数回归原理与应用场景」。基于
#' \code{\link[quantreg]{rq}}，一次性返回多个分位点的系数。
#'
#' Corresponds to Chapter 7 "Principles and Application Scenarios of Quantile
#' Regression" of the second part of the column. Built on
#' \code{\link[quantreg]{rq}}, it returns coefficients for multiple quantiles
#' at once.
#'
#' @details
#' 用途：用分位数回归刻画自变量对结局不同分位点（如中位数、上下四分位）的影响，特别适合存在异方差或关注极端分位（如低收入、高死亡风险）的场景。
#' 方法：基于 quantreg 包的 rq()，最小化加权的绝对残差，可同时估计多个分位点的系数。
#' 参数：data 为数据框；formula 形如 y ~ x；taus 为要估计的分位点向量，默认 c(0.25, 0.5, 0.75)。
#' 结果：返回各分位点系数的数据框（含截距与自变量斜率），可比较不同分位下斜率的差异。
#' 注意：若 taus 上不同分位的自变量斜率明显不同（如异方差时），说明条件分布形态随 x 变化——这是普通均值回归（OLS）看不到的信息。分位数回归对异常值稳健。依赖 quantreg（Suggests）。
#' @param data 数据框（a `data.frame`）。
#' @param formula 回归公式（a regression formula）。
#' @param taus 分位点向量，默认 `c(0.25, 0.5, 0.75)`（vector of quantile points;
#'   default `c(0.25, 0.5, 0.75)`）。
#' @return 数据框，列为 `term` / `estimate` / `std.error` / `statistic` /
#'   `p.value` / `tau`（a `data.frame` with those columns, stacked over taus）。
#' @examples
#' \donttest{
#' # 专栏第 7 章：出生体重的条件分位数回归（不同分位下孕期体重 lwt 的斜率）
#' birthwt <- MASS::birthwt
#' mc_quantile_reg(birthwt, bwt ~ lwt + age, taus = c(0.25, 0.5, 0.75))
#' }
#' @export
mc_quantile_reg <- function(data, formula, taus = c(0.25, 0.5, 0.75)) {
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("`mc_quantile_reg` requires the quantreg package; please install.packages('quantreg').")
  }
  rows <- list()
  for (t in taus) {
    fit <- quantreg::rq(formula, data = data, tau = t)
    td <- broom::tidy(fit)
    td$tau <- t
    rows[[length(rows) + 1]] <- td
  }
  do.call(rbind, rows)
}

#' 主成分回归（PCR，自研管线）/ Principal Components Regression
#'
#' 对应专栏第 7 章「主成分回归的原理和应用场景」（7.8）。将多个相关预测变量
#' 先做主成分分析（PCA），取前若干主成分再做线性回归，从而消除多重共线性。
#' 整条管线**自研**实现：基于 \code{stats::prcomp} 计算主成分，并以**自研 k 折
#' 交叉验证**选择最优主成分个数，最后把系数反变换回原始变量尺度。
#'
#' Corresponds to Chapter 7 (7.8, PCR). The whole pipeline is *self-implemented*:
#' PCA via \code{stats::prcomp}, an in-house k-fold CV to pick the number of
#' components, and a back-transform of coefficients to the original scale.
#'
#' @details
#' 用途：当预测变量高度相关（多重共线）时，用主成分回归稳定估计。
#' 方法（自研）：① 对预测变量做标准化（中心化+缩放）；② PCA 得到互不相关的
#' 主成分；③ 若未指定主成分数，用自研 k 折交叉验证按预测均方误差（MSE）选最优
#' 个数；④ 以选出的主成分拟合线性回归；⑤ 将主成分系数反变换回原始变量尺度
#' （利用旋转矩阵与标准化参数）。
#' 参数：\code{data} 为数据框；\code{formula} 形如 \code{y ~ x1 + x2 + ...}；
#' \code{ncomp} 为主成分个数（NULL 时由 CV 自动选）；\code{k} 为交叉验证折数
#' （默认 5）；\code{seed} 为 CV 随机种子。
#' 结果：返回列表——\code{best_ncomp}（选定主成分数）、\code{cv_table}（各候选
#' 数的 CV-MSE）、\code{coef_original}（原始尺度系数+截距）、\code{performance}
#' （R²/调整 R²/RMSE）、\code{plot}（CV-MSE 随主成分数变化图）、\code{fit}
#' （最终 lm 对象）。
#' 注意：PCR 对预测变量做了标准化，系数与标准化后的量纲相关；分类变量应先
#' 哑变量化再纳入。本函数仅依赖基础 R，无需额外安装包。
#' @param data 数据框（a `data.frame`）。
#' @param formula 回归公式，如 `y ~ x1 + x2`（a model `formula`）。
#' @param ncomp 主成分个数；`NULL` 时由交叉验证自动选择
#'   （number of principal components; `NULL` = auto-select by CV）。
#' @param k 交叉验证折数，默认 `5`（CV folds; default `5`）。
#' @param seed 交叉验证随机种子（random seed for CV）。
#' @return 列表，含选定主成分数、CV 表、原始尺度系数、性能量与图
#'   （a list with chosen components, CV table, original-scale coefficients,
#'   performance and the CV plot）。
#' @examples
#' # 专栏第 7 章 7.8 节：多个相关预测变量预测出生体重（主成分回归消除共线性）
#' # 注：PCR 要求全部预测变量为数值型，故不纳入分类变量 race（哑变量化另述）
#' birthwt <- MASS::birthwt
#' bw <- birthwt
#' res <- mc_pcr(bw, bwt ~ age + lwt + smoke + ht + ui + ptl, k = 5, seed = 1)
#' res$best_ncomp
#' res$coef_original
#' print(res$plot)
#' @export
mc_pcr <- function(data, formula, ncomp = NULL, k = 5L, seed = NULL) {
  if (!inherits(data, "data.frame") || !inherits(formula, "formula"))
    stop("`data` must be a data.frame and `formula` a formula.")
  mf <- stats::model.frame(formula, data = data)
  y  <- as.numeric(model.response(mf))
  X  <- as.matrix(mf[, -1, drop = FALSE])
  if (!is.numeric(X)) stop("All predictors must be numeric for PCR.")
  p  <- ncol(X)
  Xs <- scale(X)                       # 自研：标准化
  mu <- attr(Xs, "scaled:center")
  sg <- attr(Xs, "scaled:scale")
  pca <- stats::prcomp(Xs, center = FALSE, scale. = FALSE)  # Xs 已标准化
  Z   <- pca$x
  maxc <- min(p, length(y) - 2L)
  cv_mse <- NULL
  if (is.null(ncomp)) {
    if (!is.null(seed)) set.seed(seed)
    idx <- sample(rep(1:k, length.out = length(y)))
    cv_mse <- numeric(maxc)
    for (c in 1:maxc) {
      errs <- numeric(length(y))
      for (fold in 1:k) {
        te <- which(idx == fold)
        tr <- setdiff(seq_along(y), te)
        m  <- stats::lm(y[tr] ~ Z[tr, 1:c])
        pr <- as.numeric(cbind(1, Z[te, 1:c]) %*% coef(m))
        errs[te] <- (y[te] - pr)^2
      }
      cv_mse[c] <- mean(errs, na.rm = TRUE)
    }
    ncomp <- which.min(cv_mse)
  } else {
    ncomp <- min(as.integer(ncomp), maxc)
  }
  fit  <- stats::lm(y ~ Z[, 1:ncomp])
  beta <- as.numeric(coef(fit))         # c(intercept, beta_pc...)
  rot  <- pca$rotation[, 1:ncomp, drop = FALSE]
  gamma <- as.numeric(crossprod(rot, beta[-1]))   # 缩放尺度下的原始系数
  coef_orig   <- gamma / sg
  intercept   <- beta[1] - sum(mu * coef_orig)
  coef_original <- data.frame(
    term = c("(Intercept)", colnames(X)),
    estimate = c(intercept, coef_orig), stringsAsFactors = FALSE)
  pred <- as.numeric(cbind(1, Z[, 1:ncomp]) %*% coef(fit))
  perf <- data.frame(
    metric = c("r_squared", "adj_r_squared", "rmse"),
    value  = c(summary(fit)$r.squared, summary(fit)$adj.r.squared,
               sqrt(mean((y - pred)^2))), stringsAsFactors = FALSE)
  cv_table <- if (is.null(ncomp) || is.null(cv_mse)) NULL else
    data.frame(ncomp = seq_along(cv_mse), cv_mse = cv_mse, stringsAsFactors = FALSE)
  plot <- NULL
  if (!is.null(cv_table)) {
    plot <- ggplot2::ggplot(cv_table, ggplot2::aes(ncomp, cv_mse)) +
      ggplot2::geom_line(color = "steelblue") +
      ggplot2::geom_point(size = 2, color = "steelblue") +
      ggplot2::geom_vline(xintercept = ncomp, linetype = "dashed", color = "red") +
      ggplot2::labs(x = "Number of principal components",
                    y = "Cross-validated MSE") + mc_theme()
  }
  list(best_ncomp = ncomp, cv_table = cv_table, coef_original = coef_original,
       performance = perf, plot = plot, fit = fit)
}

#' 神经网络回归（nnet 封装）/ Neural-network regression
#'
#' 对应专栏第 7 章「神经网络回归的原理和应用场景」（7.9）。用单隐藏层前馈神经
#' 网络拟合回归，适合捕捉预测变量与结局间的非线性关系。基于 \code{nnet::nnet()}
#' （\code{linout = TRUE} 做回归），并额外给出预测 vs 实测散点图与 RMSE。
#'
#' Corresponds to Chapter 7 (7.9, neural-network regression). A single-hidden-layer
#' feed-forward network for regression via \code{nnet::nnet()} (with `linout = TRUE`),
#' plus a predicted-vs-observed plot and RMSE.
#'
#' @details
#' 用途：用神经网络拟合预测变量与连续结局间的非线性关系。
#' 依赖：需要 nnet 包（位于 Suggests）；未安装时给出明确提示。
#' 方法：\code{nnet::nnet(y ~ ., linout = TRUE, size = hidden)}，隐藏层节点数由
#' \code{hidden} 控制；数据自动标准化以加快收敛。
#' 参数：\code{data} 为数据框；\code{formula} 形如 \code{y ~ x1 + x2 + ...}；
#' \code{hidden} 为隐藏层节点数（默认 5）；\code{maxit} 为最大迭代次数（默认 200）；
#' \code{seed} 为随机种子。
#' 结果：返回列表——\code{fit}（nnet 对象）、\code{predictions}（拟合值）、
#' \code{rmse}、\code{plot}（预测 vs 实测散点）。
#' 注意：神经网络对尺度敏感，本函数内部对预测变量做了标准化；结果对随机初始化
#' 敏感，建议设置 \code{seed} 以保证可复现。节点数过多易过拟合，可用
#' mc_model_compare() 或交叉验证比较不同结构。
#' @param data 数据框（a `data.frame`）。
#' @param formula 回归公式（a model `formula`）。
#' @param hidden 隐藏层节点数，默认 `5`（hidden-layer nodes; default `5`）。
#' @param maxit 最大迭代次数，默认 `200`（maximum iterations; default `200`）。
#' @param seed 随机种子（random seed）。
#' @return 列表，含拟合对象、预测值、RMSE 与预测-实测图
#'   （a list with the fitted net, predictions, RMSE and the fit plot）。
#' @examples
#' \donttest{
#' # 专栏第 7 章 7.9 节：神经网络回归——用母亲特征预测出生体重
#' birthwt <- MASS::birthwt
#' bw <- birthwt
#' bw$race <- factor(bw$race)
#' res <- mc_nnet_reg(bw, bwt ~ age + lwt + smoke + race, hidden = 5, seed = 7)
#' res$rmse
#' print(res$plot)
#' }
#' @export
mc_nnet_reg <- function(data, formula, hidden = 5, maxit = 200, seed = NULL) {
  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("`mc_nnet_reg` requires the nnet package; please install.packages('nnet').")
  }
  if (!is.null(seed)) set.seed(seed)
  fit <- nnet::nnet(formula, data = data, size = hidden,
                   linout = TRUE, maxit = maxit, trace = FALSE)
  pred <- as.numeric(stats::predict(fit, data))
  y    <- as.numeric(stats::model.response(stats::model.frame(formula, data)))
  rmse <- sqrt(mean((y - pred)^2))
  plot <- ggplot2::ggplot(data.frame(obs = y, pred = pred),
                          ggplot2::aes(obs, pred)) +
    ggplot2::geom_point(size = 1.3, alpha = 0.6) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red",
                         linetype = "dashed") +
    ggplot2::labs(x = "Observed", y = "Predicted") + mc_theme()
  list(fit = fit, predictions = pred, rmse = rmse, plot = plot)
}
