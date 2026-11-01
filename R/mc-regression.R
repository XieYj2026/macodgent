# =============================================================================
# 软件名称    : macodgent 医学数据科学研究工具包
# 著作权人    : Dr.Xie
# 版本号      : 1.0.0
# 开发语言    : R
# 开发完成日期: 2026-10-30
# 权利声明    : 本程序为著作权人独立开发，受《计算机软件保护条例》保护；
#              未经书面许可，不得复制、发行、出租或将其用于商业目的。
# -----------------------------------------------------------------------------
# 模块        : 回归建模模块（regression）
# 本文件包含  : mc_logistic / mc_logistic_screen / mc_stepwise / mc_ordinal_logistic
#              mc_clogistic（条件Logistic/匹配病例对照）/ mc_model_compare（自研模型比较）
#              mc_logistic_screen（先单因素筛选、再多因素建模）为原创组合流程。
# =============================================================================

#' 二分类 logistic 回归 / Binary logistic regression
#'
#' 对应第 5 章：基于 \code{glm(family = binomial)}，返回整洁的系数表
#' （含 OR 及其 95% 置信区间，已实现指数化）。
#'
#' Corresponds to Chapter 5: based on \code{glm(family = binomial)}, returns a
#' tidy coefficient table (with OR and its 95% CI, already exponentiated).
#'
#' @details
#' 用途：拟合二分类 logistic 回归，得到各预测因子对结局发生的对数优势比及 OR。
#' 方法：基于最大似然估计的二项 logit 模型，因子/字符变量会自动处理为分类。
#' 参数：data 为数据框；formula 形如 结局 ~ 自变量（结局须为 0/1 或二水平因子）；conf 为 OR 的置信水平，默认 0.95。
#' 结果：返回整洁系数表，含估计值（logit 尺度）、标准误、z 值、p 值，以及已指数化的 OR 与 OR 的置信区间，便于直接解读。
#' 注意：OR>1 表示因子每增加一单位（或该水平相对参照）结局发生比升高；置信区间跨 1 提示该因子在 α=0.05 下不显著。建模前应确认无严重多重共线与空单元。结局非二分类时请改用多分类/有序模型。
#' @param data 数据框（a `data.frame`）。
#' @param formula 模型公式，如 \code{y ~ x1 + x2}（a model formula）。
#' @param conf 置信水平，默认 `0.95`（confidence level; default `0.95`）。
#' @return 数据框，列为 `term` / `estimate`(OR) / `std.error` / `statistic` /
#'   `p.value` / `conf.low` / `conf.high`（a tidy `data.frame`; `estimate` is the
#'   exponentiated OR）。
#' @examples
#' # 第 5 章 5.7 节：低出生体重的多因素 logistic 回归
#' birthwt <- MASS::birthwt
#' mc_logistic(birthwt, low ~ age + lwt + smoke + factor(race))
#' @export
mc_logistic <- function(data, formula, conf = 0.95) {
  fit <- eval(substitute(stats::glm(formula, data = DATA,
                                    family = stats::binomial),
                         list(DATA = data)))
  td <- broom::tidy(fit, conf.int = TRUE, conf.level = conf, exponentiate = TRUE)
  td
}

#' 单因素筛选 → 多因素 logistic 回归（旗舰）/ Univariate screen -> multivariable
#'
#' 对应第 5 章（拓展阅读）：先对每个候选变量做单因素 logistic 回归，按
#' `alpha` 筛出显著变量，再用这些变量构建多因素模型，一步完成「先单后多」。
#'
#' Corresponds to Chapter 5 (extended reading): run a univariate logistic
#' regression for each candidate variable, screen significant ones at \code{alpha},
#' then build a multivariable model from the selected variables.
#'
#' @details
#' 用途：一键实现「先单因素筛选、再多因素建模」的经典流行病学套路，省去手动分步的麻烦。
#' 方法：先对每个候选 x_vars 单独跑单因素 logistic 回归，按 alpha 阈值（默认 0.05）筛选显著变量，再把入选变量一起放进多因素模型。
#' 参数：data 为数据框；y 为二分类结局列名（字符串）；x_vars 为候选自变量列名向量；alpha 为单因素筛选的显著性阈值；conf 为置信水平。
#' 结果：返回列表，含 univariate（各变量单因素结果）、selected（入选变量名）、multivariable（多因素整洁系数表）。
#' 注意：此方法方便但不是金标准——逐步/筛选可能引入选择偏倚与过拟合，且单因素显著≠多因素显著（变量间可能相互掩盖或放大）。报告时应同时呈现单因素与多因素结果，并说明筛选规则。最终模型仍建议用 mc_logistic() 复核。
#' @param data 数据框（a `data.frame`）。
#' @param y 二分类结局变量名（字符）（binary outcome variable name, character）。
#' @param x_vars 候选预测变量名（字符向量）（candidate predictor names,
#'   character vector）。
#' @param alpha 单因素筛选显著性阈值，默认 `0.05`（univariate screening
#'   significance threshold; default `0.05`）。
#' @param conf 多因素模型置信水平，默认 `0.95`（multivariable confidence level;
#'   default `0.95`）。
#' @return 列表，含 `univariate`（单因素表）/ `selected`（入选变量）/ 
#'   `multivariable`（多因素系数表，无入选变量时为 `NULL`）（a list with
#'   \code{univariate}, \code{selected}, and \code{multivariable}）。
#' @examples
#' # 第 5 章（拓展阅读）：先单因素筛选、再多因素建模（"先单后多"）
#' birthwt <- MASS::birthwt
#' bw <- birthwt
#' bw$race  <- factor(bw$race)
#' bw$smoke <- factor(bw$smoke)
#' res <- mc_logistic_screen(bw, y = "low",
#'                           x_vars = c("age","lwt","race","smoke","ht","ui"))
#' cat("==== 单因素 / univariate ====\n"); print(res$univariate, row.names = FALSE)
#' cat("\n==== 入选变量 / selected ====\n"); print(res$selected)
#' cat("\n==== 多因素 / multivariable ====\n"); print(as.data.frame(res$multivariable), row.names = FALSE)
#' @export
mc_logistic_screen <- function(data, y, x_vars, alpha = 0.05, conf = 0.95) {
  univ <- list()
  sig <- c()
  for (xv in x_vars) {
    f <- stats::as.formula(paste(y, "~", xv))
    m <- eval(substitute(stats::glm(f, data = DATA, family = stats::binomial),
                        list(DATA = data)))
    td <- broom::tidy(m, conf.int = TRUE, exponentiate = TRUE)
    rows_x <- td[td$term != "(Intercept)", ]
    if (nrow(rows_x) == 0) next
    p_min <- min(rows_x$p.value)
    best <- rows_x[which.min(rows_x$p.value), ]
    univ[[xv]] <- data.frame(
      variable = xv, OR = best$estimate,
      CI.low = best$conf.low, CI.high = best$conf.high,
      p.value = p_min, stringsAsFactors = FALSE)
    if (p_min < alpha) sig <- c(sig, xv)
  }
  univ_df <- do.call(rbind, univ)
  if (length(sig) > 0) {
    f2 <- stats::as.formula(paste(y, "~", paste(sig, collapse = " + ")))
    m2 <- eval(substitute(stats::glm(f2, data = DATA, family = stats::binomial),
                          list(DATA = data)))
    mult <- broom::tidy(m2, conf.int = TRUE, conf.level = conf, exponentiate = TRUE)
    mult$variable <- mult$term
  } else {
    mult <- NULL
    message("Univariate screening selected no variables; multivariable model not built.")
  }
  list(univariate = univ_df, selected = sig, multivariable = mult)
}

#' 逐步变量选择 / Stepwise variable selection
#'
#' 对应第 5、7 章：基于 AIC 的逐步选择（默认 `k = 2`），返回系数表
#' （OR / 回归系数）。
#'
#' Corresponds to Chapters 5 and 7: stepwise selection based on AIC
#' (default \code{k = 2}), returning a coefficient table (OR / regression
#' coefficients).
#'
#' @details
#' 用途：基于信息准则（AIC）做逐步变量选择，自动从全模型中增删变量，得到更简洁的回归模型。
#' 方法：以 AIC 最小化为目标，支持 both（默认，双向，可进可出）、backward（从全模型删）、forward（从空模型加）；k 为惩罚项（默认 2 对应 AIC，取 log(n) 可得 BIC）。
#' 参数：data 为数据框；formula 为含全部候选自变量的公式；direction 选择方向；k 为 AIC/BIC 惩罚。
#' 结果：返回列表，含最终 coefficients（系数表）与 formula（最终入选变量的公式）。
#' 注意：逐步回归便于探索，但统计上存在缺陷——p 值被低估、置信区间偏窄、易过拟合，且结果依赖进入/剔除的随意阈值。它适合做变量筛选的「初稿」，重要结论仍需用领域知识与预先注册的分析计划支撑。
#' @param data 数据框（a `data.frame`）。
#' @param formula 完整模型公式（full model formula）。
#' @param direction `"both"`（默认）/ `"backward"` / `"forward"`（selection
#'   direction）。
#' @param k AIC 惩罚项，默认 `2`（AIC penalty; default `2`）。
#' @return 列表，含 `model`（最终 glm 对象）/ `coefficients`（系数表）/
#'   `formula`（最终公式）（a list with \code{model}, \code{coefficients},
#'   \code{formula}）。
#' @examples
#' # 第 5/7 章：基于 AIC 的逐步变量选择
#' birthwt <- MASS::birthwt
#' bw <- birthwt
#' bw$race  <- factor(bw$race)
#' bw$smoke <- factor(bw$smoke)
#' mc_stepwise(bw, low ~ age + lwt + race + smoke + ht + ui)$coefficients
#' @export
mc_stepwise <- function(data, formula, direction = c("both", "backward", "forward"),
                        k = 2) {
  direction <- match.arg(direction)
  fit <- eval(substitute(stats::glm(formula, data = DATA,
                                   family = stats::binomial),
                         list(DATA = data)))
  step_fit <- stats::step(fit, direction = direction, k = k, trace = 0)
  td <- broom::tidy(step_fit, conf.int = TRUE, exponentiate = TRUE)
  td$variable <- td$term
  list(model = step_fit, coefficients = td, formula = formula(step_fit))
}

#' 有序 logistic 回归 / Ordinal logistic regression
#'
#' 对应第 5 章「有序 logistic 回归」。基于 \pkg{MASS} 的 \code{polr}，
#' 结局须为有序因子。返回系数表（含 OR）。
#'
#' Corresponds to Chapter 5 ("ordinal logistic regression"). Based on
#' \pkg{MASS}'s \code{polr}; the outcome must be an ordered factor. Returns a
#' coefficient table with OR.
#'
#' @details
#' 用途：拟合有序 logistic 回归（比例优势模型），适用于结局为有序多分类（如轻/中/重度）的场景。
#' 方法：基于 MASS 包的 polr()，假设各分割点的自变量效应相同（比例优势假设）。
#' 参数：data 为数据框；formula 形如 有序结局 ~ 自变量，结局必须是有序因子（ordered()）。
#' 结果：返回含系数（logit 尺度）、标准误、与已指数化的 OR（优势比）的系数表。
#' 注意：若比例优势假设不成立（可用 MASS::polr 配合检验或分层的 OR 差异大），结果解释需谨慎，可改用部分比例优势模型或无序多分类 logistic。务必先把结局转成 ordered 因子，否则模型无意义。依赖 MASS（Suggests）。
#' @param data 数据框（a `data.frame`）。
#' @param formula 模型公式；结局为有序因子（the outcome is an ordered factor）。
#' @return 数据框，列为 `variable` / `Value` / `Std.Error` / `t.value` / `OR` /
#'   `p.value`（a `data.frame` of coefficients; `OR = exp(Value)`）。
#' @examples
#' \donttest{
#' # 第 5 章：有序 logistic——产前检查次数(ftv)的 ordinal 结局
#' birthwt <- MASS::birthwt
#' bw <- birthwt
#' bw$ftv_ord <- ordered(cut(bw$ftv, c(-1, 0, 2, Inf),
#'                          labels = c("0次", "1-2次", "3次+")))
#' bw$smoke <- factor(bw$smoke)
#' mc_ordinal_logistic(bw, ftv_ord ~ age + lwt + smoke)
#' }
#' @export
mc_ordinal_logistic <- function(data, formula) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("`mc_ordinal_logistic` requires the MASS package (shipped with R).")
  }
  fit <- eval(substitute(MASS::polr(formula, data = DATA, Hess = TRUE),
                        list(DATA = data)))
  co <- summary(fit)$coefficients
  cf <- as.data.frame(co)
  cf$variable <- rownames(cf)
  cf$OR <- exp(cf$Value)
  cf$p.value <- 2 * stats::pnorm(-abs(cf$`t value`))
  cf
}

#' 条件 Logistic 回归（匹配病例对照）/ Conditional logistic regression
#'
#' 对应第 5 章「条件 Logistic 回归模型的分析」（5.10）。用于 1:M 配对病例
#' 对照研究，通过 \code{survival::clogit()} 在每层（配对组）内做条件似然估计，
#' 消除配对因素混杂。相比普通 Logistic，它不估计配对截距，只估计暴露的
#' 条件优势比（OR）。
#'
#' Corresponds to Chapter 5 (5.10, conditional logistic regression). Suitable for
#' 1:M matched case-control studies; fits \code{survival::clogit()} which
#' conditions on each matched set, yielding a conditional odds ratio.
#'
#' @details
#' 用途：分析配对（匹配）病例对照数据中暴露与结局的关联，控制配对变量带来的混杂。
#' 依赖：需要 survival 包（位于 Suggests）；未安装时给出明确提示。
#' 方法：在公式内用 \code{strata(配对变量)} 指定配对层，\code{survival::clogit()}
#' 在每层内做条件 logistic 似然估计。
#' 参数：\code{data} 为数据框；\code{formula} 为模型公式，响应变量为 0/1（1=病例），
#' 必须含 \code{strata(配对变量)} 项。
#' 结果：返回系数表，含估计值、标准误、z 值、p 值、OR（exp(估计)）与 OR 的
#' 95% 置信区间。OR>1 表示该暴露为危险因素。
#' 注意：每个配对层内应同时包含病例与对照；若某层只有一种状态，该层对条件
#' 似然无贡献。模型比较可配合 mc_model_compare() 使用。
#' @param data 数据框（a `data.frame`）。
#' @param formula 含 `strata(配对变量)` 的模型公式（a formula with a
#'   `strata()` term）。
#' @return 系数表（估计 / 标准误 / z / p / OR / OR 置信区间）
#'   （a coefficient table with OR and its 95% CI）。
#' @examples
#' # 第 5 章 5.10 节：以内置配对病例对照数据集 infert 演示条件 Logistic
#' # case：是否不孕(1=病例)；stratum：配对组；spontaneous/induced：流产史
#' infert <- datasets::infert
#' mc_clogistic(infert, case ~ spontaneous + induced + strata(stratum))
#' @export
mc_clogistic <- function(data, formula) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("`mc_clogistic` requires the survival package; please install.packages('survival').")
  }
  fit <- survival::clogit(formula, data = data)
  beta <- as.numeric(fit$coefficients)
  nam  <- names(fit$coefficients)
  se   <- sqrt(diag(as.matrix(fit$var)))
  z    <- beta / se
  p    <- 2 * stats::pnorm(-abs(z))
  or   <- exp(beta)
  ci   <- exp(beta + c(-1, 1) * stats::qnorm(0.975) * se)
  data.frame(variable = nam, estimate = beta, std_error = se, z = z,
             p.value = p, OR = or, OR_lower = ci[1], OR_upper = ci[2],
             stringsAsFactors = FALSE)
}

#' 回归模型比较（自研 AIC/BIC/AICc/权重）/ Model comparison
#'
#' 对应第 5 章「多个不同逻辑回归模型的比较和找出优化模型」（5.8）。对一组
#' 已拟合的模型，计算 AIC、BIC、**自研**的 AICc（小样本校正）以及基于 AICc 的
#' Akaike 权重，并按权重排序，帮助选出相对最优模型。不依赖任何额外包，仅用
#' 基础 R 的 \code{AIC}/\code{BIC}/\code{nobs}。
#'
#' Corresponds to Chapter 5 (5.8, comparing multiple models). Computes AIC, BIC,
#' a *self-implemented* small-sample-corrected AICc and Akaike weights from AICc,
#' then ranks models by weight.
#'
#' @details
#' 用途：在多个候选回归模型间做信息准则比较，选出相对最优者。
#' 方法（自研）：AIC、BIC 取自基础 R；AICc = AIC + 2k(k+1)/(n-k-1)（k 为参数数、
#' n 为样本量，对样本量小于参数数+1 的模型记为 NA）；再计算 delta = AICc - min，
#' Akaike 权重 w = exp(-0.5*delta) / sum(exp(-0.5*delta))，权重越大模型越优。
#' 参数：\code{...} 为若干已拟合模型对象（如 lm/glm/coxph 等，需支持
#' \code{AIC()}/\code{BIC()}/\code{nobs()}）；\code{n} 为样本量（缺省时尝试用
#' \code{nobs()} 自动获取）。
#' 结果：返回按权重降序排列的模型比较表，含 k（参数数）、n、AIC、BIC、AICc、
#' delta_AICc、weight。
#' 注意：比较的模型必须基于同一数据集与同一响应变量，否则准则不可比。权重反映
#' 各模型相对最优的概率，不表示绝对拟合优度。
#' @param ... 已拟合的模型对象（fitted model objects）。
#' @param n 样本量；缺省时自动从模型获取（sample size; auto-detected if NULL）。
#' @return 按 Akaike 权重降序的模型比较表（a ranked comparison table）。
#' @examples
#' # 第 5 章 5.8 节：比较预测出生体重的多个线性模型（AIC/BIC/AICc/权重）
#' birthwt <- MASS::birthwt
#' m1 <- lm(bwt ~ lwt, data = birthwt)
#' m2 <- lm(bwt ~ lwt + age, data = birthwt)
#' m3 <- lm(bwt ~ lwt + age + smoke + factor(race), data = birthwt)
#' mc_model_compare(m1, m2, m3)
#' @export
mc_model_compare <- function(..., n = NULL) {
  models <- list(...)
  if (length(models) < 2) stop("`mc_model_compare` needs at least two models.")
  nm <- if (is.null(names(models))) {
    paste0("model", seq_along(models))
  } else {
    ifelse(nzchar(names(models)), names(models), paste0("model", seq_along(models)))
  }
  aic <- sapply(models, function(m) stats::AIC(m))
  bic <- sapply(models, function(m) stats::BIC(m))
  k   <- sapply(models, function(m) length(stats::coef(m)))
  if (is.null(n)) {
    nob <- tryCatch(sapply(models, function(m) stats::nobs(m)), error = function(e) NULL)
    n <- if (is.null(nob)) stop("无法自动获取样本量，请显式传入 n。") else nob[1]
  }
  aicc <- ifelse(n - k - 1 > 0, aic + 2 * k * (k + 1) / (n - k - 1), NA_real_)
  delta <- aicc - min(aicc, na.rm = TRUE)
  w <- if (any(is.na(aicc))) {
    rep(NA_real_, length(aicc))
  } else {
    e <- exp(-0.5 * delta)
    e / sum(e, na.rm = TRUE)
  }
  tab <- data.frame(model = nm, k = k, n = n, AIC = aic, BIC = bic,
                    AICc = aicc, delta_AICc = delta, weight = w,
                    stringsAsFactors = FALSE)
  tab[order(-tab$weight), ]
}
