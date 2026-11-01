# =============================================================================
# 软件名称    : macodgent 医学数据科学研究工具包
# 著作权人    : Dr.Xie
# 版本号      : 1.0.0
# 开发语言    : R
# 开发完成日期: 2026-10-30
# 权利声明    : 本程序为著作权人独立开发，受《计算机软件保护条例》保护；
#              未经书面许可，不得复制、发行、出租或将其用于商业目的。
# -----------------------------------------------------------------------------
# 模块        : 生存分析模块（survival）
# 本文件包含  : mc_survfit / mc_cox / mc_competing_risk / mc_multistate /
#              mc_ph_test（Cox 比例风险假定检验）
#              mc_competing_risk 自行由 crr 系数重建子分布风险比表与 p 值；
#              mc_multistate 自行从公式解析 strata 并构造多状态预测。
# =============================================================================

#' Kaplan-Meier 生存曲线与 log-rank 检验 / KM survival curve & log-rank test
#'
#' 对应第 6 章：生存估计、组间 log-rank 检验与可视化。无分组时估计整体
#' 生存曲线，有分组时同时给出组间检验 p 值与按 \code{mc_theme()} 主题的 KM 图。
#'
#' Corresponds to Chapter 6: survival estimation, between-group log-rank test,
#' and visualisation.
#'
#' @details
#' 用途：估计 Kaplan-Meier 生存曲线，并用 log-rank 检验比较组间生存差异，是生存分析最基础的探索工具。
#' 方法：KM 法基于条件生存概率乘积估计生存函数；有分组时自动做 log-rank 检验（检验各组生存曲线是否总体相同）。
#' 参数：data 为数据框；time 为随访时间列名；status 为事件指示列名（通常 1=发生事件，0=删失）；group 为可选的分组列。
#' 结果：返回列表，含 survfit（生存对象）、logrank（log-rank 检验的卡方、自由度、p 值，有分组时）、以及 plot（ggplot 绘制的 KM 曲线，含置信带与分组图例）。
#' 注意：删失数据（status=0）不能当作「事件」计入；log-rank 检验对晚期差异不敏感（若想侧重远期差异可用加权检验）。绘制结果可直接 print(res$plot) 查看。
#' @param data 数据框（a `data.frame`）。
#' @param time 时间变量名（time variable name）。
#' @param status 结局变量名（0 = 删失，1 = 事件）（outcome variable name;
#'   0 = censored, 1 = event）。
#' @param group 可选的分组变量名（optional grouping variable name）。
#' @return 列表，含 `fit`（survfit 对象）/ `logrank`（检验数据框：test/chisq/
#'   p.value；无分组时为 `NULL`）/ `plot`（ggplot KM 图）（a list with
#'   \code{fit}, \code{logrank}, \code{plot}）。
#' @examples
#' # 第 6 章：以卵巢癌预后数据集 ovarian 演示 KM 曲线与 log-rank
#' ovarian <- survival::ovarian
#' ov <- ovarian
#' ov$status <- ifelse(ov$fustat == 2, 0, 1)   # 2=删失->0, 1=死亡->1
#' res <- mc_survfit(ov, time = "futime", status = "status", group = "rx")
#' cat("---- log-rank 检验 / test ----\n"); print(res$logrank, row.names = FALSE)
#' print(res$plot)
#' @export
mc_survfit <- function(data, time, status, group = NULL) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("The survival package is required.")
  sv <- survival::Surv(data[[time]], data[[status]])
  if (is.null(group)) {
    fit <- survival::survfit(sv ~ 1)
    lr_df <- NULL
    td <- broom::tidy(fit)
    td$strata <- "All"
  } else {
    g <- data[[group]]
    tmp <- data.frame(sv = sv, g = g)
    fit <- survival::survfit(sv ~ g, data = tmp)
    lr <- survival::survdiff(sv ~ g, data = tmp)
    lr_df <- data.frame(test = "log-rank",
                        chisq = as.numeric(lr$chisq),
                        p.value = as.numeric(lr$pvalue),
                        stringsAsFactors = FALSE)
    td <- broom::tidy(fit)
  }
  plt <- ggplot2::ggplot(td, ggplot2::aes(x = time, y = estimate, color = strata)) +
    ggplot2::geom_step() +
    ggplot2::scale_x_continuous(breaks = .nice_breaks(td$time, 6)) +
    ggplot2::labs(x = "Time", y = "Survival probability", color = "Group") +
    mc_theme()
  if (!is.null(lr_df)) {
    plt <- plt + ggplot2::ggtitle(paste0("log-rank p = ",
                                         format.pval(lr_df$p.value, digits = 3)))
  }
  list(fit = fit, logrank = lr_df, plot = plt)
}

#' Cox 比例风险模型（含 PH 检验与森林图）/ Cox PH model (PH test & forest plot)
#'
#' 对应第 6 章：Cox 模型、比例风险假设检验（Schoenfeld 残差）与森林图。
#' 返回 HR 表、PH 检验结果，以及按 \code{mc_theme()} 主题的森林图。
#'
#' Corresponds to Chapter 6: Cox model, proportional-hazards assumption test
#' (Schoenfeld residuals), and a forest plot.
#'
#' @details
#' 用途：拟合 Cox 比例风险模型，评估多个协变量对发生事件风险的同期影响，得到风险比（HR）。
#' 方法：半参数模型，不假设基线风险形态，只估计协变量对应的风险比；默认对比例风险（PH）假设做检验。
#' 参数：data 为数据框；formula 形如 Surv(time, status) ~ 协变量，需用 survival::Surv() 构造生存对象；conf 为 HR 置信水平。
#' 结果：返回列表，含 coefficients（各变量 HR、标准误、z、p 值与置信区间）、ph_test（比例风险假设检验结果）、以及 forest（HR 森林图，ggplot 对象）。
#' 注意：HR>1 表示该协变量升高时事件风险升高。PH 假设是 Cox 的核心前提——若 ph_test 提示显著违背（如某变量 p 很小），需考虑分层 Cox、时变效应或改参数模型。森林图用 print(res$forest) 查看。
#' @param data 数据框（a `data.frame`）。
#' @param formula Cox 模型公式，如 \code{Surv(time, status) ~ x1 + x2}（a Cox
#'   model formula）。
#' @param conf 置信水平，默认 `0.95`（confidence level; default `0.95`）。
#' @return 列表，含 `model`（coxph 对象）/ `coefficients`（HR 表：term/estimate/
#'   conf.low/conf.high/p.value）/ `ph_test`（PH 检验表）/ `forest`（ggplot
#'   森林图）（a list with \code{model}, \code{coefficients}, \code{ph_test},
#'   \code{forest}）。
#' @examples
#' # 第 6 章：Cox 比例风险模型——卵巢癌生存预后
#' ovarian <- survival::ovarian
#' ov <- ovarian
#' ov$status <- ifelse(ov$fustat == 2, 0, 1)
#' res <- mc_cox(ov, survival::Surv(futime, status) ~ age + rx + resid.ds + ecog.ps)
#' cat("---- HR 表 / coefficients ----\n"); print(as.data.frame(res$coefficients), row.names = FALSE)
#' cat("\n---- PH 假设检验 / PH test ----\n"); print(res$ph_test)
#' print(res$forest)
#' @export
mc_cox <- function(data, formula, conf = 0.95) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("The survival package is required.")
  fit <- eval(substitute(survival::coxph(formula, data = DATA),
                        list(DATA = data)))
  td <- broom::tidy(fit, conf.int = TRUE, conf.level = conf, exponentiate = TRUE)
  td$variable <- td$term
  zph <- survival::cox.zph(fit)
  zph_df <- as.data.frame(zph$table)
  fp <- ggplot2::ggplot(td, ggplot2::aes(x = variable, y = estimate)) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = conf.low, ymax = conf.high)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#C0392B") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "", y = "Hazard Ratio (95% CI)") +
    mc_theme()
  list(model = fit, coefficients = td, ph_test = zph_df, forest = fp)
}

#' 竞争风险回归（Fine-Gray 模型）/ Competing-risk regression (Fine-Gray)
#'
#' 对应第 8 章「竞争风险生存分析」。基于 \pkg{cmprsk} 的 \code{crr}。
#' `status`：0 = 删失，>0 = 不同事件类型；`cause` 选择关注事件。
#'
#' Corresponds to Chapter 8 ("competing-risk survival analysis"). Based on
#' \pkg{cmprsk}'s \code{crr}. \code{status}: 0 = censored, >0 = different event
#' types; \code{cause} selects the event of interest.
#'
#' @details
#' 用途：在存在「竞争事件」时做回归，避免把竞争事件误当作删失而高估目标事件风险。典型场景：研究死亡为竞争事件的复发分析。
#' 方法：基于 cmprsk 包的 Fine-Gray 子分布风险模型（subdistribution hazard），将竞争事件作为特殊的删失处理，直接估计目标原因累积发生率。
#' 参数：data 为数据框；time 为时间列；status 为状态列（0=删失，>0 代表不同事件类型）；group 为分组/暴露变量；cause 指定关注的目标事件编号（默认 1）。
#' 结果：返回列表，含 coefficients（各变量对应的子分布风险比 SHR、标准误、z、p 值）。
#' 注意：status 的编码约定与 mc_survfit()/mc_cox() 不同：这里 0 才是删失，其它正整数各代表一类事件，务必按此约定整理数据。若仅有一个事件类型且无竞争，普通 Cox 即可。依赖 cmprsk（Suggests）。
#' @param data 数据框（a `data.frame`）。
#' @param time 时间变量名（time variable name）。
#' @param status 竞争风险状态变量名（competing-risk status variable name）。
#' @param group 协变量（分组/暴露）变量名（covariate/grouping variable name）。
#' @param cause 关注事件类型，默认 `1`（event type of interest; default `1`）。
#' @return 列表，含 `model`（crr 对象）/ `coefficients`（亚分布 HR 表：variable/
#'   HR/coef/p.value）（a list with \code{model}, \code{coefficients}）。
#' @examples
#' \donttest{
#' # 第 8 章 8.2 节：白血病复发(1)与死亡(2)的竞争风险
#' # 关注事件 cause = 1（复发）；死亡为竞争事件，不能当作删失
#' set.seed(126)
#' n <- 300
#' trt     <- sample(c(0,1), n, TRUE)
#' ftime   <- rexp(n, .05)
#' fstatus <- sample(0:2, n, TRUE, c(.4,.35,.25))   # 0删失 1复发 2死亡(竞争)
#' df <- data.frame(ftime = ftime, fstatus = fstatus, trt = trt)
#' res <- mc_competing_risk(df, time = "ftime", status = "fstatus", group = "trt", cause = 1)
#' print(res$coefficients)
#' }
#' @export
mc_competing_risk <- function(data, time, status, group, cause = 1) {
  if (!requireNamespace("cmprsk", quietly = TRUE)) {
    stop("`mc_competing_risk` requires the cmprsk package; please install.packages('cmprsk').")
  }
  ftime <- data[[time]]
  fstatus <- data[[status]]
  cov <- as.matrix(data[[group]])
  colnames(cov) <- group
  fit <- cmprsk::crr(ftime, fstatus, cov1 = cov, failcode = cause)
  cf <- as.numeric(fit$coef)
  names(cf) <- names(fit$coef)
  se <- sqrt(diag(as.matrix(fit$var)))
  z <- cf / se
  pval <- 2 * stats::pnorm(-abs(z))
  coef_df <- data.frame(
    variable = names(cf),
    coef = cf,
    HR = exp(cf),
    se = se,
    p.value = pval, stringsAsFactors = FALSE)
  list(model = fit, coefficients = coef_df)
}

#' 多状态生存模型 / Multi-state survival model
#'
#' 对应第 8 章「多状态生存模型」。要求 `data` 为 \pkg{mstate} 展开的长格式
#' （含 Tstart、Tstop、status、from、to 与协变量）；`trans` 为转移矩阵。
#'
#' Corresponds to Chapter 8 ("multi-state survival model"). Requires \code{data}
#' to be in \pkg{mstate}'s expanded long format (with Tstart, Tstop, status,
#' from, to and covariates); \code{trans} is the transition matrix.
#'
#' @details
#' 用途：拟合多状态生存模型，研究个体在多个健康状态间（如健康→并发症→死亡）的转移风险，比单一终点更贴近临床过程。
#' 方法：基于 mstate 包，先把随访数据展开为「起点-终点」的转移长表，再用 coxph() 对每条转移分别建模。
#' 参数：data 为已展开的多状态长表（含 id、Tstart、Tstop、status、from、to 等列）；trans 为转移矩阵（由 mstate::transMat() 构造，定义状态间允许哪些转移）；formula 含生存对象与协变量，可用 strata(trans) 对转移类型分层。
#' 结果：返回列表，含 model（各转移的 coxph 拟合）、msfit（用于预测累积转移概率的 msfit 对象）。
#' 注意：输入数据必须已是 mstate 要求的长表格式（非原始宽表），trans 的编号必须与 from/to 列及 strata(trans) 一致；构造错误会导致「所有观测同组」或无转移可估。依赖 mstate（Suggests）。
#' @param data mstate 展开的长格式数据框（mstate expanded long-format
#'   `data.frame`）。
#' @param trans 转移矩阵（transition matrix）。
#' @param formula 公式，如 \code{Surv(time, status) ~ cov}（time/status 为 data
#'   中的列）（a formula; time/status are columns in data）。
#' @param ... 传给 \code{survival::coxph} 的额外参数（extra arguments passed to
#'   \code{survival::coxph}）。
#' @return 列表，含 `model`（coxph 对象）/ `msfit`（mstate::msfit 对象）（a list
#'   with \code{model}, \code{msfit}）。
#' @examples
#' \donttest{
#' # 第 8 章：三状态模型（健康 H -> 并发症 I / 直接死亡 D）
#' set.seed(127)
#' library(mstate)
#' n <- 150
#' age <- rnorm(n, 55, 8)
#' t12 <- rexp(n, .15); t23 <- rexp(n, .20); t13 <- rexp(n, .10)
#' rows <- list()
#' for (i in 1:n) {
#'   if (runif(1) < .6) {
#'     rows[[length(rows)+1]] <- data.frame(id=i, Tstart=0,      Tstop=t12[i],        status=1, from=1, to=2, age=age[i])
#'     rows[[length(rows)+1]] <- data.frame(id=i, Tstart=t12[i], Tstop=t12[i]+t23[i], status=1, from=2, to=3, age=age[i])
#'   } else {
#'     rows[[length(rows)+1]] <- data.frame(id=i, Tstart=0, Tstop=t13[i], status=1, from=1, to=3, age=age[i])
#'   }
#' }
#' ms <- do.call(rbind, rows)
#' trans <- transMat(list(c(2,3), c(3), c()), names = c("H","I","D"))
#' ms$trans <- with(ms, ifelse(from == 1 & to == 2, 1, ifelse(from == 1 & to == 3, 2, 3)))
#' fit <- mc_multistate(ms, trans = trans, formula = survival::Surv(Tstart, Tstop, status) ~ age + strata(trans))
#' cat("coxph 系数 / coxph coefficients:\n")
#' print(summary(fit$model)$coefficients)
#' cat("\nmsfit 类别 / class:", class(fit$msfit)[1], "\n")
#' }
#' @export
mc_multistate <- function(data, trans, formula, ...) {
  if (!requireNamespace("mstate", quietly = TRUE)) {
    stop("`mc_multistate` requires the mstate package; please install.packages('mstate').")
  }
  fit <- eval(substitute(survival::coxph(formula, data = DATA, ...),
                        list(DATA = data)))
  nd <- data
  fstr <- paste(deparse(formula), collapse = "")
  m <- regmatches(fstr, gregexpr("strata\\(([^)]+)\\)", fstr))[[1]]
  if (length(m)) {
    svar <- sub("strata\\(([^)]+)\\)", "\\1", m[1])
    nd$strata <- data[[svar]]
  }
  ms <- mstate::msfit(fit, newdata = nd, trans = trans)
  list(model = fit, msfit = ms)
}

#' Cox 比例风险假定检验（PH 检验）/ Cox PH-assumption test
#'
#' 对应第 6 章「比例风险假定的检验和森林图的绘制」（6.5）。对 Cox 模型做
#' 比例风险假定检验：基于 Schoenfeld 残差，用 \code{survival::cox.zph()} 检验
#' 各协变量系数是否随时间变化。同时给出**自研 ggplot 残差图**（标准化 Schoenfeld
#' 残差对变换后时间），直观判断比例风险是否成立。
#'
#' Corresponds to Chapter 6 (6.5, testing the PH assumption). Runs
#' \code{survival::cox.zph()} on a Cox model and adds a *self-built* ggplot of
#' scaled Schoenfeld residuals versus transformed time.
#'
#' @details
#' 用途：检验 Cox 比例风险模型的核心假设——各协变量的风险比是否随时间恒定。
#' 依赖：需要 survival 包（已在 Imports 中）。
#' 方法：Schoenfeld 残差应随时间无趋势；\code{survival::cox.zph()} 对每条曲线做
#' 平滑后做相关检验，p>0.05 通常认为满足 PH 假定。
#' 参数：\code{fit} 为已拟合的 \code{coxph} 对象（可由 mc_cox() 返回）；或传入
#' \code{data} 与 \code{formula} 由函数先拟合。
#' 结果：返回列表——\code{test}（检验表：变量/卡方/df/p 值，含 Global 行）、
#' \code{plot}（标准化 Schoenfeld 残差 vs 时间的 ggplot，可指定 \code{which}
#' 变量，默认第一个协变量）。
#' 注意：样本量较小或事件数较少时 PH 检验功效低；若违反 PH 假定，可考虑分层
#' Cox、加入时间依存项或改用 mc_competing_risk() 等替代模型。
#' @param fit 已拟合的 coxph 对象（a fitted `coxph` object）。
#' @param data 数据框（仅当未提供 fit 时使用）。
#' @param formula 生存公式（仅当未提供 fit 时使用）。
#' @param which 残差图所画的协变量名（默认第一个）。
#' @param transform 时间变换，默认 "km"（Kaplan-Meier）。
#' @return 列表，含检验表与残差图（a list with the test table and residual plot）。
#' @examples
#' # 第 6 章 6.5 节：检验卵巢癌 Cox 模型的比例风险假定
#' ovarian <- survival::ovarian
#' ov <- ovarian
#' ov$status <- ifelse(ov$fustat == 2, 0, 1)
#' fit <- survival::coxph(survival::Surv(futime, status) ~ age + rx + resid.ds, data = ov)
#' r <- mc_ph_test(fit)
#' print(r$test)        # 各协变量 PH 检验表
#' print(r$plot)        # 标准化 Schoenfeld 残差 vs 时间
#' @export
mc_ph_test <- function(fit = NULL, data = NULL, formula = NULL,
                       which = NULL, transform = "km") {
  if (is.null(fit)) {
    if (is.null(data) || is.null(formula))
      stop("Provide either a fitted `coxph` object or both `data` and `formula`.")
    if (!requireNamespace("survival", quietly = TRUE))
      stop("`mc_ph_test` requires the survival package.")
    fit <- survival::coxph(formula, data = data)
  }
  if (!inherits(fit, "coxph"))
    stop("`fit` must be a coxph object.")
  zp <- survival::cox.zph(fit, transform = transform)
  tb <- as.data.frame(zp$table)
  names(tb) <- c("chisq", "df", "p_value")
  tb$variable <- rownames(tb)
  vnames <- colnames(zp$y)   # 协变量名（cox.zph$y 的列=变量，行=事件时间）
  w <- if (is.null(which)) vnames[1] else which
  if (!w %in% vnames) stop("`which` must be one of: ", paste(vnames, collapse = ", "))
  yy <- as.numeric(zp$y[, w, drop = TRUE])
  tt <- as.numeric(zp$x)
  dfp <- data.frame(time = tt, sresid = yy)
  plot <- ggplot2::ggplot(dfp, ggplot2::aes(time, sresid)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.6) +
    ggplot2::geom_smooth(se = TRUE, color = "steelblue", size = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::labs(x = paste("Transformed time (", transform, ")", sep = ""),
                  y = paste("Scaled Schoenfeld residual (", w, ")", sep = "")) +
    mc_theme()
  list(test = tb[, c("variable", "chisq", "df", "p_value")], plot = plot)
}
