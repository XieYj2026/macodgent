# macodgent

**MA Codgent Toolkit for Routine Medical Data Science Analysis**  
**医学数据科学例行性分析工具箱**

A curated collection of reusable analysis functions distilled from the second  
part ("Routine Analysis Techniques") of the routine medical-data-science tutorial. Every exported function uses the `mc_` prefix and returns a tidy  
`data.frame`, so results drop straight into reports and downstream pipelines.

一套从医学数据科学例行性分析技术提炼出的可复用  
分析函数集合。所有导出函数均以 `mc_` 为前缀，并返回整洁的 `data.frame`，  
可直接落入报告与下游流程。

---

## 设计原则 / Design principles

1. 所有导出函数以 `mc_` 前缀命名，避免与基础 R 及其它包冲突。  
   Every exported function is prefixed with `mc_` to avoid name clashes.
2. 分析函数一律返回整洁 `data.frame`（而非仅打印输出），便于报告与管道。  
   Analysis functions return clean `data.frame` objects, ready for reporting.
3. 重方法（GAM、Lasso、竞争风险、多状态等）放在 `Suggests`，通过  
   `requireNamespace()` 优雅降级；核心仅依赖 `tidyr`、`survival`、`ggplot2`、  
   `broom`、`stats`、`utils`。  
   Heavy methods live in `Suggests` with graceful namespace guards; the core  
   depends only on `tidyr`, `survival`, `ggplot2`, `broom`, `stats`, `utils`.

---

## 安装 / Installation

```r
# 从源码安装 / Install from source
# R CMD INSTALL macodgent
# 或使用 devtools / or with devtools:
# devtools::install("macodgent")

# 载入 / Load
library(macodgent)
```

---

## 模块与函数速查表 / Modules & function index

共 33 个导出函数，按章节分为 8 个模块。/ 33 exported functions in 8 modules.

### 第 1 章 认识数据 / Meet the data

| 函数 Function    | 用途 Purpose                                     |
| -------------- | ---------------------------------------------- |
| `mc_mode()`    | 众数（支持多众数）/ Mode, supports multiple modes       |
| `mc_inspect()` | 数据框逐列概览（类别/缺失/唯一值）/ Per-column data overview   |
| `mc_theme()`   | MA Codgent 浅色 ggplot2 主题 / Light ggplot2 theme |
| `mc_palette()` | MA Codgent 品牌配色 / Brand colour palette         |

### 第 2 章 预处理 / Preprocessing

| 函数 Function        | 用途 Purpose                                |
| ------------------ | ----------------------------------------- |
| `mc_impute_knn()`  | K 近邻填补（自实现，无额外依赖）/ KNN imputation         |
| `mc_impute_mice()` | 多重填补（MICE 封装）/ Multiple imputation (MICE) |
| `mc_dedup()`       | 去重（按键列或全列）/ Remove duplicates             |
| `mc_scale()`       | 特征缩放（z / minmax / 稳健）/ Feature scaling    |
| `mc_wide_long()`   | 宽表转长表（tidyr 封装）/ Wide to long             |
| `mc_long_wide()`   | 长表转宽表（tidyr 封装）/ Long to wide             |

### 第 3 章 定量描述 / Quantitative description

| 函数 Function       | 用途 Purpose                                             |
| ----------------- | ------------------------------------------------------ |
| `mc_desc_quant()` | 一站式定量描述（约 20 个统计量，支持分组）/ One-shot ~20 metrics          |
| `mc_norm_test()`  | 正态性检验（Shapiro-Wilk / KS）/ Normality test               |
| `mc_var_test()`   | 方差齐性检验（Bartlett/Fligner/Levene）/ Variance homogeneity  |
| `mc_ttest()`      | t 检验统一入口（单样本/独立/配对）/ Unified t-test                    |
| `mc_anova()`      | 单因素方差分析 / One-way ANOVA                                |
| `mc_nonparam()`   | 非参数组间检验（Wilcoxon/Kruskal-Wallis）/ Non-parametric tests |

### 第 4 章 定性描述 / Qualitative description

| 函数 Function        | 用途 Purpose                                 |
| ------------------ | ------------------------------------------ |
| `mc_contingency()` | 列联表与卡方检验（2×2 含 OR）/ Contingency table + OR |
| `mc_corr()`        | 相关分析（含偏相关）/ Correlation (with partial)     |
| `mc_tableone()`    | 基线特征表 Table One（支持分层）/ Baseline table      |

### 第 5 章 常规回归 / Conventional regression

| 函数 Function             | 用途 Purpose                                   |
| ----------------------- | -------------------------------------------- |
| `mc_logistic()`         | 二分类 logistic 回归（含 OR 与 CI）/ Binary logistic  |
| `mc_logistic_screen() ` | 单因素筛选 → 多因素（先单后多）/ Univariate→multivariable  |
| `mc_stepwise()`         | 逐步变量选择（基于 AIC）/ Stepwise selection           |
| `mc_ordinal_logistic()` | 有序 logistic 回归（MASS::polr）/ Ordinal logistic |

### 第 6 章 生存分析 / Survival analysis

| 函数 Function    | 用途 Purpose                                            |
| -------------- | ----------------------------------------------------- |
| `mc_survfit()` | Kaplan-Meier 曲线 + log-rank 检验与作图/ KM curve + log-rank |
| `mc_cox()`     | Cox 模型（PH 检验 + 森林图）/ Cox PH + forest plot             |

### 第 7 章 进阶回归 / Advanced regression

| 函数 Function         | 用途 Purpose                                   |
| ------------------- | -------------------------------------------- |
| `mc_gam()`          | 广义可加模型 GAM（mgcv）/ Generalized additive model |
| `mc_lasso()`        | Lasso 回归（glmnet, alpha=1）/ Lasso             |
| `mc_ridge()`        | 岭回归（glmnet, alpha=0）/ Ridge                  |
| `mc_elastic_net()`  | 弹性网络回归（glmnet）/ Elastic net                  |
| `mc_rcs()`          | 受限立方样条回归（自实现基函数）/ Restricted cubic spline    |
| `mc_quantile_reg()` | 分位数回归（quantreg）/ Quantile regression         |

### 第 8 章 进阶生存 / Advanced survival

| 函数 Function           | 用途 Purpose                                            |
| --------------------- | ----------------------------------------------------- |
| `mc_competing_risk()` | 竞争风险回归（Fine-Gray, cmprsk）/ Competing-risk (Fine-Gray) |
| `mc_multistate()`     | 多状态生存模型（mstate）/ Multi-state model                    |

---

## 快速上手 / Quick start

```r
library(macodgent)
library(MASS)   # 内含示例数据集 birthwt

# 1) 数据概览 / Data overview
mc_inspect(birthwt)

# 2) 一站式定量描述（约 20 个统计量）/ One-shot quantitative description
mc_desc_quant(birthwt, vars = c("age", "lwt", "bwt"))
mc_desc_quant(birthwt, vars = c("age", "lwt"), group = "smoke")

# 3) 基线特征表 Table One（按 smoke 分层）/ Baseline table, stratified
mc_tableone(birthwt, vars = c("age", "lwt", "smoke", "low"), group = "smoke")

# 4) 二分类 logistic 回归（OR 与 95% CI）/ Binary logistic regression
bw <- birthwt
bw$low <- factor(bw$low)
mc_logistic(bw, low ~ age + smoke)

# 5) 单因素筛选 → 多因素（先单后多）/ Univariate screen -> multivariable
mc_logistic_screen(bw, y = "low",
  x_vars = c("age", "lwt", "race", "smoke", "ht", "ui"))

# 6) 受限立方样条（自实现基函数，无需 rms）/ Restricted cubic spline
rcs <- mc_rcs(birthwt, x = "age", y = "bwt", nk = 3)
rcs$plot

# 7) Cox 模型 + 森林图 / Cox model + forest plot
fit <- mc_cox(survival::lung, survival::Surv(time, status) ~ age + sex + ph.ecog)
fit$forest
```

---

## 查看帮助 / Reading the help

每个函数都带有详细的**中英双语**帮助页（描述、参数、返回值、可运行示例、  
参见、参考文献）。/ Every function ships with a detailed bilingual help page.

```r
?mc_desc_quant          # 单个函数帮助 / help for one function
help(package = "macodgent")   # 包级索引 / package index
```

在 RStudio 中直接输入函数名按 F1 亦可。/ Press F1 on a function name in RStudio.

---

## 许可 / License

MIT © MA Codgent OPC
