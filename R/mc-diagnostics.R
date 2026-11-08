# =============================================================================
# macodgent - shared methodological check layer
# One definition of PASS / REVIEW used by mc_linear(), mc_logistic(),
# mc_auto_model() and the automated workflow. Internal only.
# =============================================================================

.mc_check_row <- function(check, value, threshold, ok, note = "") {
  status <- if (length(ok) != 1L || is.na(ok)) "REVIEW" else
    if (isTRUE(ok)) "PASS" else "REVIEW"
  data.frame(check = check, value = if (length(value)) as.numeric(value) else NA_real_,
             threshold = threshold, status = status, note = note,
             stringsAsFactors = FALSE)
}

.mc_max_vif <- function(fit) {
  X <- tryCatch(stats::model.matrix(fit), error = function(e) NULL)
  if (is.null(X)) return(NA_real_)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (ncol(X) < 2L) return(1)
  keep <- apply(X, 2L, function(v) length(unique(v)) > 1L)
  X <- X[, keep, drop = FALSE]
  if (ncol(X) < 2L) return(1)
  cc <- tryCatch(stats::cor(X), error = function(e) NULL)
  if (is.null(cc)) return(NA_real_)
  d <- suppressWarnings(tryCatch(diag(solve(cc)),
                                 error = function(e) rep(NA_real_, ncol(X))))
  d <- d[is.finite(d) & d > 0]
  if (!length(d)) NA_real_ else max(d)
}

.mc_model_checks <- function(fit, family = c("linear", "logistic", "cox"),
                             alpha = 0.05) {
  family <- match.arg(family)
  zc <- stats::qnorm(1 - alpha / 2)
  if (family == "linear") {
    s <- summary(fit)
    res <- as.numeric(stats::residuals(fit))
    fv <- as.numeric(stats::fitted(fit))
    sh_p <- if (length(res) >= 3L && length(res) <= 5000L)
      tryCatch(as.numeric(stats::shapiro.test(res)$p.value),
               error = function(e) NA_real_) else NA_real_
    bp_p <- tryCatch({
      aux <- stats::lm(e2 ~ fv, data = data.frame(e2 = res^2, fv = fv))
      stats::pchisq(length(res) * summary(aux)$r.squared, df = 1,
                    lower.tail = FALSE)
    }, error = function(e) NA_real_)
    vif <- .mc_max_vif(fit)
    cd <- tryCatch(as.numeric(stats::cooks.distance(fit)),
                   error = function(e) rep(NA_real_, length(res)))
    nflag <- sum(cd > 4 / length(cd), na.rm = TRUE)
    f_p <- tryCatch(as.numeric(stats::pf(s$fstatistic[1], s$fstatistic[2],
                                         s$fstatistic[3], lower.tail = FALSE)),
                    error = function(e) NA_real_)
    out <- rbind(
      .mc_check_row("Residual normality", sh_p, "Shapiro p > alpha",
                    if (is.na(sh_p)) NA else sh_p > alpha,
                    "Shapiro-Wilk; over-sensitive at large n"),
      .mc_check_row("Homoscedasticity", bp_p, "Breusch-Pagan p > alpha",
                    if (is.na(bp_p)) NA else bp_p > alpha),
      .mc_check_row("Multicollinearity", vif, "max VIF < 5",
                    if (is.na(vif)) NA else vif < 5),
      .mc_check_row("Influential observations", nflag, "flagged = 0",
                    nflag == 0L, "Cook's D > 4/n"),
      .mc_check_row("Global model fit", f_p, "F test p < alpha",
                    if (is.na(f_p)) NA else f_p < alpha))
    attr(out, "max_vif") <- vif
    return(out)
  }
  if (family == "logistic") {
    mf <- stats::model.frame(fit)
    yraw <- stats::model.response(mf)
    y <- if (is.factor(yraw) || is.character(yraw))
      as.integer(factor(yraw)) - 1L else as.numeric(yraw)
    p <- as.numeric(stats::fitted(fit))
    ok <- is.finite(p) & is.finite(y)
    p <- p[ok]; y <- y[ok]
    auc <- tryCatch(.mc_auc(p, y), error = function(e) NA_real_)
    n_coef <- length(stats::coef(fit))
    epv <- min(sum(y == 1), sum(y == 0)) / max(1L, n_coef - 1L)
    vif <- .mc_max_vif(fit)
    lrt_p <- tryCatch({
      f0_form <- stats::as.formula(
        paste(deparse(stats::formula(fit)[[2]]), "~ 1"))
      f0 <- stats::glm(f0_form, data = mf, family = stats::binomial)
      a <- stats::anova(f0, fit, test = "Chisq")
      as.numeric(a[2, "Pr(>Chi)"])
    }, error = function(e) NA_real_)
    return(rbind(
      .mc_check_row("Discrimination", auc, "AUC >= 0.70",
                    if (is.na(auc)) NA else auc >= 0.70),
      .mc_check_row("Events per variable", epv, "EPV >= 10",
                    epv >= 10),
      .mc_check_row("Multicollinearity", vif, "max VIF < 5",
                    if (is.na(vif)) NA else vif < 5),
      .mc_check_row("Global model fit", lrt_p, "LRT p < alpha",
                    if (is.na(lrt_p)) NA else lrt_p < alpha)))
  }
  ph_p <- tryCatch({
    z <- survival::cox.zph(fit)
    as.numeric(z$table["GLOBAL", "p"])
  }, error = function(e) NA_real_)
  conc <- tryCatch(as.numeric(summary(fit)$concordance[1]),
                   error = function(e) NA_real_)
  rbind(
    .mc_check_row("Proportional hazards", ph_p, "global PH p > alpha",
                  if (is.na(ph_p)) NA else ph_p > alpha),
    .mc_check_row("Discrimination", conc, "concordance >= 0.60",
                  if (is.na(conc)) NA else conc >= 0.60))
}
