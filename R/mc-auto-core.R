# =============================================================================
# Software name : macodgent medical data-science research toolkit
# Copyright holder : Dr.Xie
# Version : 1.0.2
# Language : R
# -----------------------------------------------------------------------------
# Module : automated one-call analysis
# mc_auto_analyze() turns one data frame into a complete, auditable analysis
# bundle: column audit, descriptives, Table One, group comparison, correlation,
# univariable screening, multivariable modelling, model diagnostics, written
# findings and a decision log. mc_auto_report() renders the bundle as a
# self-contained light-theme HTML report. Only base R plus the package's
# existing Imports are used; no external CSS, JavaScript or CDN is required.
# =============================================================================

# ---- internal helpers -------------------------------------------------------

# Quote a column name so that non-syntactic names work inside formulas.
.mc_bt <- function(x) paste0("`", gsub("`", "", x), "`")

.mc_esc <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

.mc_pfmt <- function(p, digits = 3L) {
  vapply(p, function(v) {
    if (is.na(v) || !is.finite(v)) return("NA")
    format.pval(v, digits = digits)
  }, character(1), USE.NAMES = FALSE)
}

.mc_ci_fmt <- function(est, lo, hi) {
  if (any(!is.finite(c(est, lo, hi)))) return("NA")
  sprintf("%.2f (%.2f-%.2f)", est, lo, hi)
}

# Coarse but clinically useful column typing.
.mc_kind <- function(x) {
  if (inherits(x, c("Date", "POSIXt"))) return("date")
  if (is.factor(x) || is.character(x) || is.logical(x)) return("categorical")
  if (is.numeric(x)) {
    u <- length(unique(x[!is.na(x)]))
    if (u <= 2L) return("binary")
    return("numeric")
  }
  "other"
}

.mc_skew <- function(z) {
  n <- length(z)
  if (n < 3L) return(c(NA_real_, NA_real_))
  m <- mean(z)
  v <- mean((z - m)^2)
  if (!is.finite(v) || v <= 0) return(c(NA_real_, NA_real_))
  c(mean((z - m)^3) / v^1.5, mean((z - m)^4) / v^2 - 3)
}

# One row per column, numeric and categorical metrics in the same tidy table.
.mc_desc_row <- function(nm, x, alpha = 0.05) {
  kind <- .mc_kind(x)
  n_miss <- sum(is.na(x))
  xx <- x[!is.na(x)]
  d <- data.frame(
    variable = nm, type = kind, n = length(xx), n_missing = n_miss,
    pct_missing = round(100 * n_miss / max(1L, length(x)), 1),
    n_unique = length(unique(xx)),
    mean = NA_real_, sd = NA_real_, median = NA_real_, q1 = NA_real_,
    q3 = NA_real_, min = NA_real_, max = NA_real_,
    skewness = NA_real_, kurtosis = NA_real_, normality_p = NA_real_,
    distribution = NA_character_, top_level = NA_character_,
    top_pct = NA_real_, stringsAsFactors = FALSE)
  if (kind %in% c("numeric", "binary") && is.numeric(x)) {
    z <- as.numeric(xx)
    q <- as.numeric(stats::quantile(z, c(0.25, 0.5, 0.75), names = FALSE))
    d$mean <- mean(z); d$sd <- stats::sd(z)
    d$median <- q[2]; d$q1 <- q[1]; d$q3 <- q[3]
    d$min <- min(z); d$max <- max(z)
    sk <- .mc_skew(z); d$skewness <- sk[1]; d$kurtosis <- sk[2]
    if (length(z) >= 3L && length(z) <= 5000L && length(unique(z)) > 2L) {
      pv <- tryCatch(as.numeric(stats::shapiro.test(z)$p.value),
                     error = function(e) NA_real_)
      d$normality_p <- pv
      d$distribution <- if (is.na(pv)) "not tested" else
        if (pv >= alpha) "approximately normal" else "non-normal"
    } else {
      d$distribution <- "not tested"
    }
  } else if (kind == "categorical") {
    tb <- sort(table(as.character(xx)), decreasing = TRUE)
    if (length(tb)) {
      d$top_level <- names(tb)[1]
      d$top_pct <- round(100 * as.numeric(tb[1]) / sum(tb), 1)
    }
    d$distribution <- paste0(length(tb), " levels")
  } else {
    d$distribution <- "unsupported"
  }
  d
}

.mc_check <- function(check, value, threshold, status, note = "") {
  data.frame(check = check, value = value, threshold = threshold,
             status = status, note = note, stringsAsFactors = FALSE)
}

.mc_check_status <- function(ok) {
  if (is.na(ok)) "UNDETERMINED" else if (isTRUE(ok)) "PASS" else "REVIEW"
}

# Unified coefficient table for lm / glm / coxph.
.mc_coef <- function(fit, family, conf, p_adjust, alpha) {
  expo <- family != "linear"
  td <- tryCatch(as.data.frame(broom::tidy(
    fit, conf.int = TRUE, conf.level = conf, exponentiate = expo)),
    error = function(e) NULL)
  if (is.null(td) || !"term" %in% names(td)) return(NULL)
  for (k in c("estimate", "std.error", "statistic", "p.value",
              "conf.low", "conf.high"))
    if (!k %in% names(td)) td[[k]] <- NA_real_
  td <- td[as.character(td$term) != "(Intercept)", , drop = FALSE]
  if (!nrow(td)) return(NULL)
  out <- data.frame(
    term = as.character(td$term),
    effect = if (family == "linear") "Beta" else
      if (family == "logistic") "OR" else "HR",
    estimate = as.numeric(td$estimate),
    std.error = as.numeric(td$std.error),
    conf.low = as.numeric(td$conf.low),
    conf.high = as.numeric(td$conf.high),
    p.value = as.numeric(td$p.value), stringsAsFactors = FALSE)
  out$p.adjusted <- as.numeric(stats::p.adjust(out$p.value, method = p_adjust))
  out$significant <- !is.na(out$p.value) & out$p.value < alpha
  out$n <- tryCatch(as.integer(stats::nobs(fit)),
                    error = function(e) NA_integer_)
  rownames(out) <- NULL
  out
}

.mc_vif <- function(fit) {
  X <- tryCatch(stats::model.matrix(fit), error = function(e) NULL)
  if (is.null(X)) return(NA_real_)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (ncol(X) < 2L) return(if (ncol(X) == 1L) 1 else NA_real_)
  cc <- tryCatch(stats::cor(X), error = function(e) NULL)
  if (is.null(cc)) return(NA_real_)
  v <- suppressWarnings(tryCatch(as.numeric(diag(solve(cc))),
                                 error = function(e) rep(NA_real_, ncol(cc))))
  v[!is.finite(v) | v < 0] <- NA_real_
  if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
}

.mc_fit <- function(f, d, family) {
  if (family == "linear") return(stats::lm(f, data = d))
  if (family == "logistic")
    return(stats::glm(f, data = d, family = stats::binomial))
  # Embed formula and data in the call. coxph stores its call, and cox.zph()
  # re-evaluates it after this helper's local bindings have disappeared; a call
  # made of bare symbols would then fail with "object 'd' not found".
  eval(substitute(survival::coxph(FORMULA, data = DATA, x = TRUE),
                  list(FORMULA = f, DATA = d)))
}

# Univariable p value from a likelihood-ratio / F comparison with the null
# model; falls back to the smallest Wald p value if that comparison fails.
.mc_screen_p <- function(d, family, fit, co) {
  p_all <- tryCatch({
    f0 <- stats::update(stats::formula(fit), . ~ 1)
    null_fit <- .mc_fit(f0, d, family)
    a <- if (family == "logistic") stats::anova(null_fit, fit, test = "Chisq") else
      stats::anova(null_fit, fit)
    as.numeric(a[nrow(a), ncol(a)])
  }, error = function(e) NA_real_)
  if (!is.finite(p_all))
    p_all <- suppressWarnings(min(co$p.value, na.rm = TRUE))
  if (!is.finite(p_all)) NA_real_ else p_all
}

# ---- internal helpers: HTML -------------------------------------------------

.mc_html_table <- function(d, max_rows = 30L, digits = 3L) {
  if (is.null(d) || !is.data.frame(d) || nrow(d) == 0L)
    return('<p class="empty">Not available for this dataset.</p>')
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  show <- utils::head(d, max_rows)
  hdr <- paste0("<thead><tr>",
                paste0("<th>", .mc_esc(names(show)), "</th>", collapse = ""),
                "</tr></thead>")
  body <- paste0("<tbody>", paste(vapply(seq_len(nrow(show)), function(i) {
    paste0("<tr>", paste(vapply(seq_along(show), function(j) {
      v <- show[[j]][i]
      txt <- if (is.numeric(v)) {
        if (is.na(v) || !is.finite(v)) "NA" else sprintf("%.*g", digits, v)
      } else .mc_esc(v)
      paste0("<td>", txt, "</td>")
    }, character(1)), collapse = ""), "</tr>")
  }, character(1)), collapse = ""), "</tbody>")
  more <- if (nrow(d) > max_rows)
    sprintf('<p class="empty">Showing %d of %d rows.</p>', max_rows, nrow(d)) else ""
  paste0('<div class="tw"><table>', hdr, body, '</table></div>', more)
}

.mc_badge <- function(s) {
  cls <- switch(as.character(s), PASS = "ok", REVIEW = "warn", "muted")
  paste0('<span class="badge ', cls, '">', .mc_esc(s), "</span>")
}

.mc_checks_html <- function(checks) {
  if (is.null(checks) || !nrow(checks))
    return('<p class="empty">No model checks available.</p>')
  rows <- vapply(seq_len(nrow(checks)), function(i) {
    v <- checks$value[i]
    vt <- if (is.numeric(v)) {
      if (is.na(v) || !is.finite(v)) "NA" else sprintf("%.3g", v)
    } else .mc_esc(v)
    paste0("<tr><td>", .mc_esc(checks$check[i]), "</td><td>", vt,
           "</td><td>", .mc_esc(checks$threshold[i]), "</td><td>",
           .mc_badge(checks$status[i]), "</td><td class=\"note\">",
           .mc_esc(checks$note[i]), "</td></tr>")
  }, character(1))
  paste0('<div class="tw"><table><thead><tr><th>Check</th><th>Value</th>',
         '<th>Threshold</th><th>Status</th><th>Note</th></tr></thead><tbody>',
         paste(rows, collapse = ""), "</tbody></table></div>")
}

# A dependency-free forest plot for the model coefficient table.
.mc_forest <- function(co, effect) {
  if (is.null(co) || !nrow(co)) return("")
  logscale <- effect != "Beta"
  tr <- function(v) if (logscale) log(pmax(as.numeric(v), 1e-8)) else as.numeric(v)
  lo <- tr(co$conf.low); hi <- tr(co$conf.high); est <- tr(co$estimate)
  rng <- suppressWarnings(range(c(lo, hi, est, 0), na.rm = TRUE, finite = TRUE))
  if (length(rng) != 2L || !all(is.finite(rng)) || diff(rng) == 0)
    rng <- c(-1, 1)
  pad <- diff(rng) * 0.08
  rng <- rng + c(-pad, pad)
  pos <- function(v) {
    p <- 100 * (v - rng[1]) / diff(rng)
    p[!is.finite(p)] <- 50
    pmax(0, pmin(100, p))
  }
  ref <- pos(0)
  rows <- vapply(seq_len(nrow(co)), function(i) {
    l <- pos(lo[i]); h <- pos(hi[i]); e <- pos(est[i])
    sig <- isTRUE(co$significant[i])
    paste0('<div class="frow"><div class="flab">', .mc_esc(co$term[i]),
           '</div><div class="ftrack"><span class="fref" style="left:',
           sprintf("%.1f", ref), '%"></span><span class="fline',
           if (sig) " sig" else "", '" style="left:',
           sprintf("%.1f", min(l, h)), "%;width:", sprintf("%.1f", abs(h - l)),
           '%"></span><span class="fdot', if (sig) " sig" else "",
           '" style="left:', sprintf("%.1f", e), '%"></span></div>',
           '<div class="fval">', .mc_ci_fmt(co$estimate[i], co$conf.low[i],
                                            co$conf.high[i]), "</div></div>")
  }, character(1))
  lab <- if (logscale) sprintf("%.2f to %.2f", exp(rng[1]), exp(rng[2])) else
    sprintf("%.2f to %.2f", rng[1], rng[2])
  ref_lab <- if (logscale) "1 (no effect)" else "0 (no effect)"
  paste0('<div class="forest">', paste(rows, collapse = ""),
         '<div class="faxis"><span>', .mc_esc(lab), '</span><span class="refnote">',
         'vertical line = ', ref_lab, '</span></div></div>')
}

.mc_css <- paste0(
  "*{box-sizing:border-box}body{margin:0;background:#f5f7fb;color:#1c2536;",
  "font:14px/1.55 -apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif}",
  ".wrap{max-width:1180px;margin:0 auto;padding:38px 26px 70px}",
  "header{margin-bottom:26px}h1{font-size:26px;font-weight:600;letter-spacing:.02em;margin:0 0 6px}",
  ".sub{color:#7a869c;font-size:13px;margin:0}",
  ".kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:22px 0 26px}",
  ".kpi{background:#fff;border:1px solid #e6eaf1;border-radius:12px;padding:16px 18px;",
  "box-shadow:0 2px 10px rgba(28,45,84,.04)}",
  ".kpi b{display:block;font-size:22px;font-weight:600;color:#1a439c;letter-spacing:-.01em}",
  ".kpi span{display:block;font-size:12px;color:#8a94a6;margin-top:3px}",
  ".grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:16px}",
  "section.card{background:#fff;border:1px solid #e6eaf1;border-radius:14px;padding:20px 22px;",
  "margin-bottom:16px;box-shadow:0 2px 12px rgba(28,45,84,.04)}",
  "section.wide{grid-column:1/-1}h2{font-size:15px;font-weight:600;margin:0 0 14px;",
  "letter-spacing:.03em;color:#1a439c;border-left:3px solid #1a439c;padding-left:10px}",
  ".tw{overflow-x:auto}table{width:100%;border-collapse:collapse;font-size:12.5px}",
  "th{text-align:left;color:#6e7a90;font-weight:600;background:#f8fafd;padding:9px 10px;",
  "border-bottom:1px solid #e6eaf1;white-space:nowrap}",
  "td{padding:8px 10px;border-bottom:1px solid #f0f2f7;vertical-align:top}",
  "tr:last-child td{border-bottom:0}tbody tr:hover{background:#fafbfe}",
  "td.note{color:#8a94a6}.empty{color:#98a2b3;font-size:13px;margin:6px 0}",
  ".badge{display:inline-block;padding:2px 9px;border-radius:20px;font-size:11px;font-weight:600}",
  ".badge.ok{background:#e7f5ec;color:#1c7331}.badge.warn{background:#fdf0dc;color:#a8680c}",
  ".badge.muted{background:#eef0f4;color:#7a869c}",
  ".findings{margin:0;padding-left:18px}.findings li{margin-bottom:7px}",
  ".warn{background:#fff8ec;border:1px solid #f3dfb8;border-radius:10px;padding:12px 16px;",
  "margin-top:12px;font-size:13px;color:#8a5a0c}",
  ".forest{margin-top:6px}.frow{display:grid;grid-template-columns:190px 1fr 150px;gap:12px;",
  "align-items:center;padding:5px 0;font-size:12.5px}",
  ".flab{overflow:hidden;text-overflow:ellipsis;white-space:nowrap;color:#3b4759}",
  ".ftrack{position:relative;height:16px;background:#f4f6fa;border-radius:4px}",
  ".fref{position:absolute;top:0;bottom:0;width:1px;background:#b8c2d4}",
  ".fline{position:absolute;top:6px;height:4px;background:#9aa9c6;border-radius:3px}",
  ".fline.sig{background:#1a439c}.fdot{position:absolute;top:2px;width:11px;height:11px;",
  "margin-left:-5px;border-radius:50%;background:#fff;border:2px solid #7a869c}",
  ".fdot.sig{border-color:#1a439c;background:#1a439c}.fval{text-align:right;color:#5a6579;",
  "font-variant-numeric:tabular-nums}",
  ".faxis{display:flex;justify-content:space-between;color:#98a2b3;font-size:11px;",
  "margin:8px 0 0 202px}.refnote{color:#8a94a6}",
  "footer{margin-top:26px;color:#98a2b3;font-size:12px;border-top:1px solid #e6eaf1;padding-top:14px}",
  "@media(max-width:700px){.frow{grid-template-columns:120px 1fr 110px}.faxis{margin-left:132px}}",
  "th.sortable{cursor:pointer;user-select:none;white-space:nowrap}",
  "th.sortable:hover{color:#1a439c;text-decoration:underline}",
  "@media print{body{background:#fff}.card,.kpi{box-shadow:none;break-inside:avoid}",
  "header .sub,footer{color:#666}table{font-size:11px}.noprint{display:none}}")

.mc_card <- function(title, body, wide = FALSE)
  paste0('<section class="card', if (wide) " wide" else "", '"><h2>',
         .mc_esc(title), "</h2>", body, "</section>")

# ---- exported: one-call analysis -------------------------------------------

