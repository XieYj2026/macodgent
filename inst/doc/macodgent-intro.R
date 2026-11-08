## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
library(macodgent)
data(sim_htn_trial)
d <- sim_htn_trial

## -----------------------------------------------------------------------------
dim(d)
attr(d, "simulated")
attr(d, "provenance")

## -----------------------------------------------------------------------------
mc_inspect(d[, c("id", "age", "bmi", "smoke", "crp")], detailed = TRUE)
mc_stat_mode(d$site)

## -----------------------------------------------------------------------------
mc_describe_numeric(d, vars = c("age", "bmi", "sbp_baseline"),
                    group = "treatment")
mc_test_normality(d$crp)
mc_correlate(d, vars = c("age", "bmi", "sbp_baseline", "ldl", "glucose"))

## -----------------------------------------------------------------------------
mc_compare_groups(d, group = "treatment",
                  vars = c("age", "bmi", "crp", "smoke", "sex"))
mc_baseline_table(d, vars = c("age", "bmi", "sex", "smoke", "sbp_baseline"),
                  group = "treatment")

## -----------------------------------------------------------------------------
sc <- mc_logistic_screen(d, y = "bp_control",
  x_vars = c("age", "bmi", "ldl", "glucose", "smoke", "diabetes"),
  p_adjust = "BH")
sc$univariate
rep <- mc_logistic_report(d,
  bp_control ~ age + bmi + treatment + diabetes, auc_boot = 0L)
rep$checks

## -----------------------------------------------------------------------------
km <- mc_kaplan_meier(d, time = "followup_days", status = "cv_event",
                      group = "treatment")
km$logrank
cox <- mc_cox(d, survival::Surv(followup_days, cv_event) ~
                age + bmi + treatment + diabetes)
cox$performance

## -----------------------------------------------------------------------------
rcs <- mc_rcs(d, x = "sbp_baseline", y = "sbp_week12", nk = 4)
rcs$nonlinearity_p
head(rcs$curve)

## -----------------------------------------------------------------------------
auto <- mc_auto_analyze(d, outcome = "bp_control", group = "treatment",
  vars = c("age", "bmi", "sbp_baseline", "ldl", "glucose", "smoke", "diabetes"),
  auc_boot = 0L)
auto$findings
auto$model$checks

