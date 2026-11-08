# =============================================================================
# Software name : macodgent medical data-science research toolkit
# Copyright holder : Dr.Xie
# Version : 1.0.2
# Language : R
# -----------------------------------------------------------------------------
# Module : backward-compatible names
# The 15 ambiguous names below were superseded by self-explanatory names.
# They remain exported so existing scripts keep working. No deprecation warning
# is emitted on purpose: these names are still used throughout the shipped
# documentation and teaching material, and silent forwarding avoids flooding
# every example and test with warnings. Migration is documented in NEWS.md.
# This file is named zzz-aliases.R so it is sourced last and can safely bind
# the new function objects directly (identical signatures, zero overhead).
# =============================================================================

#' Superseded function names (backward compatibility)
#'
#' These names are kept as exact aliases of the current, self-explanatory
#' function names. They behave identically and accept exactly the same
#' arguments. New code should use the replacement listed below.
#'
#' \tabular{ll}{
#' \strong{Superseded name} \tab \strong{Use instead} \cr
#' \code{mc_mode()} \tab \code{mc_stat_mode()} \cr
#' \code{mc_pcr()} \tab \code{mc_pca_regression()} \cr
#' \code{mc_desc_quant()} \tab \code{mc_describe_numeric()} \cr
#' \code{mc_norm_test()} \tab \code{mc_test_normality()} \cr
#' \code{mc_var_test()} \tab \code{mc_test_equal_variance()} \cr
#' \code{mc_nonparam()} \tab \code{mc_test_rank()} \cr
#' \code{mc_anova()} \tab \code{mc_anova_oneway()} \cr
#' \code{mc_survfit()} \tab \code{mc_kaplan_meier()} \cr
#' \code{mc_wide_long()} \tab \code{mc_pivot_longer()} \cr
#' \code{mc_long_wide()} \tab \code{mc_pivot_wider()} \cr
#' \code{mc_scale()} \tab \code{mc_scale_numeric()} \cr
#' \code{mc_tableone()} \tab \code{mc_baseline_table()} \cr
#' \code{mc_ph_test()} \tab \code{mc_test_proportional_hazards()} \cr
#' \code{mc_corr()} \tab \code{mc_correlate()} \cr
#' \code{mc_contingency()} \tab \code{mc_crosstab_test()} \cr
#' }
#'
#' @details
#' The renaming removes four classes of ambiguity: \code{mc_pcr} read as
#' polymerase chain reaction in a medical context; \code{mc_mode} collided with
#' base \code{mode()} (storage type); \code{mc_desc_quant}, \code{mc_norm_test}
#' and \code{mc_var_test} used abbreviations that could mean quantile,
#' normalisation or "variable"; \code{mc_survfit} duplicated the name of
#' \code{survival::survfit} while returning a different object.
#'
#' @name mc_legacy_aliases
#' @keywords internal
#' @examples
#' # SIMULATED DATA -- 3000 synthetic hypertension-trial subjects.
#' # Not real patients; see ?sim_htn_trial before interpreting any result.
#' data(sim_htn_trial)
#' d <- sim_htn_trial
#' mc_mode(d$site)                      # identical to mc_stat_mode(d$site)
#' mc_norm_test(d$crp)                  # identical to mc_test_normality(d$crp)
NULL

#' @rdname mc_legacy_aliases
#' @export
mc_mode <- mc_stat_mode

#' @rdname mc_legacy_aliases
#' @export
mc_pcr <- mc_pca_regression

#' @rdname mc_legacy_aliases
#' @export
mc_desc_quant <- mc_describe_numeric

#' @rdname mc_legacy_aliases
#' @export
mc_norm_test <- mc_test_normality

#' @rdname mc_legacy_aliases
#' @export
mc_var_test <- mc_test_equal_variance

#' @rdname mc_legacy_aliases
#' @export
mc_nonparam <- mc_test_rank

#' @rdname mc_legacy_aliases
#' @export
mc_anova <- mc_anova_oneway

#' @rdname mc_legacy_aliases
#' @export
mc_survfit <- mc_kaplan_meier

#' @rdname mc_legacy_aliases
#' @export
mc_wide_long <- mc_pivot_longer

#' @rdname mc_legacy_aliases
#' @export
mc_long_wide <- mc_pivot_wider

#' @rdname mc_legacy_aliases
#' @export
mc_scale <- mc_scale_numeric

#' @rdname mc_legacy_aliases
#' @export
mc_tableone <- mc_baseline_table

#' @rdname mc_legacy_aliases
#' @export
mc_ph_test <- mc_test_proportional_hazards

#' @rdname mc_legacy_aliases
#' @export
mc_corr <- mc_correlate

#' @rdname mc_legacy_aliases
#' @export
mc_contingency <- mc_crosstab_test
