#' Extended two-way fixed-effects difference-in-differences
#'
#' @description
#' `jwdid` estimates cohort- and time-specific treatment effects for staggered
#' interventions using Wooldridge's extended two-way fixed-effects (ETWFE)
#' framework. It provides the Stata `jwdid` feature set in R, including
#' heterogeneity restrictions, nonlinear models, continuous treatment,
#' anticipation, panel corrections, and post-estimation aggregation.
#'
#' @details
#' The package complements `did`, which implements semiparametric group-time
#' estimators, and `etwfe`, which implements the core ETWFE estimator. `jwdid`
#' adds Stata-compatible sample construction and the complete `estat`-style
#' aggregation layer. See `vignette("jwdid")` for a workflow,
#' `vignette("stata-migration")` for option mappings, and
#' `vignette("gravity")` for Poisson flow models.
#'
#' @references
#' Wooldridge, J. M. (2021). Two-Way Fixed Effects, the Two-Way Mundlak
#' Regression, and Difference-in-Differences Estimators.
#' \doi{10.2139/ssrn.3906345}.
#'
#' Wooldridge, J. M. (2023). Simple approaches to nonlinear
#' difference-in-differences with panel data. *The Econometrics Journal*, 26,
#' C31-C66. \doi{10.1093/ectj/utad016}.
#'
#' Callaway, B. and Sant'Anna, P. H. C. (2021). Difference-in-Differences with
#' multiple time periods. *Journal of Econometrics*, 225, 200-230.
#' \doi{10.1016/j.jeconom.2020.12.001}.
#'
#' @seealso [jwdid()], [jwdid_aggte()], [plot.jwdid_aggte()]
#' @keywords internal
"_PACKAGE"
