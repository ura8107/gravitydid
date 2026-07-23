
<!-- README.md is generated from README.Rmd. Please edit that file. -->

# jwdid

<!-- badges: start -->

[![R-CMD-check](https://github.com/ura8107/gravitydid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ura8107/gravitydid/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/license/mit)
<!-- badges: end -->

`jwdid` estimates extended two-way fixed-effects (ETWFE)
difference-in-differences models for staggered interventions. It brings
the released Stata `jwdid` workflow to R: flexible treatment-effect
heterogeneity, never- or not-yet-treated controls, anticipation,
nonlinear outcomes, continuous treatment intensity, and `estat`-style
aggregation.

## Installation

Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("ura8107/gravitydid")
```

## Quick start

The example uses the public `mpdta` county-level employment panel
distributed with the `did` package.

``` r
data("mpdta", package = "did")

fit <- jwdid(
  lemp ~ 1,
  data = mpdta,
  ivar = countyreal,
  tvar = year,
  gvar = first.treat,
  never = TRUE
)

jwdid_simple(fit)
#>  term    estimate  std.error    conf.low conf.high
#>   ATT -0.03995128 0.01179628 -0.06307155 -0.016831
```

Event-study aggregation returns one estimate per relative period:

``` r
event <- jwdid_event(fit, window = c(-3, 3), pretrend = TRUE)
event
#>  term    estimate  std.error     conf.low    conf.high
#>    -3  0.02502183 0.01815434 -0.010560032  0.060603691
#>    -2  0.02445874 0.01426679 -0.003503654  0.052421144
#>    -1  0.00000000 0.00000000  0.000000000  0.000000000
#>     0 -0.01993182 0.01185754 -0.043172166  0.003308533
#>     1 -0.05095737 0.01687068 -0.084023289 -0.017891445
#>     2 -0.13725874 0.03658948 -0.208972794 -0.065544684
#>     3 -0.10081136 0.03450427 -0.168438493 -0.033184233
event$pretrend
#> $chi2
#> [1] 3.017003
#>
#> $df
#> [1] 2
#>
#> $p
#> [1] 0.2212413
```

``` r
plot(
  event,
  style = "pointrange",
  pre = list(colour = "#0072B2"),
  post = list(colour = "#D55E00")
)
```

<img src="man/figures/README-event-plot-1.png" alt="" width="672" />

## Choosing an R DiD implementation

| Package | Primary approach | `jwdid`-style option parity | Nonlinear ETWFE | Aggregation |
|----|----|---:|---:|----|
| `jwdid` | Wooldridge ETWFE | Broad | Poisson, logit, probit | simple, cohort, calendar, event, ATT(g,t) |
| `etwfe` | Wooldridge ETWFE | Core estimator | Via package interface | marginal effects |
| `did` | Callaway–Sant’Anna group-time ATT | Different interface | No | group, calendar, event |

The packages answer related questions with different estimands and
implementations. Matching option names does not by itself make estimates
from different methods numerically identical.

## Validation

The implementation is covered by regression, error-contract,
nonlinear-model, aggregation, plotting, and Stata-parity tests. CI runs
`R CMD check` on Linux, macOS, and Windows. Numerical validation is
reported only for fixtures and reference outputs included in the test
suite; it is not a claim that every possible Stata model specification
is identical.

## References

- Wooldridge, J. M. (2021). *Two-Way Fixed Effects, the Two-Way Mundlak
  Regression, and Difference-in-Differences Estimators*.
  <https://doi.org/10.2139/ssrn.3906345>
- Wooldridge, J. M. (2023). Simple approaches to nonlinear
  difference-in-differences with panel data. *The Econometrics Journal*,
  26, C31–C66. <https://doi.org/10.1093/ectj/utad016>
- Callaway, B. and Sant’Anna, P. H. C. (2021). Difference-in-Differences
  with multiple time periods. *Journal of Econometrics*, 225, 200–230.
  <https://doi.org/10.1016/j.jeconom.2020.12.001>
