# `jwdid` for R — Design Specification

**Target:** an R package that ports the Stata command `jwdid` (Fernando Rios-Avila) — Wooldridge's Extended TWFE (ETWFE) estimator for staggered DiD — including its post-estimation aggregation command `jwdid_estat` and its plotting layer.

**Reference implementation:** https://github.com/friosavila/stpackages `jwdid/` (`jwdid.ado` v2.201, `jwdid_estat.ado` v2.2, `jwdid_plot.ado`)
**Reference documentation:** https://friosavila.github.io/app_metrics/app_metrics11.html
**Underlying method:** Wooldridge (2021, 2023); gravity extensions from Nagengast, Rios-Avila & Yotov (2024).

---

## 0. Why a new package (positioning)

`etwfe` (Grant McDermott, CRAN 0.6.2) already implements the core ETWFE idea in R:

```r
etwfe(fml, tvar, gvar, data, ivar, xvar, tref, gref, cgroup = c("notyet","never"), fe, family, ...)
emfx(object, type = c("simple","group","calendar","event"), ...)
```

It does **not** cover the following, all of which are `jwdid` features and all of which are the reason this package exists:

| Missing in `etwfe` | `jwdid` feature |
|---|---|
| Heterogeneity restrictions | `hettype(time | cohort | timecohort | event | eventcohort | twfe)` |
| Anticipation | `anticipation(#)` |
| Unbalanced-panel correction | `corr` |
| Correlated random effects / Mundlak | `cre` |
| Separated covariate roles | `exogvar()`, `xtvar()`, `xgvar()`, `xattvar()` |
| Extra high-dimensional FE | `fevar()` |
| Continuous / intensity treatment | `trtvar()` + `estat ..., asis` |
| ATT(g,t) grid and "any period" aggregation | `estat attgt`, `estat any` |
| Sub-population aggregation restrictions | `estat ..., orestriction()` |
| Aggregation split by a covariate | `estat ..., over()` / `over2()` |
| Event-window selection and censoring | `estat event, window()` / `cwindow()` |
| Formal parallel-trends test | `estat event, pretrend` |
| Covariate demeaning control | `xasis` |

So: **`jwdid` (R) = a faithful, feature-complete port**, not a re-derivation. Numerical parity with the Stata command is the primary acceptance criterion.

**Package name:** `jwdid`. (Not taken on CRAN. Keeps the Stata name so users transfer knowledge directly.)

---

## 1. Technology mapping

| Stata | R |
|---|---|
| `reghdfe` | `fixest::feols()` |
| `ppmlhdfe` | `fixest::fepois()` |
| `poisson` / `logit` / `probit` / `fracreg` | `fixest::feglm()` |
| `margins, at() contrast(atcontrast(r)) subpop() over()` | `marginaleffects::avg_comparisons()` |
| `hdfe y x, abs() gen()` | `fixest::demean()` |
| `test` (joint Wald) | manual Wald from aggregated `b`, `V` |
| `estat plot` | `ggplot2` via a `plot()` S3 method |

**Dependencies.** Imports: `fixest`, `marginaleffects`, `stats`, `utils`, `generics`. Suggests: `ggplot2`, `testthat (>= 3.0)`, `did`, `etwfe`, `knitr`, `rmarkdown`, `modelsummary`, `data.table`.
Keep the hard dependency set small — `fixest` + `marginaleffects` is the whole engine.

**Environment on this machine (verified):** R 4.5.3; `fixest`, `marginaleffects`, `etwfe`, `did`, `devtools`, `testthat`, `roxygen2`, `ggplot2`, `pkgdown`, `data.table` all installed.

---

## 2. The central mechanism — the `.tr` switch

This is the single most important thing to port correctly, and it is what makes the whole post-estimation layer work.

In `jwdid.ado`, **every** treatment term in the regression is multiplied by a data column `__tr__ ∈ {0,1}`:

```stata
local xvar `xvar' c.__tr__#i`i'.`gvar'#i`j'.`tvar'
local xvar `xvar' c.__tr__#i`i'.`gvar'#i`j'.`tvar'#c.(`xxvar')
```

Then `jwdid_estat` recovers the ATT by moving that one switch:

```stata
margins, subpop(if __etr__==1) at(__tr__=(0 1)) contrast(atcontrast(r)) over(...)
```

`at(__tr__=0)` zeroes **all** treatment terms at once — including the covariate-interacted ones — giving \(\hat Y(0)\), the never-treated counterfactual. `at(__tr__=1)` gives \(\hat Y(obs)\). The contrast is \(\widehat{ATT}_i\), and `over()`/`subpop()` produce the aggregations. Under a nonlinear link this correctly returns the ATT **on the response scale**, not the latent scale.

### Porting requirement (non-negotiable)

The treatment terms must stay **symbolic** in the model formula so that `marginaleffects` can recompute the design matrix with `.tr = 0`. **Do not pre-materialize the interaction columns into `data`** — if you do, setting `.tr <- 0` will not propagate to the derived columns and every ATT will silently come back as zero or garbage.

Recommended construction:

```r
# cellf: factor whose levels are exactly the retained heterogeneity cells,
#        with every excluded observation and the reference cell collapsed
#        into a single level ".ref"
fixest::i(cellf, .tr, ref = ".ref")
```

`fixest::i(f, var, ref=)` builds `var * 1{f == level}` for every level except `ref`, and is evaluated at `model.matrix` time, so `predict(m, newdata)` with `.tr = 0` correctly returns \(\hat Y(0)\).

**Verify this round-trip before building anything else.** A one-off smoke test:

```r
m  <- feols(y ~ i(cellf, .tr, ref = ".ref") | id + t, data = d)
d0 <- transform(d, .tr = 0)
stopifnot(!isTRUE(all.equal(predict(m, d), predict(m, d0))))  # must differ
```

If `i()` does not round-trip through `marginaleffects` for some model class, fall back to **plan B**: compute the ATT analytically. For observation \(i\) with treatment-term design row \(x_i\) and coefficients \(\theta\),

- linear: \(\widehat{ATT}_i = x_i'\theta\)
- nonlinear: \(\widehat{ATT}_i = H(\eta_i) - H(\eta_i - x_i'\theta)\)

and aggregate as \(AGGTE_r = \sum_i w_i R_i \widehat{ATT}_i / \sum_i w_i R_i\) with delta-method SEs from \(\partial AGGTE_r/\partial\beta\) and the model `vcov`. Plan B is more code but fully deterministic; keep it behind an internal switch so it can be used to cross-check plan A in tests.

---

## 3. Sample construction — must match `jwdid.ado` exactly

Order matters. Reproduce it literally.

1. **Complete cases.** `touse` = non-missing on `y`, all covariates, `ivar`, `tvar`, `gvar`; intersected with `subset`.

2. **`gvar` semantics.** `gvar = 0` ⇒ never treated (within the window). `gvar = g > 0` ⇒ treatment starts at period `g`. Treatment is absorbing. If `trtvar` (a 0/1 or intensity variable in \([0,1]\)) is supplied instead, derive `gvar` as in `_gjwgvar`: per unit, the minimum `tvar` at which `trtvar > 0`, and `0` for units never positive. Error if `trtvar` has values outside \([0,1]\).

3. **`gap`.** Take the union of the strictly-positive unique values of `gvar` and `tvar`, sort, and set `gap` = the minimum first difference. This is the period spacing (usually 1, but 5 for a quinquennial panel, etc.). All anticipation arithmetic is in units of `gap`.

4. **Anticipation.** Stata: `anti = 1` by default, `antigap = gap * anti`, `antigap0 = gap * (anti - 1)`.
   **R exposes `anticipation = 0` as the default** (number of pre-treatment periods assumed already affected — matching the `did` package convention), and internally sets `anti = anticipation + 1`.
   ⚠️ Document this mapping loudly: R `anticipation = 0` ↔ Stata no option; R `anticipation = 1` ↔ Stata `anticipation(2)`.

5. **Drop always-treated units.** Drop observations where the unit's first observed period is already at or past the (anticipation-adjusted) treatment date:
   ```
   touse = 0  if  min(tvar | unit) >= gvar - antigap0  &  gvar != 0  &  tvar >= gvar - antigap0
   ```
   Such units have no usable pre-period.

6. **No never-treated fallback.** If no observation has `gvar == 0`:
   set `gvarmax = max(gvar)`, drop all observations with `tvar >= gvarmax`, and treat the last cohort as the control group. Otherwise `gvarmax = Inf`.

7. **`.tr` (the switch).** With `trtvar` absent:
   ```
   .tr = 0
   .tr = 1   if  tvar >= gvar - antigap  &  gvar > 0
   .tr = 1   (everywhere)                             if never = TRUE
   .tr = 0   if  gvar >= gvarmax
   ```
   Setting `.tr = 1` globally under `never` is harmless and intended: never-treated units have `gvar == 0`, which is not in `glist`, so no interaction term ever matches them.

   With `trtvar` present:
   ```
   .tr = trtvar
   .tr = 1   if  never = TRUE  &  trtvar == 0  &  gvar != 0
   .tr = 0   if  gvar >= gvarmax
   .tr = 0   if  tvar <  gvar - antigap
   ```

8. **`.etr` (effectively treated — the aggregation subpopulation).**
   ```
   .etr = 1  if  tvar > gvar - antigap  &  gvar > 0
   ```
   Note the **strict** `>` here versus `>=` for `.tr`. With no anticipation this makes `.etr` equal to `tvar >= gvar`, while `.tr` also switches on at `t = g - 1`. That one-period offset is the reference-period handling and must be preserved verbatim.

9. **`glist`** = sorted unique `gvar` with `0 < gvar < gvarmax`.
   **`tlist`** = sorted unique `tvar` with `tvar >= min(tvar)` in the estimation sample.

---

## 4. Heterogeneity cells (`hettype`)

Define `.post` (used by `cohort`, `time`, `twfe`):

```
.post = 0  "Base"      (gvar == 0, i.e. never treated)
.post = 1  "Pre-Trt"   tvar <  gvar - antigap   & gvar > 0
.post = 2  "Post-Trt"  tvar >= gvar - antigap0  & gvar > 0
```

Define `.evnt` (relative event time; used by `event`, `eventcohort`):

```
.evnt = (tvar - gvar) * (gvar > 0) - gap * (gvar == 0)
```
i.e. never-treated are parked one `gap` below zero. Then apply optional binning: `.evnt <- pmin(.evnt, ul)` and (only when `never = TRUE`) `.evnt <- pmax(.evnt, ll)`. Finally shift so the minimum is 1 (`cev = 1 - min(.evnt)`); keep `cev` on the object so labels can be mapped back to true event time.

### Cell definition and retained terms per `hettype`

| `hettype` | cell | terms retained (`never = TRUE`) | terms retained (`never = FALSE`) |
|---|---|---|---|
| `timecohort` (default) | `gvar × tvar`, `g ∈ glist`, `t ∈ tlist` | all non-empty cells except `t == g - antigap` | non-empty cells with `t >= g - antigap0` |
| `cohort` | `gvar × .post`, `.post ∈ {1,2}` | both `.post` levels | `.post == 2` only |
| `time` | `tvar × .post` | both | `.post == 2` only |
| `event` | `.evnt` | all levels except `e == -antigap` | levels with `e > -antigap` |
| `eventcohort` | `gvar × .evnt` (**forces `never = TRUE`**) | all cells except `.evnt == evbase` | n/a |
| `twfe` | `.post` | both | `.post == 2` only |

Only **non-empty** cells get a term (Stata does `count if ...` before adding). Everything not retained — including all never-treated observations — is folded into the `.ref` level of `cellf`.

`eventcohort` additionally takes a recode map and a base event value (Stata's `parse_hettype_evco` with `evbase()` and `erecode()`). In R, expose `hettype_recode` (a named list or a function applied to `.evnt`) and `hettype_evbase`. Error if any post-recode value is negative.

`hettype_ll` / `hettype_ul` correspond to Stata's `hettype(event, ll() ul())`; error if `ll >= ul`.

⚠️ **Upstream quirk — `ll` is not gated on `never`.** `jwdid.ado:381` reads

```stata
if "never"!="" & "`rll'"!="" qui: replace __evnt__=`rll' if __evnt__<`rll' & __evnt__!=.
```

`"never"` is a *literal string*, not the macro `` "`never'" `` (contrast line 380, which correctly writes `` "`rul'" ``, and line 526, which correctly writes `` "`never'" ``). The condition is therefore always true, so shipped Stata applies `ll` under `notyet` as well. This is almost certainly a missing pair of backticks upstream, but since numerical parity is the acceptance criterion, **match the shipped behaviour**: apply `ll` regardless of `never`, document it as an upstream quirk, and pin it with a test so the choice stays deliberate. It is a one-line change if the intended semantics are ever wanted instead.

---

## 5. Covariates

### 5.1 Roles

| Argument | Enters as | Stata name |
|---|---|---|
| RHS of `fml` (`x`) | main effects, interacted with treatment cells **and** with time and cohort | `varlist` after `y` |
| `xattvar` | interacted with treatment cells only | `xattvar()` |
| `exogvar` | additively only, no interactions | `exogvar()` / `exovar()` |
| `xtvar` | interacted with `tvar` only | `xtvar()` |
| `xgvar` | interacted with `gvar` only | `xgvar()` |
| `fevar` | additional absorbed fixed effects (linear / Poisson only) | `fevar()` |

### 5.2 Demeaning (`xasis = FALSE`, the default)

Stata: `hdfe y x xattvar, abs(cell) gen(_x_)` — covariates are demeaned **within the heterogeneity cell** (`gvar × tvar` for `timecohort`, etc.), and `_x_*` replaces the raw covariates inside the treatment interactions. This is Wooldridge's \(\tilde x\): it makes the raw \(\theta_{g,t}\) coefficients directly readable as ATT(g,t).

R: `fixest::demean(X, f = cellf_raw, weights = w)` where `cellf_raw` is the *un-collapsed* cell factor (never-treated form their own cells). With `xasis = TRUE`, skip demeaning and use raw covariates — aggregated ATTs are unchanged, only the raw coefficients lose their direct interpretation.

Note: the dependent variable is passed to `hdfe` in Stata but its demeaned version is immediately dropped (`capture drop _x_`y''). Don't demean `y`.

### 5.3 Time-constant vs time-varying covariates (`is_x_fix`)

When unit FE are used (`ivar` supplied and `group = FALSE`) and covariates are present, split `x` into:

- `xvarcons` — constant within unit (collinear with unit FE ⇒ excluded from the additive and cohort-interacted terms)
- `xvarvar` — varies within unit

Implementation: for each expanded covariate column, check whether `max - min` within unit is zero across all units. Use the *expanded* design (factor levels / interactions expanded, omitted levels dropped) as Stata's `ms_fvstrip ..., expand dropomit` does. When group FE are used, `xvarvar = x` (no split).

### 5.4 Non-treatment interaction blocks

- **`ogxvar`** — cohort-specific slopes: `i(gvar_f, z, ref = 0)` for every `z ∈ {xvarvar, xgvar}`. Reference cohort is `0` (never treated). Omitted entirely if all covariates are time-constant and the estimator is `feols`/`fepois` with unit FE.
- **`otxvar`** — time-specific slopes: `i(tvar_f, z, ref = first(tlist))` for every `z ∈ {x, xtvar}`. Stata skips the first time level (`if cj > 1`) to avoid collinearity — mirror that with `ref =`.
- **`exogvar`** — plain additive terms.

---

## 6. Estimation backends

| Condition | R call | Absorbed FE |
|---|---|---|
| `method = NULL`, `ivar` given, `group = FALSE` | `fixest::feols()` | `ivar + tvar + fevar` |
| `method = NULL`, group FE | `fixest::feols()` | `gvar + tvar + fevar` |
| `method = "ppmlhdfe"` / `"fepois"` | `fixest::fepois()` | `ivar + tvar + fevar` |
| `method = "poisson"` | `fixest::feglm(family = poisson())` | `gvar + tvar` |
| `method = "logit"` | `fixest::feglm(family = binomial("logit"))` | `gvar + tvar` |
| `method = "probit"` | `fixest::feglm(family = binomial("probit"))` | `gvar + tvar` |
| user-supplied function | passed through | as given |

**Group FE is forced** whenever `ivar` is missing, or `method` is set to anything other than `ppmlhdfe`/`fepois` — exactly as in the `.ado`. Document why: with non-Poisson nonlinear models, unit FE create an incidental-parameters problem. For linear models on a **balanced** panel, group FE and unit FE give numerically identical treatment coefficients; on an **unbalanced** panel they differ, which is what `corr` is for.

**Nonlinear models: do NOT absorb the fixed effects.** Stata's generic branch (`jwdid.ado:682`) enters `i.gvar i.tvar` as explicit regressors. Absorbing them in `feglm` gives the identical MLE *for the point estimates*, but `marginaleffects` then refuses to compute standard errors at all:

> For this model type, `marginaleffects` cannot take into account the uncertainty in fixed-effects parameters. … Set `vcov=FALSE` to compute estimates without standard errors.

Since the whole product is inference on aggregated ATTs, absorbing here would ship a nonlinear branch with no inference. Enter the FE as explicit dummies instead — this reproduces Stata exactly and is cheap, because every nonlinear method except `ppmlhdfe` is already forced onto **group** FE, so the dummy count is `#cohorts + #periods`, not `#units`. Verified on `mpdta` with `method = "poisson"`: absorbed 12 coefficients / SE unavailable; dummies 20 coefficients / identical point estimate (diff `4.7e-14`) and a usable SE (`0.01225`).

**`fepois` / `ppmlhdfe` keep absorbed unit FE**, where dummies are infeasible. Do not fall back to `NA` standard errors — use the generalised analytic delta method, conditioning on the estimated fixed effects, which is what Stata does via `ppmlhdfe …, d` (note that `jwdid.ado:664-666` passes `d` precisely so that `margins` can predict). With `ATT_i = H(a_i + x_i'θ) − H(a_i)`:

```
J = colMeans( H'(η1) * X1  −  H'(η0) * X0 )      over the aggregation subpopulation
Var(AGGTE) = J' V J
```

`predict(fit, newdata = , type = "link")` returns the full linear predictor including the absorbed FE, so this is directly implementable. Verified on `mpdta`: `fepois` ATT `-26.49`, delta-method SE `14.21`.

Also check the balanced/unbalanced flag (Stata's `is_balanced` mata routine) so the package can warn when `corr` would matter.

**Clustering.** Default `cluster = ivar` when `ivar` is given; none (heteroskedasticity-robust) otherwise. Explicit `cluster` always wins. Stata falls back to `vce(robust)` when there is no cluster variable and no `ivar` — mirror that.

### `corr` — unbalanced-panel correction

Stata's `myhdmean`: regress each regressor on the unit FE and keep the fitted values as extra controls, dropping any that end up constant. R: `fixest::demean()` to get residuals, then `fitted = x - resid`; drop columns with zero variance and run `_rmcoll`-style collinearity pruning (`fixest` handles the latter, but drop constants explicitly first). Only meaningful with group FE + `ivar` present.

### `cre` — Mundlak / correlated random effects

Stata's `cre_jwdid`: for each expanded regressor, compute its (weighted) within-`ivar` mean, add it as `_cre_*`, drop those identical to the original, then prune collinear ones. R: group means via `stats::ave()` or a `data.table`/`collapse` group-mean, then a QR-based collinearity drop. Then fit with `_cre_*` plus `i(tvar)` and no unit FE.

---

## 7. Aggregation — `jwdid_estat` equivalents

### 7.1 API

```r
aggte(object,
      type         = c("simple", "group", "calendar", "event", "attgt", "any"),
      weights      = NULL,        # aggregation weights; default = estimation weights
      orestriction = NULL,        # expression evaluated in the model frame, e.g. ~ dx == 0
      over         = NULL,        # split `simple`/`any` by a variable
      over2        = NULL,        # split group/calendar/event by a variable
      window       = NULL,        # event only: keep e in [a, b]
      cwindow      = NULL,        # event only: censor e at a and b
      pretrend     = FALSE,       # event only: joint pre-trend Wald test
      asis         = FALSE,       # continuous treatment: use observed intensity
      predict      = c("response", "link"),
      ...)
```

Ship thin wrappers `jwdid_simple()`, `jwdid_group()`, `jwdid_calendar()`, `jwdid_event()`, `jwdid_attgt()`, `jwdid_any()`, plus an `estat()` alias so Stata users can type what they already know. Returns an object of class `jwdid_aggte`.

### 7.2 Subpopulation and grouping per type — copy from the `.ado`

| `type` | subpopulation | grouped by |
|---|---|---|
| `simple` | `.etr == 1 & tosel` | none (or `over`, restricted to non-zero/non-missing where `.etr == 1`) |
| `group` | `.etr == 1 & tosel` | `gvar`, defined only where `.etr == 1 & first_obs_period < gvar` |
| `calendar` | `.etr == 1 & tosel` | `tvar`, same `first_obs_period < gvar` restriction |
| `event` | `.etr == 1 & tosel` if `type = "notyet"`; `.tr != 0 & tosel` if `type = "never"` | `tvar - gvar`, defined only where `gvar != 0` |
| `attgt` | `tosel` only — **no `.etr` restriction** | `gvar × tvar` for `gvar != 0` (the full ATT(g,t) grid, pre-periods included) |
| `any` | `gvar != 0 & tosel` — **no `.etr` restriction** | none (or `over`) |

`tosel` is `TRUE` everywhere unless `orestriction` is supplied, in which case it is the evaluated expression. The `first_obs_period < gvar` guard for `group`/`calendar` excludes units with no pre-period; it uses the min of `tvar` within `gvar × ivar` on the estimation sample.

### 7.3 The contrast

```r
marginaleffects::avg_comparisons(
  object$model,
  variables = list(.tr = c(0, 1)),   # lo = 0, hi = 1
  newdata   = subpop_rows,
  by        = grouping_var,
  wts       = agg_weights,
  type      = predict_type
)
```

With `asis = TRUE` (continuous/intensity treatment), the "hi" arm must be the **observed** `.tr`, not `1`. Implement by building two `newdata` frames (`lo` with `.tr = 0`, `hi` unchanged) and differencing `avg_predictions()` with the delta method, or via `marginaleffects`' custom-comparison interface — whichever round-trips cleanly. Same for the `link` scale.

### 7.4 `pretrend`

Joint Wald test that all event-time estimates with \(e < -1\) are zero:

\[H_0:\ AGGTE_e = 0 \ \ \forall e < -1\]

Compute from the aggregated `b` and `V`: \(W = b_{pre}' V_{pre}^{-} b_{pre}\), \(\chi^2\) with df = rank of \(V_{pre}\) (use a generalized inverse — the pre-period block can be rank-deficient). Only valid when the model was fit with `never = TRUE`; error otherwise. Return `chi2`, `df`, `p` on the object.

Document the caveat from the source: this tests the *event-aggregated* pre-treatment ATTs, which is **not** the same as Callaway–Sant'Anna's test over all group/time-specific ATTs.

### 7.5 `window` / `cwindow`

Mutually exclusive; error if both given. `window = c(a, b)` keeps only event times in `[a, b]`. `cwindow = c(a, b)` censors: `e <- pmin(pmax(e, a), b)`, so endpoint bins pool everything beyond them. Require `a < b`, integers.

---

## 8. Plot method

```r
plot(x, style = c("ribbon", "errorbar", "pointrange", "bar"),
     level = 0.95, tight = FALSE, ref_line = TRUE, ...)
```

Returns a `ggplot`. Maps Stata's `rarea → ribbon`, `rspike → errorbar` (the Stata default), `rcap → pointrange`, `rbar → bar`. For `event`, style pre- and post-treatment periods separately (Stata's `pstyle1/color1/lwidth1/barwidth1` vs `...2`) — expose as `pre` and `post` list arguments rather than numbered scalars. `tight = TRUE` recodes the x-axis to consecutive integers so sparse cohort/time values don't leave gaps. Only defined for `group`, `calendar`, `event`, `attgt`.

---

## 9. Object model and S3 methods

**`jwdid`** (from `jwdid()`) — fields:
`model` (the `fixest` fit), `call`, `data` (model frame incl. `.tr`, `.etr`, `.post`/`.evnt`, `cellf`), `ivar`, `tvar`, `gvar`, `hettype`, `type` (`"never"` / `"notyet"`), `gap`, `anticipation`, `antigap`, `glist`, `tlist`, `cev`, `balanced`, `method`, `cluster`, `weights`, `xasis`, `cre`, `corr`, `demean_info`.

**`jwdid_aggte`** (from `aggte()`) — fields:
`estimate`, `std.error`, `conf.low`, `conf.high`, `statistic`, `p.value`, `term`, `type`, `b`, `V`, `pretrend` (list or `NULL`), `call`, `parent` (the `jwdid` call).

**Methods:** `print`, `summary`, `coef`, `vcov`, `nobs`, `predict`, `plot`, and `tidy`/`glance` registered against `generics` so `modelsummary` and `broom` work out of the box.

---

## 10. Validation strategy

This is the acceptance criterion, not an afterthought. Parity with Stata is the product.

1. **Golden reference from Stata.** Put a `.do` file and a CSV of expected results in `inst/stata/`. Generate reference output for `mpdta` across the full option grid (`never` × `hettype` × `method` × covariates × `anticipation`). Tests compare against the CSV with `tolerance = 1e-6`. This makes parity reproducible and auditable rather than a claim.

2. **`mpdta` baseline.** The Callaway–Sant'Anna county teen-employment panel, available as `did::mpdta`. The canonical example:
   ```r
   m <- jwdid(lemp ~ 1, data = mpdta, ivar = countyreal, tvar = year, gvar = first, never = TRUE)
   aggte(m, "simple"); aggte(m, "calendar"); aggte(m, "group"); aggte(m, "event")
   ```

3. **Cross-package equivalence.**
   - `never = TRUE`, no covariates ⇒ `aggte(m, "attgt")` must equal `did::att_gt()` exactly (this equivalence is stated in the source documentation, and is the strongest available check).
   - Overlapping feature set ⇒ must equal `etwfe::etwfe()` + `etwfe::emfx()`.

4. **Internal consistency.**
   - `hettype = "event"` vs `hettype = "timecohort"`. ⚠️ **These are different models and do not agree in general — do not test for equality.** `timecohort` estimates the unrestricted ATT(g,t) grid, and its event aggregation is an explicit observation-count-weighted mean of those cells. `hettype = "event"` instead *imposes* cross-cohort homogeneity within each event time, so OLS returns a precision-weighted combination determined by the design's leverage structure. The two weighting schemes coincide only in special cases. (Stata's documentation is explicit that `hettype()` "reduces the heterogeneity of the treatment effects" — it is a restriction, not a reparameterization.)
     The valid tests are:
     - **Single treated cohort** ⇒ `e = t - g` is a bijection with `t`, so the two are the *same* model and must agree to machine precision. Verified on a synthetic single-cohort panel: max diff `0`.
     - **Hand-built restricted model** ⇒ `feols(y ~ i(.ef, .tr, ref = ".ref") | ivar + tvar)` built directly from `.evnt` must reproduce `hettype = "event"`. Verified on `mpdta`: max diff `1.0e-17`.
     - **Aggregation is a no-op** ⇒ `aggte(m_event, "event")` must return the event model's own coefficients.
     - On multi-cohort data the two *will* differ (on `mpdta`, by up to ~0.009). Assert that they differ rather than that they match, so a future refactor that silently collapses one into the other is caught.
   - Balanced panel: group FE ≡ unit FE for linear models (treatment coefficients identical).
   - `xasis = TRUE` and `xasis = FALSE` must give identical *aggregated* ATTs (only raw coefficients differ). This is a sharp test of the demeaning code.
   - Plan A (`marginaleffects`) vs plan B (analytic delta method) must agree.

5. **Nonlinear.** Poisson results must match a hand-built `fepois()` fit with the same design matrix.

6. **Edge cases with explicit tests:** no never-treated units; unbalanced panel; single treated cohort; `gap > 1` (quinquennial panel); anticipation > 0; always-treated units present (must be dropped); missing values scattered through covariates; weights supplied.

---

## 11. Implementation phases

Each phase ends green — `devtools::test()` passes and the phase's tests are written *before* the code (the repo has `r-skills:tdd` available; use it).

| Phase | Scope | Done when |
|---|---|---|
| **1. Skeleton + core** | `DESCRIPTION`/`NAMESPACE`/roxygen; data prep (§3); `hettype = "timecohort"`; linear `feols`; `.tr` round-trip verified; `aggte()` for `simple`/`group`/`calendar`/`event` | `mpdta` + `never = TRUE` matches `did::att_gt` and `etwfe::emfx` |
| **2. Heterogeneity** | remaining `hettype` values, `.post`/`.evnt`, `ll`/`ul` binning, `anticipation`, `never`/`notyet`, no-never-treated fallback | event aggregation from `hettype = "event"` ≡ from `"timecohort"` |
| **3. Covariates** | demeaning, `xasis`, `exogvar`/`xtvar`/`xgvar`/`xattvar`, `ogxvar`/`otxvar`, `is_x_fix` split | `xasis` TRUE/FALSE give identical aggregated ATTs |
| **4. Nonlinear + panel corrections** | `fepois`/`feglm` backends, forced group FE, `corr`, `cre`, balanced-panel detection | Poisson matches hand-built `fepois`; `corr` closes the unbalanced gap |
| **5. Full estat** | `attgt`, `any`, `over`/`over2`, `orestriction`, `window`/`cwindow`, `pretrend`, aggregation weights, `asis` continuous treatment | full option grid matches the Stata golden CSV |
| **6. Polish** | `plot()`, all S3 + `tidy`/`glance`, README, vignettes (getting started; Stata→R migration table; gravity/trade application), pkgdown, `R CMD check --as-cran` clean | zero ERRORs/WARNINGs/NOTEs |

---

## 12. Repository layout

```
jwdid/
├── DESCRIPTION
├── NAMESPACE                     # roxygen-generated
├── LICENSE                       # match upstream jwdid's license
├── README.md
├── R/
│   ├── jwdid.R                   # main entry point, argument handling
│   ├── prep.R                    # §3 sample construction, gvar/trtvar, gap, .tr, .etr
│   ├── cells.R                   # §4 hettype -> cellf, .post, .evnt
│   ├── covariates.R              # §5 demeaning, is_x_fix, role assignment
│   ├── formula.R                 # assembles the fixest formula
│   ├── fit.R                     # §6 backend dispatch, corr, cre
│   ├── aggte.R                   # §7 all estat equivalents
│   ├── pretrend.R                # §7.4 Wald test
│   ├── plot.R                    # §8
│   ├── methods.R                 # §9 print/summary/tidy/glance/...
│   ├── utils.R
│   └── data.R                    # dataset documentation
├── man/                          # roxygen-generated
├── data/                         # mpdta copy (or reference did::mpdta from Suggests)
├── inst/stata/                   # §10.1 golden reference .do + expected-results CSV
├── tests/testthat/
│   ├── test-prep.R
│   ├── test-cells.R
│   ├── test-covariates.R
│   ├── test-equiv-did.R
│   ├── test-equiv-etwfe.R
│   ├── test-equiv-stata.R
│   ├── test-aggte.R
│   ├── test-nonlinear.R
│   └── test-edge-cases.R
└── vignettes/
    ├── jwdid.Rmd                 # getting started
    ├── stata-migration.Rmd       # option-by-option Stata -> R table
    └── gravity.Rmd               # trade / gravity application
```

---

## 13. Stata → R argument mapping (ship this as a vignette table)

| Stata | R | Note |
|---|---|---|
| `jwdid y x, ...` | `jwdid(y ~ x, data = , ...)` | `y ~ 1` for no covariates |
| `ivar(i)` | `ivar = i` | omit ⇒ repeated cross-section, group FE |
| `tvar(t)` / `time(t)` | `tvar = t` | |
| `gvar(g)` | `gvar = g` | `0` = never treated |
| `trtvar(d)` | `trtvar = d` | 0/1 or intensity in \([0,1]\) |
| `cluster(c)` | `cluster = c` | default `ivar` |
| `never` | `never = TRUE` | not-yet-treated excluded as controls |
| `group` | `group = TRUE` | force group FE |
| `method(poisson)` | `method = "poisson"` | |
| `hettype(event, ll(-5) ul(5))` | `hettype = "event", hettype_ll = -5, hettype_ul = 5` | |
| `anticipation(2)` | `anticipation = 1` | ⚠️ **off by one** — see §3.4 |
| `xasis` | `xasis = TRUE` | |
| `corr` | `corr = TRUE` | |
| `cre` | `cre = TRUE` | |
| `exogvar(z)` | `exogvar = ~ z` | |
| `xtvar(z)` / `xgvar(z)` / `xattvar(z)` | `xtvar = ~ z` etc. | |
| `fevar(f)` | `fevar = ~ f` | linear / Poisson only |
| `[pw = w]` | `weights = ~ w` | |
| `estat simple` | `aggte(m, "simple")` | |
| `estat group` / `calendar` / `event` / `attgt` / `any` | `aggte(m, "...")` | |
| `estat event, pretrend` | `aggte(m, "event", pretrend = TRUE)` | |
| `estat event, window(-3 3)` | `aggte(m, "event", window = c(-3, 3))` | |
| `estat ..., orestriction(dx==0)` | `aggte(m, ..., orestriction = ~ dx == 0)` | |
| `estat ..., over(v)` | `aggte(m, ..., over = ~ v)` | |
| `estat plot` | `plot(aggte(m, "event"))` | |
| `estat ..., post` / `estore()` / `esave()` | not ported — assign the returned object | R idiom |
| `vce(unconditional)` | not ported initially | note the limitation; Stata itself cannot do this after `reghdfe`/`ppmlhdfe` |

---

## 14. Known risks

1. **`i()` × `marginaleffects` round-trip.** The whole post-estimation layer rests on `.tr` propagating symbolically. Verify in phase 1 before anything else; keep the analytic delta-method path (plan B) as an implemented, tested fallback rather than a hypothetical one.
2. **`demean()` centering convention.** Stata's `hdfe` may re-add the grand mean where `fixest::demean()` does not. Irrelevant for aggregated ATTs, visible in raw coefficients. Pin down against Stata in phase 3.
3. **Collinear-term pruning.** `reghdfe` drops collinear terms with specific naming (`o.` prefix); `fixest` drops them too but reports differently. Coefficient *names* will not match Stata; coefficient *values* must. Test on values only.
4. **`feglm` with absorbed FE vs Stata's explicit dummies.** Same MLE in theory; confirm numerically for logit and probit.
5. **Empty-cell handling.** Stata skips cells with zero observations. Factor levels in R must be dropped, not carried as all-zero columns, or `fixest` will drop them itself and the coefficient vector will misalign with expectations.
6. **`gap > 1` panels.** All anticipation arithmetic is in `gap` units. Get this wrong and quinquennial-panel results are silently off. Explicit test required.

---

## 15. Out of scope for v0.1

- `vce(unconditional)` standard errors (Stata cannot do this after `reghdfe`/`ppmlhdfe` either)
- `index` option (Stata result-storage bookkeeping)
- `estat ..., post` / `estore()` / `esave()` (R users assign objects)
- Stata `.ster` file interop
