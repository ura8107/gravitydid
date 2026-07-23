#' Estimate an extended two-way fixed-effects DiD model
#'
#' `anticipation = 0` corresponds to omitting Stata's `anticipation()` option.
#' R value `k` corresponds to Stata `anticipation(k + 1)`.
#' @param fml Model formula.
#' @param data Data frame.
#' @param ivar,tvar,gvar Unit, time, and first-treatment columns.
#' @param trtvar Optional continuous treatment column.
#' @param never Use never-treated rather than not-yet-treated controls.
#' @param anticipation Number of affected pre-treatment periods.
#' @param hettype Heterogeneity restriction.
#' @param group Use cohort rather than unit fixed effects.
#' @param hettype_ll,hettype_ul Event-time binning limits.
#' @param hettype_recode,hettype_evbase Event-cohort recode and base.
#' @param xattvar,exogvar,xtvar,xgvar,fevar One-sided role formulas.
#' @param xasis Use raw rather than cell-demeaned treatment covariates.
#' @param method Estimation backend.
#' @param corr Apply the unbalanced-panel correction.
#' @param cre Apply correlated-random-effects controls.
#' @param cluster,weights Cluster and estimation-weight columns.
#' @param ... Named arguments passed to the fixest backend.
#' @return A `jwdid` object.
#' @details
#' `hettype` controls which treatment-effect cells enter the regression:
#' `"timecohort"` (the default) allows cohort-by-time effects, `"cohort"` and
#' `"time"` impose the corresponding single-index restriction, `"event"`
#' indexes effects by event time, `"eventcohort"` allows event-by-cohort
#' effects, and `"twfe"` imposes one common effect. `hettype_ll` and
#' `hettype_ul` bin event times outside the supplied limits;
#' `hettype_recode` and `hettype_evbase` provide explicit event-time recoding.
#'
#' By default, not-yet-treated observations form the comparison group.
#' `never = TRUE` instead restricts comparisons to never-treated observations.
#' `anticipation = k` shifts treatment status back by `k` periods. Thus the R
#' default, `anticipation = 0`, corresponds to omitting Stata's
#' `anticipation()` option, while R value `k` corresponds to Stata value
#' `k + 1`.
#'
#' A missing `method` fits the linear ETWFE model with [fixest::feols()].
#' `"fepois"` and its Stata-compatible alias `"ppmlhdfe"` use
#' [fixest::fepois()]; `"poisson"`, `"logit"`, and `"probit"` use
#' [fixest::feglm()]. `corr` adds the unbalanced-panel correction when group
#' fixed effects are used. `cre` adds correlated-random-effects controls and is
#' available for the non-absorbed nonlinear models. Supplying `trtvar` enables
#' continuous treatment intensity; use `asis = TRUE` in [jwdid_aggte()] to
#' aggregate the observed intensity contrast.
#'
#' @examplesIf requireNamespace("did", quietly = TRUE)
#' data("mpdta", package = "did")
#' fit <- jwdid(
#'   lemp ~ 1, mpdta,
#'   ivar = countyreal, tvar = year, gvar = first.treat,
#'   never = TRUE
#' )
#' fit
#' jwdid_simple(fit)
#' @references
#' Wooldridge, J. M. (2021). Two-Way Fixed Effects, the Two-Way Mundlak
#' Regression, and Difference-in-Differences Estimators.
#' \doi{10.2139/ssrn.3906345}.
#'
#' Wooldridge, J. M. (2023). Simple approaches to nonlinear
#' difference-in-differences with panel data. *The Econometrics Journal*, 26,
#' C31-C66. \doi{10.1093/ectj/utad016}.
#' @seealso [jwdid_aggte()], [plot.jwdid_aggte()], [fixest::feols()]
#' @export
jwdid <- function(fml, data, ivar = NULL, tvar, gvar = NULL, trtvar = NULL,
                  never = FALSE, anticipation = 0,
                  hettype = "timecohort", group = FALSE,
                  hettype_ll = NULL, hettype_ul = NULL,
                  hettype_recode = NULL, hettype_evbase = NULL,
                  xattvar = NULL, exogvar = NULL, xtvar = NULL, xgvar = NULL,
                  fevar = NULL, xasis = FALSE,
                  method = NULL, corr = FALSE, cre = FALSE,
                  cluster = NULL, weights = NULL, ...) {
  data <- as.data.frame(data)
  caller <- parent.frame()
  tvar_nm <- .jwdid_name(substitute(tvar), data, "tvar", caller)
  gvar_nm <- if (missing(gvar) || is.null(substitute(gvar))) NULL else .jwdid_name(substitute(gvar), data, "gvar", caller)
  ivar_nm <- if (missing(ivar) || is.null(substitute(ivar))) NULL else .jwdid_name(substitute(ivar), data, "ivar", caller)
  trt_nm <- if (missing(trtvar) || is.null(substitute(trtvar))) NULL else .jwdid_name(substitute(trtvar), data, "trtvar", caller)
  if (is.null(gvar_nm) && is.null(trt_nm)) stop("Supply either `gvar` or `trtvar`.",call.=FALSE)
  if (is.null(gvar_nm)) {
    gvar_nm <- ".jwdid_gvar"
    data[[gvar_nm]] <- 0
  }
  cl_nm <- if (missing(cluster) || is.null(substitute(cluster))) ivar_nm else .jwdid_name(substitute(cluster), data, "cluster", caller)
  wt_nm <- if (missing(weights) || is.null(substitute(weights))) NULL else .jwdid_name(substitute(weights), data, "weights", caller)
  dots <- list(...)
  if (length(dots)) {
    dn <- names(dots)
    if (is.null(dn) || any(!nzchar(dn))) stop("All arguments in `...` must be named.", call.=FALSE)
    allowed <- setdiff(names(formals(fixest::feols)), "...")
    bad <- setdiff(dn, allowed)
    if (length(bad)) stop("Unknown argument(s): ", paste(bad,collapse=", "), call.=FALSE)
  }
  hettype <- match.arg(hettype, c("timecohort", "cohort", "time", "event",
                                  "eventcohort", "twfe"))
  if (!is.null(method)) method <- match.arg(method,c("ppmlhdfe","fepois","poisson","logit","probit"))
  if (cre && (is.null(method) || method %in% c("ppmlhdfe","fepois"))) {
    stop("`cre` is available with poisson, logit, or probit methods.",call.=FALSE)
  }
  group_eff <- group || is.null(ivar_nm) || (!is.null(method) && !method %in% c("ppmlhdfe","fepois"))
  if (hettype == "eventcohort") {
    if (!never) warning("`hettype = 'eventcohort'` forces `never = TRUE`.", call.=FALSE)
    if (!is.null(hettype_ll) || !is.null(hettype_ul)) {
      warning("`hettype_ll` and `hettype_ul` are ignored for eventcohort.", call.=FALSE)
      hettype_ll <- hettype_ul <- NULL
    }
    never <- TRUE
  }
  role_vars <- unique(c(all.vars(xattvar),all.vars(exogvar),all.vars(xtvar),
                        all.vars(xgvar),all.vars(fevar)))
  d <- .jwdid_prep(fml, data, ivar_nm, tvar_nm, gvar_nm, trt_nm,
                   anticipation, never, wt_nm, role_vars)
  pp <- attr(d, "prep")
  d <- .jwdid_cells(d, gvar_nm, tvar_nm, hettype, never, pp$antigap,
                    hettype_ll, hettype_ul, hettype_recode, hettype_evbase)
  x <- .jwdid_mm(.jwdid_rhs_formula(fml),d,".jwdid_x")
  d <- x$data; xvars <- x$names
  xa <- .jwdid_add_formula_vars(xattvar,d,".jwdid_xa"); d <- xa$data
  ex <- .jwdid_add_formula_vars(exogvar,d,".jwdid_ex"); d <- ex$data
  xt <- .jwdid_add_formula_vars(xtvar,d,".jwdid_xt"); d <- xt$data
  xg <- .jwdid_add_formula_vars(xgvar,d,".jwdid_xg"); d <- xg$data
  treat_raw <- c(xvars,xa$names)
  if (!xasis && length(treat_raw)) {
    dm <- .jwdid_demean(d,treat_raw,d$.cell_raw,wt_nm)
    d <- dm$data; treat_vars <- dm$names
  } else treat_vars <- treat_raw

  lhs <- deparse1(fml[[2]])
  tr_terms <- c("fixest::i(.cellf, .tr, ref = '.ref')",
    vapply(treat_vars,function(v) sprintf("fixest::i(.cellf, .tr * %s, ref = '.ref')",v),character(1)))
  varying <- .jwdid_is_varying(d,xvars,ivar_nm)
  main_x <- if (!is.null(ivar_nm) && !group_eff) xvars[varying] else xvars
  gref <- if (0 %in% d[[gvar_nm]]) 0 else pp$gvarmax
  ogvars <- c(if (!is.null(ivar_nm) && !group_eff) xvars[varying] else xvars,
              xg$names)
  og <- vapply(ogvars,function(v)sprintf("fixest::i(%s, %s, ref = %s)",gvar_nm,v,gref),character(1))
  otvars <- c(xvars,xt$names)
  ot <- vapply(otvars,function(v)sprintf("fixest::i(%s, %s, ref = %s)",tvar_nm,v,pp$tlist[1]),character(1))
  rhs <- paste(c(tr_terms,main_x,og,ot,ex$names),collapse=" + ")
  correction_names <- character()
  if (corr && group_eff && !is.null(ivar_nm)) {
    corr_rhs <- paste(rhs,sprintf("fixest::i(%s, ref = %s)",tvar_nm,pp$tlist[1]),sep=" + ")
    cc <- .jwdid_unit_means(d,corr_rhs,ivar_nm,wt_nm,prefix=".jwdid_corr")
    d <- cc$data; correction_names <- cc$names
    if (length(correction_names)) rhs <- paste(rhs,paste(correction_names,collapse=" + "),sep=" + ")
  }
  cre_names <- character()
  if (cre) {
    cre_abs <- if (!is.null(ivar_nm)) ivar_nm else gvar_nm
    cre_rhs <- paste(rhs,sprintf("fixest::i(%s, ref = %s)",tvar_nm,pp$tlist[1]),sep=" + ")
    cm <- .jwdid_unit_means(d,cre_rhs,cre_abs,wt_nm,prefix=".jwdid_cre",
                            drop_identical=TRUE)
    d <- cm$data; cre_names <- cm$names
  }
  fe <- c(if (!is.null(ivar_nm) && !group_eff) ivar_nm else gvar_nm, tvar_nm)
  if (!is.null(fevar)) fe <- c(fe,all.vars(fevar))
  explicit_nonlinear <- !is.null(method) && method %in% c("poisson","logit","probit")
  if (explicit_nonlinear) {
    if (!is.null(fevar)) warning("`fevar` is ignored by non-Poisson nonlinear methods.",call.=FALSE)
    rhs <- if (cre) paste(c(rhs,cre_names,
      sprintf("fixest::i(%s, ref = %s)",tvar_nm,pp$tlist[1])),collapse=" + ") else
      paste(c(rhs,sprintf("fixest::i(%s, ref = %s)",gvar_nm,gref),
              sprintf("fixest::i(%s, ref = %s)",tvar_nm,pp$tlist[1])),collapse=" + ")
    ff <- stats::as.formula(paste(lhs,"~",rhs))
  } else {
    ff <- stats::as.formula(paste(lhs, "~", rhs, "|", paste(fe, collapse = "+")))
  }
  vc <- if (!is.null(cl_nm)) stats::as.formula(paste0("~", cl_nm)) else "hetero"
  w <- if (!is.null(wt_nm)) stats::as.formula(paste0("~", wt_nm)) else NULL
  fit_args <- utils::modifyList(list(fml=ff,data=d,vcov=vc,weights=w,
                                      notes=FALSE,fixef.rm="none"),dots)
  if (is.null(method)) fit <- do.call(fixest::feols,fit_args)
  else if (method %in% c("ppmlhdfe","fepois")) fit <- do.call(fixest::fepois,fit_args)
  else {
    fit_args$family <- switch(method,poisson=stats::poisson(),
      logit=stats::binomial("logit"),probit=stats::binomial("probit"))
    fit <- do.call(fixest::feglm,fit_args)
  }
  balanced <- if (is.null(ivar_nm)) NA else {
    counts <- table(d[[ivar_nm]])
    length(unique(as.integer(counts))) == 1L && unique(as.integer(counts)) == length(unique(d[[tvar_nm]]))
  }
  out <- list(model = fit, call = match.call(), data = d, ivar = ivar_nm,
    tvar = tvar_nm, gvar = gvar_nm, hettype = hettype,
    type = if (never) "never" else "notyet", gap = pp$gap,
    anticipation = anticipation, antigap = pp$antigap, glist = pp$glist,
    tlist = pp$tlist, cluster = cl_nm, weights = wt_nm, xasis=xasis,
    method=if(is.null(method))"feols" else method, balanced=balanced,
    group=group_eff, corr=corr, cre=cre, trtvar=trt_nm,
    demean_info=list(raw=treat_raw,used=treat_vars,labels=c(x$labels,xa$labels)),
    correction_names=correction_names, cre_names=cre_names)
  class(out) <- "jwdid"
  out
}
