.jwdid_group_rows <- function(object, type, window=NULL, cwindow=NULL,
                              select=rep(TRUE,nrow(object$data)), split=NULL) {
  d <- object$data
  keep <- d$.etr == 1 & select
  if (type == "simple") by <- NULL
  else if (type == "any") { keep <- d[[object$gvar]] != 0 & select; by <- NULL }
  else if (type == "attgt") {
    keep <- d[[object$gvar]] != 0 & select
    d$.attgt <- paste(d[[object$gvar]],d[[object$tvar]],sep="::")
    by <- ".attgt"
  }
  else if (type %in% c("group", "calendar")) {
    block <- if (is.null(object$ivar)) rep.int(1L, nrow(d)) else
      interaction(d[[object$gvar]], d[[object$ivar]], drop = TRUE)
    first_t <- ave(d[[object$tvar]], block, FUN = min)
    keep <- keep & first_t < d[[object$gvar]] & select
    by <- if (type == "group") object$gvar else object$tvar
  }
  else {
    keep <- if (object$type == "never") d$.tr != 0 else d$.etr == 1
    d$.event <- ifelse(d[[object$gvar]] == 0, NA,
                       d[[object$tvar]] - d[[object$gvar]])
    keep <- keep & !is.na(d$.event) & select
    if (!is.null(window)) keep <- keep & d$.event >= window[1] & d$.event <= window[2]
    if (!is.null(cwindow)) d$.event <- pmin(pmax(d$.event,cwindow[1]),cwindow[2])
    by <- ".event"
  }
  if (!is.null(split)) {
    d$.jwdid_by <- interaction(if(is.null(by))factor("ATT") else factor(d[[by]]),
                               factor(d[[split]]),drop=TRUE,sep="::")
    by <- ".jwdid_by"
  }
  list(data = d, rows = which(keep), by = by)
}

.jwdid_analytic <- function(object, dat, rows, by = NULL, agg_weights=NULL,
                            asis=FALSE) {
  fit <- object$model
  hi <- dat; lo <- dat; if (!asis) hi$.tr <- 1; lo$.tr <- 0
  X1 <- stats::model.matrix(fit, data = hi, type = "rhs")
  X0 <- stats::model.matrix(fit, data = lo, type = "rhs")
  cn <- names(stats::coef(fit))
  if (!identical(colnames(X1), cn) || !identical(colnames(X0), cn)) {
    stop("Treatment design columns do not match fitted coefficients.", call. = FALSE)
  }
  D <- X1[, cn, drop = FALSE] - X0[, cn, drop = FALSE]
  D <- D[rows, , drop = FALSE]
  dat <- dat[rows, , drop = FALSE]
  w <- if (!is.null(agg_weights)) agg_weights[rows] else
    if (is.null(object$weights)) rep(1, nrow(dat)) else dat[[object$weights]]
  groups <- if (is.null(by)) factor(rep("ATT", nrow(dat))) else factor(dat[[by]])
  lev <- levels(groups)
  nonlinear <- !is.null(fit$family) && fit$family$family != "gaussian"
  if (nonlinear) {
    eta1 <- stats::predict(fit,newdata=hi,type="link")[rows]
    eta0 <- stats::predict(fit,newdata=lo,type="link")[rows]
    individual <- fit$family$linkinv(eta1)-fit$family$linkinv(eta0)
    A <- fit$family$mu.eta(eta1)*X1[rows,cn,drop=FALSE] -
      fit$family$mu.eta(eta0)*X0[rows,cn,drop=FALSE]
  } else {
    individual <- as.vector(D %*% stats::coef(fit)); A <- D
  }
  G <- t(vapply(lev, function(z) {
    ii <- groups == z; colSums(A[ii, , drop = FALSE] * w[ii]) / sum(w[ii])
  }, numeric(ncol(A))))
  b <- vapply(lev,function(z) {
    ii <- groups == z; sum(individual[ii]*w[ii])/sum(w[ii])
  },numeric(1))
  V <- G %*% stats::vcov(fit) %*% t(G)
  se <- sqrt(pmax(diag(V), 0)); z <- ifelse(se == 0, NA_real_, b / se)
  tab <- data.frame(term = lev, estimate = b, std.error = se,
    conf.low = b - stats::qnorm(.975) * se,
    conf.high = b + stats::qnorm(.975) * se,
    statistic = z, p.value = 2 * stats::pnorm(-abs(z)), check.names = FALSE)
  list(table = tab, b = setNames(b, lev), V = V)
}

.jwdid_marginaleffects <- function(object, dat, rows, by = NULL, agg_weights=NULL) {
  cn <- names(stats::coef(object$model))
  mm <- stats::model.matrix(object$model, data = dat, type = "rhs")
  if (!identical(colnames(mm), cn)) {
    stop("Treatment design columns do not match fitted coefficients.", call. = FALSE)
  }
  dat <- dat[rows, , drop = FALSE]
  nonlinear <- !is.null(object$model$family) && object$model$family$family != "gaussian"
  absorbed_nonlinear <- nonlinear && object$method %in% c("fepois","ppmlhdfe")
  args <- list(model = object$model, variables = list(.tr = c(0, 1)),
               newdata = dat, vcov = if (absorbed_nonlinear) FALSE else stats::vcov(object$model))
  if (!is.null(agg_weights)) {
    dat$.jwdid_agg_w <- agg_weights[rows]; args$newdata <- dat; args$wts <- ".jwdid_agg_w"
  } else if (!is.null(object$weights)) args$wts <- object$weights
  if (!is.null(by)) args$by <- by
  z <- do.call(marginaleffects::avg_comparisons, args)
  tab <- as.data.frame(z)
  term <- if (is.null(by)) "ATT" else as.character(tab[[by]])
  keep <- c("estimate", "std.error", "conf.low", "conf.high", "statistic", "p.value")
  for (nm in setdiff(keep,names(tab))) tab[[nm]] <- NA_real_
  out <- data.frame(term = term, tab[, keep, drop = FALSE], check.names = FALSE)
  V <- if (absorbed_nonlinear) matrix(NA_real_,nrow(out),nrow(out)) else stats::vcov(z)
  list(table = out, b = setNames(out$estimate, term), V = V)
}

#' Aggregate treatment effects
#'
#' @param object A fitted `jwdid` object.
#' @param type Aggregation type.
#' @param weights Optional aggregation weights.
#' @param orestriction Optional logical restriction evaluated in the model frame.
#' @param over,over2 Optional grouping variables.
#' @param window,cwindow Event-time selection or censoring endpoints.
#' @param pretrend Compute the joint pre-trend Wald test.
#' @param asis Use observed continuous treatment intensity.
#' @param engine Contrast engine.
#' @param ... Reserved; unknown arguments are rejected.
#' @return A `jwdid_aggte` object.
#' @details
#' The supported aggregation types are:
#'
#' * `"simple"`: one average treatment effect on the treated (ATT);
#' * `"group"`: ATT by first-treatment cohort;
#' * `"calendar"`: ATT by calendar period;
#' * `"event"`: ATT by event time;
#' * `"attgt"`: cohort-by-calendar-time ATT; and
#' * `"any"`: an average over all treated observations.
#'
#' `window` drops event times outside two inclusive endpoints, whereas
#' `cwindow` censors them into the endpoint bins. `pretrend = TRUE` appends a
#' joint Wald test of available event coefficients before event time -1 and
#' requires a model fitted with `never = TRUE`. `orestriction` restricts the
#' aggregation sample. `over` splits `"simple"` or `"any"` estimates and
#' `over2` splits `"group"`, `"calendar"`, or `"event"` estimates.
#'
#' `engine = "analytic"` obtains standard errors by the delta method directly
#' from the fitted model. `"marginaleffects"` uses
#' [marginaleffects::avg_comparisons()], and `"auto"` selects an appropriate
#' engine for the fitted family.
#'
#' @examplesIf requireNamespace("did", quietly = TRUE)
#' data("mpdta", package = "did")
#' fit <- jwdid(
#'   lemp ~ 1, mpdta,
#'   ivar = countyreal, tvar = year, gvar = first.treat,
#'   never = TRUE
#' )
#' jwdid_simple(fit)
#' event <- jwdid_event(fit, window = c(-3, 3), pretrend = TRUE)
#' event
#' \dontrun{
#' plot(event)
#' }
#' @references
#' Wooldridge, J. M. (2021). Two-Way Fixed Effects, the Two-Way Mundlak
#' Regression, and Difference-in-Differences Estimators.
#' \doi{10.2139/ssrn.3906345}.
#'
#' Callaway, B. and Sant'Anna, P. H. C. (2021).
#' Difference-in-Differences with multiple time periods.
#' *Journal of Econometrics*, 225, 200-230.
#' \doi{10.1016/j.jeconom.2020.12.001}.
#' @seealso [jwdid()], [plot.jwdid_aggte()]
#' @export
jwdid_aggte <- function(object, type = c("simple", "group", "calendar", "event"),
                  weights = NULL, orestriction = NULL, over = NULL, over2 = NULL,
                  window = NULL, cwindow = NULL, pretrend = FALSE,
                  asis = FALSE,
                  engine = c("auto", "analytic", "marginaleffects"), ...) {
  raw_names <- names(as.list(sys.call()))[-1L]
  raw_names <- raw_names[nzchar(raw_names)]
  exact <- setdiff(names(formals(jwdid_aggte)),"...")
  partial <- setdiff(raw_names,exact)
  if (length(partial)) stop("Unknown argument(s): ",paste(partial,collapse=", "),call.=FALSE)
  if (!inherits(object, "jwdid")) stop("`object` must be a jwdid model.", call. = FALSE)
  dots <- list(...)
  if (length(dots)) stop("Unknown argument(s): ",paste(names(dots),collapse=", "),call.=FALSE)
  type <- match.arg(type,c("simple","group","calendar","event","attgt","any")); engine <- match.arg(engine)
  if (!is.null(window) && !is.null(cwindow)) stop("`window` and `cwindow` are mutually exclusive.",call.=FALSE)
  for (z in list(window,cwindow)) if (!is.null(z) && (length(z)!=2L || z[1]>=z[2]))
    stop("Event windows must contain two increasing endpoints.",call.=FALSE)
  if ((!is.null(window)||!is.null(cwindow)||pretrend) && type!="event")
    stop("Event windows and pretrend are only available for event aggregation.",call.=FALSE)
  if (!is.null(over) && !type %in% c("simple","any")) stop("`over` is only for simple or any.",call.=FALSE)
  if (!is.null(over2) && !type %in% c("group","calendar","event")) stop("`over2` is only for group, calendar, or event.",call.=FALSE)
  if (!is.null(over) && !is.null(over2)) stop("Specify only one of `over` and `over2`.",call.=FALSE)
  d <- object$data
  select <- rep(TRUE,nrow(d))
  if (!is.null(orestriction)) {
    expr <- if(inherits(orestriction,"formula")) orestriction[[2L]] else orestriction
    select <- eval(expr,d,parent.frame())
    if (!is.logical(select)||length(select)!=nrow(d)) stop("`orestriction` must evaluate to one logical value per row.",call.=FALSE)
    select[is.na(select)] <- FALSE
  }
  split_fml <- if(!is.null(over)) over else over2
  split <- if(is.null(split_fml)) NULL else {
    nm <- if(inherits(split_fml,"formula")) all.vars(split_fml) else as.character(split_fml)
    if(length(nm)!=1L||!nm%in%names(d)) stop("Grouping variable must name one model-frame column.",call.=FALSE)
    nm
  }
  agg_w <- NULL
  if (!is.null(weights)) {
    if (inherits(weights,"formula")) {
      nm <- all.vars(weights); if(length(nm)!=1L) stop("`weights` formula must name one column.",call.=FALSE)
      agg_w <- d[[nm]]
    } else if(is.character(weights)&&length(weights)==1L) agg_w <- d[[weights]]
    else if(is.numeric(weights)&&length(weights)==nrow(d)) agg_w <- weights
    else stop("Invalid aggregation `weights`.",call.=FALSE)
  }
  nonlinear <- !is.null(object$model$family) && object$model$family$family != "gaussian"
  if (asis && is.null(object$trtvar)) warning("`asis` has no effect without `trtvar`.",call.=FALSE)
  if (asis) engine <- "analytic"
  if (engine == "marginaleffects" && nonlinear && object$method %in% c("fepois","ppmlhdfe")) {
    stop("marginaleffects cannot compute SEs with absorbed nonlinear fixed effects; use `engine = 'analytic'`.",call.=FALSE)
  }
  if (engine == "auto") engine <- if (nonlinear && !object$method %in% c("fepois","ppmlhdfe")) "marginaleffects" else "analytic"
  gr <- .jwdid_group_rows(object,type,window,cwindow,select,split)
  if (!length(gr$rows)) stop("The requested aggregation has no observations.", call. = FALSE)
  a <- if (engine == "analytic") .jwdid_analytic(object,gr$data,gr$rows,gr$by,agg_w,asis) else
    .jwdid_marginaleffects(object,gr$data,gr$rows,gr$by,agg_w)
  out <- c(as.list(a$table), list(type = type, b = a$b, V = a$V,
    call = match.call(), parent = object$call))
  class(out) <- "jwdid_aggte"
  if (pretrend) {
    if (object$type != "never") stop("`pretrend` requires `never = TRUE`.",call.=FALSE)
    ev <- suppressWarnings(as.numeric(out$term)); ii <- which(ev < -1)
    if (!length(ii)) stop("No pre-treatment event estimates are available.",call.=FALSE)
    VV <- out$V[ii,ii,drop=FALSE]; bb <- out$b[ii]
    ee <- eigen(VV,symmetric=TRUE); tol <- max(dim(VV))*max(abs(ee$values))*.Machine$double.eps
    pos <- ee$values > tol; rank <- sum(pos)
    Vinv <- if(rank) ee$vectors[,pos,drop=FALSE] %*% diag(1/ee$values[pos],rank) %*% t(ee$vectors[,pos,drop=FALSE]) else matrix(0,nrow(VV),ncol(VV))
    chi2 <- as.numeric(t(bb)%*%Vinv%*%bb)
    out$pretrend <- list(chi2=chi2,df=rank,p=stats::pchisq(chi2,rank,lower.tail=FALSE))
  } else out$pretrend <- NULL
  out
}

# Kept internal for source compatibility; not exported because did::aggte masks it.
aggte <- jwdid_aggte
#' @rdname jwdid_aggte
#' @export
estat <- jwdid_aggte
#' @rdname jwdid_aggte
#' @export
jwdid_simple <- function(object, ...) jwdid_aggte(object, "simple", ...)
#' @rdname jwdid_aggte
#' @export
jwdid_group <- function(object, ...) jwdid_aggte(object, "group", ...)
#' @rdname jwdid_aggte
#' @export
jwdid_calendar <- function(object, ...) jwdid_aggte(object, "calendar", ...)
#' @rdname jwdid_aggte
#' @export
jwdid_event <- function(object, ...) jwdid_aggte(object, "event", ...)
#' @rdname jwdid_aggte
#' @export
jwdid_attgt <- function(object, ...) jwdid_aggte(object, "attgt", ...)
#' @rdname jwdid_aggte
#' @export
jwdid_any <- function(object, ...) jwdid_aggte(object, "any", ...)
