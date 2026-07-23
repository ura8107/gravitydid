#' Methods for fitted jwdid models
#'
#' @description
#' Standard extraction, prediction, summary, and printing methods for objects
#' returned by [jwdid()].
#'
#' @param x,object A fitted `jwdid` object.
#' @param ... Additional arguments passed to the corresponding method for the
#'   underlying `fixest` model.
#' @return `print()` returns `x` invisibly; `coef()` and `vcov()` return model
#'   coefficients and their covariance matrix; `nobs()` returns the number of
#'   observations; `predict()` returns fitted-model predictions; and
#'   `summary()` returns the underlying `fixest` summary.
#' @name jwdid-methods
NULL

#' @rdname jwdid-methods
#' @export
print.jwdid <- function(x, ...) { cat("jwdid ETWFE model\n"); print(x$model); invisible(x) }
#' @rdname jwdid-methods
#' @export
coef.jwdid <- function(object, ...) stats::coef(object$model, ...)
#' @rdname jwdid-methods
#' @export
vcov.jwdid <- function(object, ...) stats::vcov(object$model, ...)
#' @rdname jwdid-methods
#' @export
nobs.jwdid <- function(object, ...) stats::nobs(object$model, ...)
#' @rdname jwdid-methods
#' @export
predict.jwdid <- function(object, ...) stats::predict(object$model, ...)

#' Methods for aggregated jwdid estimates
#'
#' @description
#' Printing and summary methods for objects returned by [jwdid_aggte()].
#'
#' @param x,object A `jwdid_aggte` object.
#' @param ... Additional arguments, currently unused.
#' @return `print()` returns `x` invisibly. `summary()` returns a data frame of
#'   estimates, standard errors, test statistics, p-values, and confidence
#'   intervals.
#' @name jwdid-aggte-methods
NULL

#' @rdname jwdid-aggte-methods
#' @export
print.jwdid_aggte <- function(x, ...) {
  tab <- data.frame(term=x$term, estimate=x$estimate, std.error=x$std.error,
                    conf.low=x$conf.low, conf.high=x$conf.high)
  print(tab, row.names = FALSE); invisible(x)
}

#' @rdname jwdid-methods
#' @export
summary.jwdid <- function(object, ...) summary(object$model, ...)
#' @rdname jwdid-aggte-methods
#' @export
summary.jwdid_aggte <- function(object, ...) {
  out <- data.frame(term=object$term,estimate=object$estimate,
    std.error=object$std.error,statistic=object$statistic,p.value=object$p.value,
    conf.low=object$conf.low,conf.high=object$conf.high)
  class(out) <- c("summary.jwdid_aggte","data.frame"); out
}
#' Tidy jwdid results
#'
#' @param x A `jwdid` or `jwdid_aggte` object.
#' @param ... Additional arguments, currently unused.
#' @return A data frame containing one row per coefficient or aggregate.
#' @name tidy-jwdid
NULL

#' @rdname tidy-jwdid
#' @export
tidy.jwdid <- function(x, ...) {
  z <- stats::coef(x); se <- sqrt(diag(stats::vcov(x)))
  data.frame(term=names(z),estimate=unname(z),std.error=unname(se),
    statistic=unname(z/se),p.value=2*stats::pnorm(-abs(z/se)),row.names=NULL)
}
#' @rdname tidy-jwdid
#' @export
tidy.jwdid_aggte <- function(x, ...) summary.jwdid_aggte(x)

#' Model-level summaries for jwdid results
#'
#' @param x A `jwdid` or `jwdid_aggte` object.
#' @param ... Additional arguments, currently unused.
#' @return A one-row data frame of model or aggregation diagnostics.
#' @name glance-jwdid
NULL

#' @rdname glance-jwdid
#' @export
glance.jwdid <- function(x, ...) data.frame(nobs=stats::nobs(x),
  AIC=stats::AIC(x$model),BIC=stats::BIC(x$model),logLik=as.numeric(stats::logLik(x$model)))
#' @rdname glance-jwdid
#' @export
glance.jwdid_aggte <- function(x, ...) data.frame(type=x$type,n_terms=length(x$estimate),
  pretrend_p=if(is.null(x$pretrend))NA_real_ else x$pretrend$p)
