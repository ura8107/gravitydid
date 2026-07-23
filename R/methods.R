#' @export
print.jwdid <- function(x, ...) { cat("jwdid ETWFE model\n"); print(x$model); invisible(x) }
#' @export
coef.jwdid <- function(object, ...) stats::coef(object$model, ...)
#' @export
vcov.jwdid <- function(object, ...) stats::vcov(object$model, ...)
#' @export
nobs.jwdid <- function(object, ...) stats::nobs(object$model, ...)
#' @export
predict.jwdid <- function(object, ...) stats::predict(object$model, ...)
#' @export
print.jwdid_aggte <- function(x, ...) {
  tab <- data.frame(term=x$term, estimate=x$estimate, std.error=x$std.error,
                    conf.low=x$conf.low, conf.high=x$conf.high)
  print(tab, row.names = FALSE); invisible(x)
}

#' @export
summary.jwdid <- function(object, ...) summary(object$model, ...)
#' @export
summary.jwdid_aggte <- function(object, ...) {
  out <- data.frame(term=object$term,estimate=object$estimate,
    std.error=object$std.error,statistic=object$statistic,p.value=object$p.value,
    conf.low=object$conf.low,conf.high=object$conf.high)
  class(out) <- c("summary.jwdid_aggte","data.frame"); out
}
#' @export
tidy.jwdid <- function(x, ...) {
  z <- stats::coef(x); se <- sqrt(diag(stats::vcov(x)))
  data.frame(term=names(z),estimate=unname(z),std.error=unname(se),
    statistic=unname(z/se),p.value=2*stats::pnorm(-abs(z/se)),row.names=NULL)
}
#' @export
tidy.jwdid_aggte <- function(x, ...) summary.jwdid_aggte(x)
#' @export
glance.jwdid <- function(x, ...) data.frame(nobs=stats::nobs(x),
  AIC=stats::AIC(x$model),BIC=stats::BIC(x$model),logLik=as.numeric(stats::logLik(x$model)))
#' @export
glance.jwdid_aggte <- function(x, ...) data.frame(type=x$type,n_terms=length(x$estimate),
  pretrend_p=if(is.null(x$pretrend))NA_real_ else x$pretrend$p)
