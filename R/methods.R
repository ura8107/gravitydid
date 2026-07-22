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
