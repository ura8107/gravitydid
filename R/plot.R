#' Plot aggregated treatment effects
#'
#' @param x A `jwdid_aggte` object.
#' @param style Confidence-interval style.
#' @param level Confidence level.
#' @param tight Replace sparse x values by consecutive positions.
#' @param ref_line Draw a horizontal zero line.
#' @param pre,post Lists of graphical parameters for event periods.
#' @param ... Additional arguments (currently unused).
#' @return A `ggplot2` object.
#' @details
#' `"errorbar"` draws points with vertical intervals, `"ribbon"` draws a line
#' and shaded interval, `"pointrange"` uses a compact point-range geometry, and
#' `"bar"` draws bars with interval whiskers. For event studies, `pre` and
#' `post` can each contain a `colour` value. `tight = TRUE` places sparse
#' periods at consecutive x positions, and `ref_line = FALSE` removes the
#' horizontal zero reference.
#'
#' @examplesIf requireNamespace("did", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
#' data("mpdta", package = "did")
#' fit <- jwdid(
#'   lemp ~ 1, mpdta,
#'   ivar = countyreal, tvar = year, gvar = first.treat,
#'   never = TRUE
#' )
#' event <- jwdid_event(fit, window = c(-3, 3))
#' plot(event, style = "pointrange")
#' @seealso [jwdid_aggte()], [ggplot2::ggplot()]
#' @export
plot.jwdid_aggte <- function(x,style=c("errorbar","ribbon","pointrange","bar"),
  level=.95,tight=FALSE,ref_line=TRUE,pre=list(),post=list(),...) {
  if (!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2 to plot results.",call.=FALSE)
  if (!x$type %in% c("group","calendar","event","attgt"))
    stop("Plotting is available for group, calendar, event, and attgt results.",call.=FALSE)
  style <- match.arg(style); d <- as.data.frame(tidy.jwdid_aggte(x))
  crit <- stats::qnorm(1-(1-level)/2); d$conf.low <- d$estimate-crit*d$std.error
  d$conf.high <- d$estimate+crit*d$std.error
  num <- suppressWarnings(as.numeric(d$term))
  d$.x <- if(all(is.finite(num))) num else seq_len(nrow(d))
  if (tight) d$.x <- seq_len(nrow(d))
  d$.period <- if(x$type=="event" && all(is.finite(num))) ifelse(num<0,"pre","post") else "estimate"
  p <- ggplot2::ggplot(d,ggplot2::aes(x=.x,y=estimate,colour=.period,fill=.period))
  if (ref_line) p <- p + ggplot2::geom_hline(yintercept=0,linetype=2,colour="grey50")
  p <- switch(style,
    errorbar=p+ggplot2::geom_errorbar(ggplot2::aes(ymin=conf.low,ymax=conf.high),width=.12)+ggplot2::geom_point(),
    ribbon=p+ggplot2::geom_ribbon(ggplot2::aes(ymin=conf.low,ymax=conf.high),alpha=.2,colour=NA)+ggplot2::geom_line()+ggplot2::geom_point(),
    pointrange=p+ggplot2::geom_pointrange(ggplot2::aes(ymin=conf.low,ymax=conf.high)),
    bar=p+ggplot2::geom_col(position="identity",alpha=.65)+ggplot2::geom_errorbar(ggplot2::aes(ymin=conf.low,ymax=conf.high),width=.12))
  if (x$type=="event" && any(num==-1))
    p <- p+ggplot2::geom_point(data=d[num==-1,,drop=FALSE],shape=21,fill="white",size=3,show.legend=FALSE)
  if (x$type=="event" && (!is.null(pre$colour) || !is.null(post$colour))) {
    cols <- c(pre=if(is.null(pre$colour))"#0072B2" else pre$colour,
              post=if(is.null(post$colour))"#D55E00" else post$colour)
    p <- p+ggplot2::scale_colour_manual(values=cols)+ggplot2::scale_fill_manual(values=cols)
  }
  p+ggplot2::labs(x=x$type,y="ATT",colour=NULL,fill=NULL)+ggplot2::theme_minimal()
}
