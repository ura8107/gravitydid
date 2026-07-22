.jwdid_prep <- function(fml, data, ivar, tvar, gvar, trtvar = NULL,
                        anticipation = 0, never = FALSE, weights = NULL,
                        extra_cols = NULL) {
  if (!is.numeric(anticipation) || length(anticipation) != 1L ||
      anticipation < 0 || anticipation != as.integer(anticipation)) {
    stop("`anticipation` must be one non-negative integer.", call. = FALSE)
  }
  keep <- .jwdid_complete(fml, data, c(ivar, tvar, gvar, trtvar,extra_cols), weights)
  original_row <- seq_len(nrow(data))[keep]
  d <- data[keep, , drop = FALSE]
  if (!nrow(d)) stop("No complete observations remain.", call. = FALSE)

  if (!is.null(trtvar)) {
    z <- d[[trtvar]]
    if (any(z < 0 | z > 1)) stop("`trtvar` must lie in [0, 1].", call. = FALSE)
    ids <- d[[ivar]]
    tt <- d[[tvar]]
    first <- ave(seq_along(z), ids, FUN = function(ii) {
      hit <- ii[z[ii] > 0]
      rep(if (length(hit)) min(tt[hit]) else 0, length(ii))
    })
    d[[gvar]] <- first
  }
  # The ado first removes units which are always treated under the raw gvar,
  # before detecting the no-never-treated fallback and computing the gap.
  block <- if (is.null(ivar)) rep.int(1L, nrow(d)) else d[[ivar]]
  first_t <- ave(d[[tvar]], block, FUN = min)
  always0 <- d[[gvar]] != 0 & first_t >= d[[gvar]] & d[[tvar]] >= d[[gvar]]
  original_row <- original_row[!always0]
  d <- d[!always0, , drop = FALSE]
  has_never <- any(d[[gvar]] == 0)
  gvarmax <- Inf
  if (!has_never) {
    gvarmax <- max(d[[gvar]])
    trim <- d[[tvar]] < gvarmax
    original_row <- original_row[trim]
    d <- d[trim, , drop = FALSE]
  }
  gap <- .jwdid_gap(d[[gvar]], d[[tvar]])
  anti <- anticipation + 1L
  antigap <- gap * anti
  antigap0 <- gap * (anti - 1L)

  block <- if (is.null(ivar)) rep.int(1L, nrow(d)) else d[[ivar]]
  first_t <- ave(d[[tvar]], block, FUN = min)
  always <- d[[gvar]] != 0 & first_t >= d[[gvar]] - antigap0 &
    d[[tvar]] >= d[[gvar]] - antigap0
  original_row <- original_row[!always]
  d <- d[!always, , drop = FALSE]
  g <- d[[gvar]]; tt <- d[[tvar]]
  if (is.null(trtvar)) {
    tr <- as.numeric(tt >= g - antigap & g > 0)
    if (never) tr[] <- 1
  } else {
    tr <- d[[trtvar]]
    if (never) tr[tr == 0 & g != 0] <- 1
    tr[tt < g - antigap] <- 0
  }
  tr[g >= gvarmax] <- 0
  d$.tr <- tr
  d$.etr <- as.numeric(tt > g - antigap & g > 0)
  d$.jwdid_row <- original_row
  attr(d, "prep") <- list(gap = gap, antigap = antigap,
    gvarmax = gvarmax, glist = sort(unique(g[g > 0 & g < gvarmax])),
    tlist = sort(unique(tt)), has_never = has_never)
  d
}
