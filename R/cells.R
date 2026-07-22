.jwdid_cells_timecohort <- function(d, gvar, tvar, never, antigap) {
  g <- d[[gvar]]; tt <- d[[tvar]]
  retained <- g %in% attr(d, "prep")$glist
  if (never) retained <- retained & tt != g - antigap
  else retained <- retained & tt >= g - antigap + attr(d, "prep")$gap
  raw <- interaction(g, tt, drop = TRUE, lex.order = TRUE)
  labels <- paste(g, tt, sep = "::")
  vals <- ifelse(retained, labels, ".ref")
  pairs <- unique(data.frame(g = g[retained], t = tt[retained]))
  pairs <- pairs[order(pairs$g, pairs$t), , drop = FALSE]
  lev <- c(".ref", paste(pairs$g, pairs$t, sep = "::"))
  d$.cell_raw <- raw
  d$.cellf <- factor(vals, levels = lev)
  d
}

.jwdid_cells <- function(d, gvar, tvar, hettype, never, antigap,
                         ll = NULL, ul = NULL, recode = NULL, evbase = NULL) {
  if (hettype == "timecohort") {
    return(.jwdid_cells_timecohort(d, gvar, tvar, never, antigap))
  }
  pp <- attr(d, "prep"); gap <- pp$gap
  g <- d[[gvar]]; tt <- d[[tvar]]; antigap0 <- antigap - gap
  d$.post <- ifelse(g == 0, 0L,
    ifelse(tt < g - antigap, 1L,
      ifelse(tt >= g - antigap0, 2L, 0L)))
  d$.evnt <- (tt - g) * (g > 0) - gap * (g == 0)
  if (!is.null(ll) && !is.null(ul) && ll >= ul) stop("`hettype_ll` must be less than `hettype_ul`.", call.=FALSE)
  if (!is.null(ul)) d$.evnt <- pmin(d$.evnt, ul)
  # Stata jwdid.ado 2.201 applies ll for both control-group modes because
  # line 381 tests the literal "never". Preserve shipped behavior for parity.
  if (!is.null(ll)) d$.evnt <- pmax(d$.evnt, ll)

  retained <- g %in% pp$glist
  if (hettype == "cohort") {
    retained <- retained & if (never) d$.post %in% c(1L, 2L) else d$.post == 2L
    a <- g; b <- d$.post
  } else if (hettype == "time") {
    retained <- retained & if (never) d$.post %in% c(1L, 2L) else d$.post == 2L
    a <- tt; b <- d$.post
  } else if (hettype == "event") {
    retained <- retained & if (never) d$.evnt != -antigap else d$.evnt > -antigap
    a <- d$.evnt; b <- NULL
  } else if (hettype == "twfe") {
    retained <- retained & if (never) d$.post %in% c(1L, 2L) else d$.post == 2L
    a <- d$.post; b <- NULL
  } else if (hettype == "eventcohort") {
    if (!never) stop("Internal error: eventcohort must force never controls.", call.=FALSE)
    ev <- d$.evnt
    if (is.function(recode)) ev <- recode(ev)
    if (is.list(recode) && length(recode)) {
      ev0 <- ev
      assigned <- rep(FALSE, length(ev))
      for (nm in names(recode)) {
        hit <- !assigned & ev0 %in% recode[[nm]]
        ev[hit] <- as.numeric(nm)
        assigned[hit] <- TRUE
      }
    }
    if (any(ev < 0, na.rm=TRUE)) stop("Recoded eventcohort values must be non-negative.", call.=FALSE)
    if (is.null(evbase)) stop("`hettype_evbase` is required for eventcohort.", call.=FALSE)
    base <- evbase
    retained <- retained & ev != base
    a <- g; b <- ev
    d$.evnt <- ev
  } else stop("Unsupported `hettype`.", call.=FALSE)

  vals <- if (is.null(b)) as.character(a) else paste(a, b, sep="::")
  vals[!retained | is.na(retained)] <- ".ref"
  if (is.null(b)) {
    lev <- c(".ref", as.character(sort(unique(a[retained]))))
  } else {
    pairs <- unique(data.frame(a=a[retained], b=b[retained]))
    pairs <- pairs[order(pairs$a, pairs$b),,drop=FALSE]
    lev <- c(".ref", paste(pairs$a, pairs$b,sep="::"))
  }
  d$.cell_raw <- if (is.null(b)) factor(a) else interaction(a,b,drop=TRUE,lex.order=TRUE)
  d$.cellf <- factor(vals, levels=lev)
  attr(d, "prep") <- pp
  d
}
