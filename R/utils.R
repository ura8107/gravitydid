.jwdid_name <- function(expr, data, arg, env = parent.frame()) {
  if (is.null(expr)) return(NULL)
  val <- tryCatch(eval(expr, env), error = function(e) NULL)
  if (is.character(val) && length(val) == 1L && val %in% names(data)) return(val)
  if (inherits(val, "formula")) {
    vars <- all.vars(val)
    if (length(vars) == 1L && vars %in% names(data)) return(vars)
  }
  nm <- deparse(expr)
  if (nm %in% names(data)) return(nm)
  stop("`", arg, "` must name one column in `data`.", call. = FALSE)
}

.jwdid_gap <- function(g, t) {
  z <- sort(unique(c(g[g > 0], t[t > 0])))
  dz <- diff(z)
  if (!length(dz)) 1 else min(dz[dz > 0])
}

.jwdid_complete <- function(fml, data, cols, weights = NULL) {
  vars <- unique(c(all.vars(fml), cols, weights))
  stats::complete.cases(data[, vars, drop = FALSE])
}
