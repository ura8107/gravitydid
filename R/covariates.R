.jwdid_mm <- function(fml, data, prefix) {
  if (is.null(fml)) return(list(data=data, names=character()))
  mm <- stats::model.matrix(fml, data=data)
  mm <- mm[, colnames(mm) != "(Intercept)", drop=FALSE]
  if (!ncol(mm)) return(list(data=data, names=character()))
  nm <- paste0(prefix, seq_len(ncol(mm)))
  data[nm] <- unclass(as.data.frame(mm))
  list(data=data, names=nm, labels=colnames(mm))
}

.jwdid_rhs_formula <- function(fml) {
  rhs <- fml[[3L]]
  if (deparse1(rhs) %in% c("1", "0")) NULL else stats::as.formula(call("~", rhs), env=environment(fml))
}

.jwdid_is_varying <- function(data, vars, ivar) {
  if (!length(vars)) return(logical())
  if (is.null(ivar)) return(setNames(rep(TRUE,length(vars)),vars))
  setNames(vapply(vars,function(v) {
    any(vapply(split(data[[v]],data[[ivar]]),function(z)
      diff(range(z,na.rm=TRUE)) != 0,logical(1)))
  },logical(1)),vars)
}

.jwdid_demean <- function(data, vars, cell, weights=NULL) {
  if (!length(vars)) return(list(data=data,names=character()))
  X <- as.matrix(data[,vars,drop=FALSE])
  w <- if (is.null(weights)) NULL else data[[weights]]
  Z <- fixest::demean(X, f=list(cell), weights=w)
  nm <- paste0(".jwdid_dm",seq_len(ncol(Z)))
  data[nm] <- unclass(as.data.frame(Z))
  list(data=data,names=nm)
}

.jwdid_add_formula_vars <- function(fml, data, prefix) {
  if (is.null(fml)) return(list(data=data,names=character()))
  .jwdid_mm(fml,data,prefix)
}

.jwdid_unit_means <- function(data, rhs, ivar, weights=NULL, prefix=".jwdid_corr",
                              drop_identical=FALSE) {
  mm <- stats::model.matrix(stats::as.formula(paste("~",rhs)),data=data)
  mm <- mm[,colnames(mm)!="(Intercept)",drop=FALSE]
  if (!ncol(mm)) return(list(data=data,names=character()))
  ids <- data[[ivar]]; w <- if(is.null(weights)) rep(1,nrow(data)) else data[[weights]]
  means <- vapply(seq_len(ncol(mm)),function(j)
    ave(seq_len(nrow(mm)),ids,FUN=function(ii)
      rep(sum(mm[ii,j]*w[ii])/sum(w[ii]),length(ii))),numeric(nrow(mm)))
  if (is.null(dim(means))) means <- matrix(means,ncol=1)
  keep <- vapply(seq_len(ncol(means)),function(j) {
    stats::var(means[,j]) > .Machine$double.eps &&
      (!drop_identical || !isTRUE(all.equal(means[,j],mm[,j],tolerance=1e-12)))
  },logical(1))
  means <- means[,keep,drop=FALSE]
  if (!ncol(means)) return(list(data=data,names=character()))
  # Deterministic QR pruning mirrors _rmcoll before fixest's own collinearity pass.
  q <- qr(means,tol=1e-10); take <- sort(q$pivot[seq_len(q$rank)])
  means <- means[,take,drop=FALSE]
  nm <- paste0(prefix,seq_len(ncol(means)))
  data[nm] <- unclass(as.data.frame(means))
  list(data=data,names=nm)
}
