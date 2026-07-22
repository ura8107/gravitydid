test_that("attgt returns full grid and reproduces retained cell coefficients", {
  skip_if_not_installed("did")
  m <- jwdid(lemp~1,did::mpdta,ivar=countyreal,tvar=year,gvar="first.treat",never=TRUE)
  a <- jwdid_attgt(m)
  expect_true(any(as.numeric(sub(".*::","",a$term)) < as.numeric(sub("::.*","",a$term))))
  b <- coef(m); nm <- sub(":\\.tr$","",sub("^.*\\.cellf::","",names(b)))
  hit <- intersect(nm,a$term)
  expect_equal(setNames(a$estimate,a$term)[hit],setNames(as.numeric(b),nm)[hit],tolerance=1e-10)
})

test_that("event window selects while cwindow pools", {
  skip_if_not_installed("did")
  m <- jwdid(lemp~1,did::mpdta,ivar=countyreal,tvar=year,gvar="first.treat",never=TRUE)
  full <- jwdid_event(m); win <- jwdid_event(m,window=c(-2,2)); cens <- jwdid_event(m,cwindow=c(-2,2))
  ii <- as.numeric(full$term)>=-2 & as.numeric(full$term)<=2
  expect_equal(win$estimate,full$estimate[ii],tolerance=1e-10)
  expect_equal(as.numeric(cens$term),-2:2)
  expect_lt(length(cens$term),length(full$term))
})

test_that("pretrend Wald equals a direct eigen generalized inverse", {
  skip_if_not_installed("did")
  m <- jwdid(lemp~1,did::mpdta,ivar=countyreal,tvar=year,gvar="first.treat",never=TRUE)
  a <- jwdid_event(m,pretrend=TRUE); ii <- as.numeric(a$term) < -1
  V <- a$V[ii,ii,drop=FALSE]; b <- a$b[ii]; e <- eigen(V,symmetric=TRUE)
  keep <- e$values > max(dim(V))*max(abs(e$values))*.Machine$double.eps
  Vi <- e$vectors[,keep,drop=FALSE]%*%diag(1/e$values[keep],sum(keep))%*%t(e$vectors[,keep,drop=FALSE])
  expect_equal(a$pretrend$chi2,as.numeric(t(b)%*%Vi%*%b),tolerance=1e-10)
})

test_that("restriction, over validation, any and aggregation weights work", {
  skip_if_not_installed("did")
  d <- transform(did::mpdta,region=countyreal%%2,w2=1+(countyreal%%3))
  m <- jwdid(lemp~1,d,ivar=countyreal,tvar=year,gvar="first.treat",never=TRUE)
  r <- jwdid_simple(m,orestriction=~region==0)
  manual <- jwdid:::.jwdid_analytic(m,m$data,which(m$data$.etr==1 & m$data$region==0),NULL)
  expect_equal(r$estimate,manual$table$estimate,tolerance=1e-10)
  expect_gt(length(jwdid_any(m,over=~region)$term),1)
  expect_true(is.finite(jwdid_simple(m,weights=~w2)$estimate))
  expect_error(jwdid_group(m,over=~region),"only for")
  expect_error(jwdid_simple(m,over2=~region),"only for")
})

test_that("continuous treatment derives gvar and asis uses observed intensity", {
  set.seed(20)
  d <- expand.grid(id=1:30,t=1:5)
  d$z <- ifelse(d$id<=10 & d$t>=3,.5,ifelse(d$id<=20 & d$t>=4,.8,0))
  d$y <- .3*d$z + .1*d$t + rnorm(nrow(d),sd=.1)
  m <- jwdid(y~1,d,ivar=id,tvar=t,trtvar=z,never=TRUE)
  expect_equal(sort(unique(m$data[[m$gvar]])),c(0,3,4))
  unit <- jwdid_simple(m,asis=FALSE)
  obs <- jwdid_simple(m,asis=TRUE)
  expect_true(all(is.finite(c(obs$estimate,obs$std.error))))
  expect_false(isTRUE(all.equal(unit$estimate,obs$estimate)))
  d$z[1] <- 1.2
  expect_error(jwdid(y~1,d,ivar=id,tvar=t,trtvar=z),"\\[0, 1\\]")
})
