test_that("unknown arguments error instead of leaking into feols", {
  d <- data.frame(id=rep(1:4,each=4),t=rep(1:4,4),g=rep(c(3,3,0,0),each=4),y=rnorm(16))
  expect_error(jwdid(y~1,d,ivar=id,tvar=t,gvar=g,nevr=TRUE),"Unknown")
  expect_true("xasis" %in% names(formals(jwdid)))
})

test_that("xasis changes raw coefficients but not aggregated ATT", {
  skip_if_not_installed("did")
  d <- did::mpdta
  a <- jwdid(lemp~lpop,d,ivar=countyreal,tvar=year,gvar="first.treat",never=TRUE,xasis=FALSE)
  b <- jwdid(lemp~lpop,d,ivar=countyreal,tvar=year,gvar="first.treat",never=TRUE,xasis=TRUE)
  expect_false(isTRUE(all.equal(unname(coef(a)),unname(coef(b)))))
  expect_equal(jwdid_simple(a)$estimate,jwdid_simple(b)$estimate,tolerance=1e-8)
})

test_that("covariate roles build distinct formula blocks", {
  skip_if_not_installed("did")
  d <- transform(did::mpdta,z=lpop^2)
  m <- jwdid(lemp~lpop,d,ivar=countyreal,tvar=year,gvar="first.treat",never=TRUE,
             xattvar=~z,exogvar=~z,xtvar=~z,xgvar=~z)
  expect_true(length(m$demean_info$used)==2)
  expect_true(is.finite(jwdid_simple(m)$estimate))
})

test_that("mpdta covariate ATT agrees with etwfe", {
  skip_if_not_installed("did"); skip_if_not_installed("etwfe")
  d <- did::mpdta
  m <- jwdid(lemp~lpop,d,ivar=countyreal,tvar=year,gvar="first.treat",never=TRUE)
  ref <- etwfe::etwfe(lemp~lpop,tvar=year,gvar=first.treat,data=d,
                      ivar=countyreal,cgroup="never",vcov=~countyreal)
  refa <- etwfe::emfx(ref,"simple")
  expect_equal(jwdid_simple(m)$estimate,refa$estimate,tolerance=1e-8)
})

test_that("role and fixed-effect variables participate in complete cases", {
  d <- expand.grid(id=1:8,t=1:5); d$g <- ifelse(d$id<=4,3,0)
  d$y <- rnorm(nrow(d)); d$z <- rnorm(nrow(d)); d$f <- rep(1:2,length.out=nrow(d))
  d$z[1] <- NA
  m <- jwdid(y~1,d,ivar=id,tvar=t,gvar=g,xattvar=~z,fevar=~f,never=TRUE)
  expect_false(1 %in% m$data$.jwdid_row)
})

test_that("group FE retains cohort slopes for unit-constant covariates", {
  d <- expand.grid(id=1:20,t=1:5); d$g <- ifelse(d$id<=10,3,0)
  d$zconst <- d$id; d$y <- sin(d$id+d$t)+(d$t>=d$g & d$g>0)
  mu <- jwdid(y~zconst,d,ivar=id,tvar=t,gvar=g,never=TRUE,group=FALSE)
  mg <- jwdid(y~zconst,d,ivar=id,tvar=t,gvar=g,never=TRUE,group=TRUE)
  expect_false(any(grepl("^g::.*jwdid_x",names(coef(mu)))))
  expect_true(any(grepl("g::.*jwdid_x",names(coef(mg)))))
})

test_that("weights and cluster accept bare string and formula forms", {
  d <- expand.grid(id=1:12,t=1:5); d$g <- ifelse(d$id<=6,3,0)
  d$w <- 1+d$id/100; d$y <- rnorm(nrow(d))+(d$t>=d$g & d$g>0)
  a <- jwdid(y~1,d,ivar=id,tvar=t,gvar=g,never=TRUE,weights=w,cluster=id)
  b <- jwdid(y~1,d,ivar="id",tvar="t",gvar="g",never=TRUE,
             weights="w",cluster="id")
  c <- jwdid(y~1,d,ivar=id,tvar=t,gvar=g,never=TRUE,weights=~w,cluster=~id)
  expect_equal(coef(a),coef(b),tolerance=1e-12)
  expect_equal(coef(a),coef(c),tolerance=1e-12)
})
