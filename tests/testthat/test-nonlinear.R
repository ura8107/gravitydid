test_that("nonlinear methods dispatch and force group FE where required", {
  set.seed(10)
  d <- expand.grid(id=1:30,t=1:5); d$g <- ifelse(d$id<=15,3,0)
  eta <- .2*d$t + .3*(d$g>0 & d$t>=d$g)
  d$y <- rpois(nrow(d),exp(eta))
  mp <- jwdid(y~1,d,ivar=id,tvar=t,gvar=g,never=TRUE,method="fepois")
  mg <- jwdid(y~1,d,ivar=id,tvar=t,gvar=g,never=TRUE,method="poisson")
  expect_false(mp$group); expect_true(mg$group)
  expect_s3_class(mp$model,"fixest"); expect_s3_class(mg$model,"fixest")
  ap <- jwdid_simple(mp,engine="analytic")
  ag <- jwdid_simple(mg)
  expect_true(all(is.finite(c(ap$estimate,ap$std.error,ag$estimate,ag$std.error))))
  expect_gt(length(coef(mg)),length(coef(mp)))
  expect_error(jwdid_simple(mp,engine="marginaleffects"),"cannot compute SEs")
})

test_that("fepois backend equals a hand refit of the assembled model", {
  set.seed(11)
  d <- expand.grid(id=1:24,t=1:5); d$g <- ifelse(d$id<=12,3,0)
  d$y <- rpois(nrow(d),exp(.1*d$t+.2*(d$g>0 & d$t>=d$g)))
  m <- jwdid(y~1,d,ivar=id,tvar=t,gvar=g,never=TRUE,method="fepois")
  ref <- fixest::fepois(stats::formula(m$model),data=m$data,vcov=~id,
                        fixef.rm="none",notes=FALSE)
  expect_equal(coef(m),coef(ref),tolerance=1e-10)
})

test_that("balanced panel flag is recorded", {
  d <- expand.grid(id=1:10,t=1:4); d$g <- ifelse(d$id<=5,3,0); d$y <- rnorm(nrow(d))
  expect_true(jwdid(y~1,d,ivar=id,tvar=t,gvar=g,never=TRUE)$balanced)
  du <- d[-1,]
  expect_false(jwdid(y~1,du,ivar=id,tvar=t,gvar=g,never=TRUE)$balanced)
})

test_that("corr adds nonconstant unit fitted-value controls", {
  set.seed(12)
  d <- expand.grid(id=1:30,t=1:6); d$g <- ifelse(d$id<=10,3,ifelse(d$id<=20,5,0))
  d <- d[!(d$id %% 4 == 0 & d$t == 2),]
  d$x <- d$id/10 + rnorm(nrow(d)); d$y <- .2*d$x + rnorm(nrow(d)) + (d$g>0 & d$t>=d$g)
  a <- jwdid(y~x,d,ivar=id,tvar=t,gvar=g,never=TRUE,group=TRUE,corr=FALSE)
  b <- jwdid(y~x,d,ivar=id,tvar=t,gvar=g,never=TRUE,group=TRUE,corr=TRUE)
  expect_length(a$correction_names,0)
  expect_gt(length(b$correction_names),0)
  expect_true(all(b$correction_names %in% names(b$data)))
  expect_false(isTRUE(all.equal(coef(a),coef(b))))
})

test_that("cre adds Mundlak means and omits cohort fixed effects", {
  set.seed(13)
  d <- expand.grid(id=1:40,t=1:5); d$g <- ifelse(d$id<=20,3,0)
  d$x <- d$id/20+rnorm(nrow(d)); eta <- -.5+.1*d$t+.15*d$x+.2*(d$g>0 & d$t>=d$g)
  d$y <- rpois(nrow(d),exp(eta))
  m <- jwdid(y~x,d,ivar=id,tvar=t,gvar=g,never=TRUE,method="poisson",cre=TRUE)
  expect_gt(length(m$cre_names),0)
  expect_true(any(grepl("jwdid_cre",names(coef(m)))))
  expect_false(any(grepl("^fixest::g::[^:]+$",names(coef(m)))))
  expect_true(all(is.finite(c(jwdid_simple(m)$estimate,jwdid_simple(m)$std.error))))
  expect_error(jwdid(y~1,d,ivar=id,tvar=t,gvar=g,cre=TRUE),"available")
})
