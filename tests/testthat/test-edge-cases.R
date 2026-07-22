test_that("no-never fallback, gap anticipation, and repeated cross section run", {
  d <- expand.grid(id=1:8,t=1:6)
  d$g <- rep(c(3,3,4,4,5,5,6,6),each=6)
  d$y <- sin(d$id + d$t) + (d$t>=d$g) * (d$t-d$g+1)
  p <- jwdid:::.jwdid_prep(y~1,d,"id","t","g")
  expect_equal(attr(p,"prep")$glist,c(3,4,5))
  expect_equal(attr(p,"prep")$gvarmax,6)

  q <- expand.grid(id=1:4,t=seq(2000,2020,5)); q$g <- c(2010,2015,0,0)[q$id]
  q$y <- seq_len(nrow(q))
  mq <- jwdid(y~1,q,ivar=id,tvar=t,gvar=g,anticipation=1,never=TRUE)
  expect_equal(mq$gap,5); expect_equal(mq$antigap,10)

  rc <- transform(q, id=NULL)
  expect_s3_class(jwdid(y~1,rc,tvar=t,gvar=g,never=TRUE),"jwdid")
})
