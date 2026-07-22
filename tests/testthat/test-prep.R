test_that("gap anticipation and switch boundary reproduce the ado", {
  d <- expand.grid(id = 1:3, t = seq(2000, 2015, 5))
  d$g <- c(2005, 2010, 0)[d$id]
  d$y <- rnorm(nrow(d))
  p <- jwdid:::.jwdid_prep(y ~ 1, d, "id", "t", "g", anticipation = 0)
  expect_equal(attr(p, "prep")$gap, 5)
  expect_equal(p$.tr[p$id == 1], c(1, 1, 1, 1))
  expect_equal(p$.etr[p$id == 1], c(0, 1, 1, 1))
})

test_that("always-treated units are removed and invalid intensity errors", {
  d <- data.frame(id=rep(1:2, each=3), t=rep(1:3,2), g=rep(c(1,0),each=3), y=1:6, z=0)
  p <- jwdid:::.jwdid_prep(y ~ 1, d, "id", "t", "g")
  expect_false(any(p$id == 1))
  d$z[1] <- 2
  expect_error(jwdid:::.jwdid_prep(y ~ 1, d, "id", "t", "g", "z"), "\\[0, 1\\]")
})
