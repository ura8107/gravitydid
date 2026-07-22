test_that("symbolic .tr round-trips through fixest prediction", {
  skip_if_not_installed("did")
  m <- jwdid(lemp ~ 1, did::mpdta, ivar=countyreal, tvar=year,
             gvar="first.treat", never=TRUE)
  d1 <- m$data; d0 <- d1; d0$.tr <- 0
  expect_false(isTRUE(all.equal(predict(m$model, d1), predict(m$model, d0))))
  expect_s3_class(jwdid_aggte(m, "simple"), "jwdid_aggte")
  expect_true(all(is.finite(jwdid_aggte(m, "event")$estimate)))
  a <- jwdid_aggte(m, "simple", engine="analytic")
  b <- jwdid_aggte(m, "simple", engine="marginaleffects")
  expect_equal(a$estimate, b$estimate, tolerance=1e-8)
  expect_equal(a$std.error, b$std.error, tolerance=1e-6)
})

test_that("mpdta simple ATT agrees with etwfe", {
  skip_if_not_installed("did")
  skip_if_not_installed("etwfe")
  d <- did::mpdta
  m <- jwdid(lemp ~ 1, d, ivar=countyreal, tvar=year,
             gvar="first.treat", never=TRUE)
  ref <- etwfe::etwfe(lemp ~ 1, tvar=year, gvar=first.treat, data=d,
                      ivar=countyreal, cgroup="never", vcov=~countyreal)
  refa <- etwfe::emfx(ref, "simple")
  expect_equal(jwdid_aggte(m, "simple")$estimate, refa$estimate, tolerance=1e-8)
  expect_equal(jwdid_aggte(m, "simple")$std.error, refa$std.error, tolerance=1e-6)
})

test_that("mpdta treatment cells and full event path match universal-base did", {
  skip_if_not_installed("did")
  m <- jwdid(lemp ~ 1, did::mpdta, ivar=countyreal, tvar=year,
             gvar="first.treat", never=TRUE)
  expected_cells <- c(
    "2004::2004"=-0.010503246221, "2004::2005"=-0.070423158103,
    "2004::2006"=-0.137258738889, "2004::2007"=-0.100811363085,
    "2006::2003"=-0.003769293674, "2006::2004"= 0.002750818751,
    "2006::2006"=-0.004594606953, "2006::2007"=-0.041224471546,
    "2007::2003"= 0.003306356693, "2007::2004"= 0.033813012276,
    "2007::2005"= 0.031087119390, "2007::2007"=-0.026054410719)
  got <- coef(m)
  names(got) <- sub("^.*\\.cellf::", "", names(got))
  names(got) <- sub(":\\.tr$", "", names(got))
  expect_equal(got[names(expected_cells)], expected_cells, tolerance=1e-8)

  ev <- jwdid_aggte(m, "event")
  expected_event <- c("-4"=0.003306356693, "-3"=0.025021829598,
    "-2"=0.024458744971, "-1"=0, "0"=-0.019931816789,
    "1"=-0.050957367065, "2"=-0.137258738889, "3"=-0.100811363085)
  expect_equal(setNames(ev$estimate, ev$term), expected_event, tolerance=1e-8)
  ref <- which(ev$term == "-1")
  expect_equal(ev$estimate[ref], 0)
  expect_true(is.na(ev$statistic[ref]))
  expect_true(is.na(ev$p.value[ref]))
})
