test_that("timecohort levels use numeric pair ordering", {
  d <- expand.grid(id=1:3, t=c(1,3,10,11,12))
  d$g <- c(3,10,0)[d$id]; d$y <- seq_len(nrow(d))
  p <- jwdid:::.jwdid_prep(y~1,d,"id","t","g",never=TRUE)
  z <- jwdid:::.jwdid_cells(p,"g","t","timecohort",TRUE,attr(p,"prep")$antigap)
  lev <- levels(z$.cellf)[-1]
  expect_equal(lev, lev[order(as.numeric(sub("::.*","",lev)),
                              as.numeric(sub(".*::","",lev)))])
})

test_that("all phase 2 hettypes fit and event restriction agrees", {
  skip_if_not_installed("did")
  args <- list(fml=lemp~1, data=did::mpdta, ivar="countyreal", tvar="year",
               gvar="first.treat", never=TRUE)
  fits <- lapply(c("cohort","time","event","twfe"), function(h)
    do.call(jwdid, c(args, list(hettype=h))))
  expect_true(all(vapply(fits, inherits, logical(1), "jwdid")))
  direct <- jwdid_event(fits[[3]])
  expect_true(all(is.finite(direct$estimate)))
  # An event-restricted model aggregates its own event coefficients as a no-op.
  b <- coef(fits[[3]])
  event_names <- sub(":\\.tr$", "", sub("^.*\\.cellf::", "", names(b)))
  lookup <- setNames(as.numeric(b), event_names)
  expected <- setNames(vapply(direct$term, function(e)
    if (e == "-1") 0 else lookup[[e]], numeric(1)), direct$term)
  expect_equal(setNames(direct$estimate, direct$term), expected, tolerance=1e-10)

  # With multiple cohorts, the restricted and unrestricted models are distinct.
  unrestricted <- do.call(jwdid, args)
  expect_gt(max(abs(jwdid_event(unrestricted)$estimate - direct$estimate)), 1e-4)
})

test_that("event and timecohort are identical with one treated cohort", {
  set.seed(42)
  d <- expand.grid(id=1:80, t=1:6)
  d$g <- ifelse(d$id <= 40, 4, 0)
  d$y <- rnorm(nrow(d)) + 0.2*d$t + ifelse(d$g>0 & d$t>=d$g, d$t-d$g+1, 0)
  mt <- jwdid(y~1,d,ivar=id,tvar=t,gvar=g,never=TRUE,hettype="timecohort")
  me <- jwdid(y~1,d,ivar=id,tvar=t,gvar=g,never=TRUE,hettype="event")
  expect_equal(jwdid_event(mt)$estimate, jwdid_event(me)$estimate, tolerance=1e-10)
})

test_that("event binning validates bounds", {
  skip_if_not_installed("did")
  expect_error(jwdid(lemp~1,did::mpdta,ivar=countyreal,tvar=year,
    gvar="first.treat",never=TRUE,hettype="event",hettype_ll=2,hettype_ul=1),
    "less than")
})

test_that("event ll and ul censor retained event cells", {
  skip_if_not_installed("did")
  m <- jwdid(lemp~1,did::mpdta,ivar=countyreal,tvar=year,
    gvar="first.treat",never=TRUE,hettype="event",hettype_ll=-2,hettype_ul=2)
  lev <- as.numeric(levels(m$data$.cellf)[-1])
  expect_true(all(lev >= -2 & lev <= 2))
  expect_false(-1 %in% lev)

  mn <- jwdid(lemp~1,did::mpdta,ivar=countyreal,tvar=year,
    gvar="first.treat",never=FALSE,hettype="event",hettype_ll=0,hettype_ul=2)
  # Preserve the shipped Stata behavior: a macro typo applies ll to notyet too.
  expect_true(all(mn$data$.evnt >= 0))
  expect_true(all(as.numeric(levels(mn$data$.cellf)[-1]) > -mn$antigap))
})

test_that("eventcohort forces never and supports recode plus evbase", {
  skip_if_not_installed("did")
  recode <- function(e) e + 5
  expect_warning(m <- jwdid(lemp~1,did::mpdta,ivar=countyreal,tvar=year,
    gvar="first.treat",never=FALSE,hettype="eventcohort",
    hettype_recode=recode,hettype_evbase=4), "forces")
  expect_equal(m$type,"never")
  expect_false(any(grepl("::4$",levels(m$data$.cellf))))
  expect_true(all(m$data$.evnt >= 0))
  expect_true(is.finite(jwdid_simple(m)$estimate))

  expect_error(suppressWarnings(jwdid(lemp~1,did::mpdta,ivar=countyreal,tvar=year,
    gvar="first.treat",hettype="eventcohort")), "non-negative")
})

test_that("eventcohort affine recode is exactly timecohort", {
  skip_if_not_installed("did")
  mt <- jwdid(lemp~1,did::mpdta,ivar=countyreal,tvar=year,
    gvar="first.treat",never=TRUE,hettype="timecohort")
  me <- jwdid(lemp~1,did::mpdta,ivar=countyreal,tvar=year,
    gvar="first.treat",never=TRUE,hettype="eventcohort",
    hettype_recode=function(e)e+5,hettype_evbase=4)
  expect_equal(unname(coef(me)),unname(coef(mt)),tolerance=1e-12)
  expect_equal(jwdid_simple(me)$estimate,jwdid_simple(mt)$estimate,tolerance=1e-12)
})

test_that("eventcohort list recode uses original values and first match", {
  skip_if_not_installed("did")
  rc <- list("5"=c(-4,-3), "6"=c(5), "0"=c(-2), "1"=c(-1),
             "2"=c(0), "3"=c(1), "4"=c(2), "7"=c(3))
  m <- jwdid(lemp~1,did::mpdta,ivar=countyreal,tvar=year,
    gvar="first.treat",never=TRUE,hettype="eventcohort",
    hettype_recode=rc,hettype_evbase=1)
  expect_true(5 %in% m$data$.evnt)
  expect_false(6 %in% m$data$.evnt)

  rc_first <- list("0"=c(-4), "9"=c(-4), "1"=c(-3,-2,-1,0,1,2,3))
  m2 <- jwdid(lemp~1,did::mpdta,ivar=countyreal,tvar=year,
    gvar="first.treat",never=TRUE,hettype="eventcohort",
    hettype_recode=rc_first,hettype_evbase=1)
  expect_true(0 %in% m2$data$.evnt)
  expect_false(9 %in% m2$data$.evnt)
})

test_that("eventcohort requires evbase, permits zero, and ignores ll ul with warning", {
  skip_if_not_installed("did")
  expect_error(suppressWarnings(jwdid(lemp~1,did::mpdta,ivar=countyreal,tvar=year,
    gvar="first.treat",hettype="eventcohort",hettype_recode=function(e)e+5)),
    "required")
  expect_warning(jwdid(lemp~1,did::mpdta,ivar=countyreal,tvar=year,
    gvar="first.treat",never=TRUE,hettype="eventcohort",
    hettype_recode=function(e)e+5,hettype_evbase=4,hettype_ll=0), "ignored")
})

test_that("post reference period is Base rather than missing", {
  d <- expand.grid(id=1:3,t=1:5); d$g <- c(3,4,0)[d$id]; d$y <- seq_len(nrow(d))
  p <- jwdid:::.jwdid_prep(y~1,d,"id","t","g",never=TRUE)
  z <- jwdid:::.jwdid_cells(p,"g","t","cohort",TRUE,attr(p,"prep")$antigap)
  expect_true(all(z$.post[z$g>0 & z$t==z$g-1] == 0))
  expect_false(anyNA(z$.cell_raw))
})

test_that("retained conditions follow the hettype table", {
  d <- expand.grid(id=1:3,t=1:5)
  d$g <- c(3,4,0)[d$id]; d$y <- seq_len(nrow(d))
  for (never in c(FALSE,TRUE)) {
    p <- jwdid:::.jwdid_prep(y~1,d,"id","t","g",never=never)
    for (h in c("timecohort","cohort","time","event","twfe")) {
      z <- jwdid:::.jwdid_cells(p,"g","t",h,never,attr(p,"prep")$antigap)
      expect_true(".ref" %in% levels(z$.cellf))
      expect_false(anyNA(z$.cellf))
      if (!never) expect_true(all(z$.etr[z$.cellf != ".ref"] == 1))
    }
  }
})
