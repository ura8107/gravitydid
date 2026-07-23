test_that("tidy glance summary and plot methods work", {
  skip_if_not_installed("did")
  m <- jwdid(lemp~1,did::mpdta,ivar=countyreal,tvar=year,gvar="first.treat",never=TRUE)
  a <- jwdid_event(m)
  expect_s3_class(generics::tidy(m),"data.frame")
  expect_s3_class(generics::tidy(a),"data.frame")
  expect_s3_class(generics::glance(m),"data.frame")
  expect_s3_class(generics::glance(a),"data.frame")
  expect_s3_class(summary(a),"data.frame")
  skip_if_not_installed("ggplot2")
  for(s in c("errorbar","ribbon","pointrange","bar")) expect_s3_class(plot(a,style=s),"ggplot")
  expect_s3_class(plot(a,pre=list(colour="navy"),post=list(colour="firebrick")),"ggplot")
  if (requireNamespace("modelsummary",quietly=TRUE) &&
      requireNamespace("broom",quietly=TRUE))
    expect_s3_class(modelsummary::modelsummary(m,output="data.frame"),"data.frame")
})

test_that("aggte rejects misspelled options", {
  d <- data.frame(id=rep(1:4,each=4),t=rep(1:4,4),g=rep(c(3,3,0,0),each=4),y=rnorm(16))
  m <- jwdid(y~1,d,ivar=id,tvar=t,gvar=g,never=TRUE)
  expect_error(jwdid_aggte(m,"event",windo=c(-2,2)),"Unknown")
  expect_error(jwdid_aggte(m,"simple",orestrict=~id==1),"Unknown")
})
