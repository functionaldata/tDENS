library(testthat)

test_that('dens2quantile works on a trivial example, U(0,2)', {
  x = seq(0,2,length.out =512)
  y = rep(0.5,length.out =512)
  expect_equal( 2*seq(0, 1, length.out = 512), dens2quantile(dens=y, dSup=x) ) 
})