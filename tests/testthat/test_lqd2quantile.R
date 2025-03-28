library(testthat)

test_that('lqd2quantile works on two trivial examples, U(0,2) and U(1,3)', {
  x = seq(0,1,length.out =512)
  y.lqd = rep(log(2),length.out =512)
  expect_equal( seq(0, 2, length.out = 512), lqd2quantile(lqd=y.lqd, lqdSup=x))
  expect_equal( seq(1, 3, length.out = 512), lqd2quantile(lqd=y.lqd, lqdSup=x, lb = 1))
})