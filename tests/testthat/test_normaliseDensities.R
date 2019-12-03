library(testthat)

test_that('normaliseDensities works on a trivial example of truncated normal', {
  x = seq(-2, 2, length.out = 101)
  y = matrix(dnorm(x), nrow = 1)
  z = normaliseDensities(dmatrix = y, dSup = x)

  expect_equal( 1 , trapzRcpp(X = x, Y = z) )
})
