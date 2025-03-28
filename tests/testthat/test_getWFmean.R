library(testthat)

test_that('getWFmean works for simple example of linear densities', {
  x = seq(0,1,length.out =512)
  y = t(sapply(seq(0.5, 1.5, length.out = 4), function(b) b + 2*(1 - b)*x))
  y.qd = t(sapply(seq(0.5, 1.5, length.out = 4), function(b) (b^2 + 4*(1-b)*x)^(-1/2)))
  expect_equal( getWFmean(dmatrix = y, dSup = x), qd2dens(qd = colMeans(y.qd), qdSup = x, dSup = x) , tol = 1e-2) 
})