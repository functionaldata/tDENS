library(testthat)

test_that('normaliseDensities works on a trivial example, U(0,2)', {
  x = seq(0,2,length.out =512)
  y = matrix(rep(1,length.out =512), nrow = 1)
  
  z = normaliseDensities(dmatrix = y, dSup = x)
  expect_equal(rep(0.5, 512), as.vector(z)) 
  
})
