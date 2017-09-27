library(testthat)

test_that('DeegulariseByAlpha works on a trivial example, Beta(2,1))', {
  x = seq(0,1,length.out =122)
  y = seq(0,2,length.out =122)
  z = RegulariseByAlpha(x=x, y=y, alpha = 0.1)
  yNew = DeregulariseByAlpha(x=x, y=z, alpha = 0)
  
  expect_equal( 1 , trapzRcpp(x,Y = yNew) )
  expect_equal( y, yNew , tol = 1e-10)
})
