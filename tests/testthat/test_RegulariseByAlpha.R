library(testthat)

test_that('RegulariseByAlpha works on a trivial example, Beta(2,1))', {
  x = seq(0,1,length.out =122)
  y = seq(0,2,length.out =122)
  z = RegulariseByAlpha(x=x, y=y, alpha = 0.1)

  expect_equal( 1 , fdapace:::trapzRcpp(x,Y = z) )
  expect_equal( 0.1, min(z) ) 
})