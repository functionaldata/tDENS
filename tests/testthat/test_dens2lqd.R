library(testthat)

test_that('dens2lqd works on a trivial example, U(0,2)', {
  x = seq(0,2,length.out =512)
  y = rep(0.5,length.out =512)
  expect_equal( rep(-log(0.5), 512), dens2lqd(dens=y, dSup=x) ) 
})

test_that('dens2lqd works on a second trivial example, U(0,2)', {
  x =  seq(0,2,length.out = 512)
  y = (seq(0,2,length.out =512) + 0.1)/2.2
  yOnLQ <- dens2lqd( dens=y, dSup = x)
  expect_equal( yOnLQ[c(1,512)], -log( y[c(1,512)]))  
}) 
