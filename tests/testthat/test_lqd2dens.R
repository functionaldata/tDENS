library(testthat)

test_that('lqd2dens works on a trivial example, U(0,2)', {
  x = seq(0,1,length.out =512)
  y.lqd = rep(log(2),length.out =512)
  expect_equal( rep(0.5, 512), lqd2dens(lqd=y, lqdSup=x, dSup = seq(0, 1, length.out = 512) )) 
})

test_that('lqd2dens works on a second trivial example', {
  x =  seq(0,2,length.out = 512)
  y = (seq(0,2,length.out =512) + 0.1)/2.2
  yOnLQ <- dens2lqd( dens=y, dSup = x)
  expect_equal( lqd2dens(lqd = yOnLQ, lqdSup = x, dSup = seq(0, 1, length.out = 512))[c(1,512)], exp( -y[c(1,512)]))  
}) 