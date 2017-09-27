library(testthat)

test_that('qd2dens works on a trivial example, U(0,2)', {
  x = seq(0,1,length.out =512)
  y.qd = rep(2,length.out =512)
  expect_equal( rep(0.5, 512), qd2dens(qd=y.qd, qdSup=x, dSup = seq(0, 2, length.out = 512) )) 
})
