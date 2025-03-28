library(testthat)

test_that('dens2qd works on a trivial example, U(0,2)', {
    x = seq(0,2,length.out =512)
    y = rep(0.5,length.out =512)
    expect_equal( rep(2, 512), dens2qd(dens=y, dSup=x) ) 
})