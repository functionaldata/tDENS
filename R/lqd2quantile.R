#' Function for converting log quantile densities to quantile functions
#' 
#' @param  lqd log quantile density on lqdSup
#' @param  lqdSup support for lqd domain - must begin at 0 and end at 1
#' @param  lb lower bound of support for Density domain - default is 0. 
#' @param  useSplines fit spline to the lqd when doing the numerical integration (default: TRUE)
#' 
#' @return quantile values on lqdSup
#' 
#' @examples
#' x <- seq(1,3,length.out =512)
#' y.lqd <- rep(log(2), times = 512)
#' y <- lqd2quantile(lqd = y.lqd, lb = 1) # should equate # seq(1, 3, length.out = 512)
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2016} 
#' @export

lqd2quantile = function(lqd, lqdSup = seq(0, 1, length.out = length(lqd)), lb = 0, useSplines = TRUE){
  
  if(!all.equal( range(lqdSup),c(0,1) )){
    stop("Please check the support of the LQD domain's boundaries.")
  }
  
  if( useSplines ){
    # Could fit spline if this yields more accurate numerical integration
    lqd_sp = splinefun(lqdSup, lqd, method = 'natural')
    lqd_exp = function(t) exp(lqd_sp(t))
    # Get grid and function for density space
    Q = lb + c(0, cumsum(sapply(2:length(lqdSup), function(i) integrate(lqd_exp, lqdSup[i - 1], lqdSup[i])$value)))
  } else {
    # Get grid and function for density space
    Q = lb + cumtrapzRcpp(lqdSup, exp(lqd))
  }
  
  return(Q)
}
