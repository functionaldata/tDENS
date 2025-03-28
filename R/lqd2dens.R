#' Function for converting log quantile densities to densities
#' 
#' @param  lqd log quantile density on lqdSup
#' @param  lqdSup support forlqd domain - must begin at 0 and end at 1
#' @param  dSup support for Density domain - max and min values mark the boundary of the support.
#' @param  useSplines fit spline to the lqd when doing the numerical integration (default: TRUE)
#' 
#' @return dens density values on dSup
#' 
#' @examples
#' x <- seq(0,2,length.out =512)
#' y.lqd <- rep(log(2), times = 512)
#' y <- lqd2dens(dSup=x, lqd = y.lqd) # should equate # 1/2
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2016} 
#' @export

lqd2dens = function(lqd, lqdSup = seq(0, 1, length.out = length(lqd)), dSup, useSplines = TRUE){

  if(!all.equal( range(lqdSup),c(0,1) )){
    stop("Please check the support of the LQD domain's boundaries.")
  }
  
  if( useSplines ){
    # Could fit spline if this yields more accurate numerical integration
    lqd_sp = splinefun(lqdSup, lqd, method = 'natural')
    lqd_exp = function(t) exp(lqd_sp(t))
    # Get grid and function for density space
    dtemp = dSup[1] + c(0, cumsum(sapply(2:length(lqdSup), function(i) integrate(lqd_exp, lqdSup[i - 1], lqdSup[i])$value)))
  } else {
    # Get grid and function for density space
    dtemp = dSup[1] + cumtrapzRcpp(lqdSup, exp(lqd))
  }

  dens_temp = exp(-lqd);
  
  # Fix the support of the density to match dSup
  
  r1 = diff(range(dtemp)); r2 = diff(range(dSup))
  dtemp = (dtemp - dtemp[1])*r2/r1 + dtemp[1]; 

  # Remove duplicates
  ind = duplicated(dtemp)
  dtemp = unique(dtemp)
  
  # Interpolate to dSup and normalize
  dens = approx(x = dtemp, y = dens_temp[!ind], xout = dSup, rule = c(2,2))[[2]]
  dens = dens/trapzRcpp(X = dSup,Y = dens); # Normalize

  return(dens)
}
