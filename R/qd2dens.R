#' Function for converting Quantile Densities to Densities
#' 
#' @param  qd quantile density on qdSup
#' @param  qdSup support for quantile domain - must begin at 0 and end at 1 (default = seq(0, 1, length.out = length(qd)))
#' @param  dSup support for Density domain - max and min values mark the boundary of the support.
#' @param  useSplines fit spline to the qd when doing the numerical integration (default: TRUE)
#' 
#' @return dens density values on dSup
#' 
#' @examples 
#' 
#' x <- seq(0,1,length.out =512)
#' y <- rep(2,length.out =512)
#' y.dens <- qd2dens(qd=y, qdSup = x, dSup = seq(0, 2, length.out = 512)) # should equate # 1/2
#' 
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2016} 
#' @export

qd2dens = function(qd, qdSup = seq(0, 1, length.out = length(qd)), dSup, useSplines = TRUE){

  if(!all.equal( range(qdSup),c(0,1) )){
    stop("Please check the support of the QD domain's boundaries.")
  }
  
  if(abs(fdapace:::trapzRcpp(qdSup, qd) - diff(range(dSup)) > 1e-5)){
    stop("Quantile Density should integrate to the range of dSup.")
  }
  
  if( useSplines ){
    # Could fit spline if this yields more accurate numerical integration
    qd_sp = splinefun(qdSup, qd, method = 'natural')
    # Get grid and function for density space
    dtemp = dSup[1] + c(0, cumsum(sapply(2:length(qdSup), function(i) integrate(qd_sp, qdSup[i - 1], qdSup[i])$value)))
  } else {
    # Get grid and function for density space
    dtemp = dSup[1] + fdapace:::cumtrapzRcpp(qdSup, qd)
  }

  dens_temp = 1/qd;
  
  # Remove duplicates
  ind = duplicated(dtemp)
  dtemp = unique(dtemp)
  
  # Interpolate to dSup and normalize
  dens = approx(x = dtemp, y = dens_temp[!ind], xout = dSup, rule = c(2,2))[[2]]
  dens = dens/fdapace:::trapzRcpp(X = dSup,Y = dens); # Normalize

  return(dens)
}