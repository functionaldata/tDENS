#' Function for converting Densities to Quantile Densities
#' 
#' @param  dens density on dSup
#' @param  dSup support for Density domain - max and min values mark the boundary of the support.
#' @param  qdSup support for quantile density domain - must begin at 0 and end at 1
#' @param  useSplines fit spline to the qd when doing the numerical integration (default: TRUE)
#' 
#' @return qd quantile density values on qdSup
#' 
#' @examples 
#' 
#' x <- seq(0,2,length.out =512)
#' y <- rep(0.5,length.out =512)
#' y.qd <- dens2qd(dens=y, dSup = x) # should equate # 2
#' 
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2016} 
#' @export

dens2qd = function(dens, dSup = seq(0, 1, length.out = length(dens)), qdSup = seq(0, 1, length.out = length(dens)), useSplines = TRUE){

  if(any(dens<=0)){
    stop('Please correct negative or zero probability density estimates.')
  } 
  
  if(!all.equal( range(qdSup),c(0,1) )){
    warning("Problem with support of the QD domain's boundaries - resetting to default.")
    qdSup = seq(0, 1, length.out = length(dens))
  }
  
  if(abs( trapzRcpp(X = dSup, dens) - 1) > 1e-5){
    warning('Density does not integrate to 1 with tolerance of 1e-5 - renormalizing now.')
    dens = dens/trapzRcpp(X = dSup, Y = dens)
  }
  
  if( useSplines ){
    # Could fit spline if this yields more accurate numerical integration
    dens_sp = splinefun(dSup, dens, method = 'natural')
    # Get grid and function for density space
    qdtemp = c(0, cumsum(sapply(2:length(dSup), function(i) integrate(dens_sp, dSup[i - 1], dSup[i])$value)))
  } else {
    # Get grid and function for density space
    qdtemp = cumtrapzRcpp(dSup, dens)
  }

  qdens_temp = 1/dens;
  
  # Remove duplicates
  ind = duplicated(qdtemp)
  qdtemp = unique(qdtemp)
  
  # Interpolate to qdSup
  qd = approx(x = qdtemp, y = qdens_temp[!ind], xout = qdSup, rule = c(2,2))[[2]]
  qd = qd*diff(range(dSup))/trapzRcpp(X = qdSup,Y = qd); # Normalize

  return(qd)
}
