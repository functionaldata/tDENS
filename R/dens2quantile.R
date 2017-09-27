#' Function for converting Densities to Quantile Functions
#' 
#' @param  dens density on dSup
#' @param  dSup support for Density domain - max and min values mark the boundary of the support.
#' @param  qSup support for quantile domain - must begin at 0 and end at 1
#' @param  useSplines fit spline to the qd when doing the numerical integration (default: TRUE)
#' 
#' @return Q quantile function on qSup
#' 
#' @examples 
#' 
#' x <- seq(0,2,length.out =512)
#' y <- rep(0.5,length.out =512)
#' y.quantile <- dens2quantile(dens=y, dSup = x) # should equate # 2*seq(0, 1, length.out = 512)
#' 
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2016} 
#' @export

dens2quantile = function(dens, dSup = seq(0, 1, length.out = length(dens)), qSup = seq(0, 1, length.out = length(dens)), useSplines = TRUE){

  if(!all.equal( range(qSup),c(0,1) )){
    stop("Please check the support of the Q domain's boundaries.")
  }
  
  if(abs(trapzRcpp(dSup, dens) - 1) > 1e-5){
    stop("Density should integrate to one.")
  }
  
  if( useSplines ){
    # Could fit spline if this yields more accurate numerical integration
    dens_sp = splinefun(dSup, dens, method = 'natural')
    # Get grid and function for density space
    QSuptmp = c(0, cumsum(sapply(2:length(dSup), function(i) integrate(dens_sp, dSup[i - 1], dSup[i])$value)))
  } else {
    # Get grid and function for density space
    QSuptmp = cumtrapzRcpp(dSup, dens)
  }

  Qtmp = dSup
  
  # Remove duplicates
  ind = duplicated(QSuptmp)
  QSuptmp = unique(QSuptmp)
  
  # Interpolate to qdSup
  Q = approx(x = QSuptmp, y = Qtmp[!ind], xout = qSup, rule = c(2,2))[[2]]

  return(Q)
}
