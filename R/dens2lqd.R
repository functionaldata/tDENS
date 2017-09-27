#' Function for converting densities to log quantile density functions 
#' 
#' @param dens density values on dSup - must be strictly positive and integrate to 1
#' @param dSup support (grid) for Density domain
#' @param lqdSup support for lqd domain - must begin at 0 and end at 1; default [0,1] with N-equidistant support points
#' @param N desired number of points on a [0,1] grid for lqd function; default length(dSup)
#' 
#' @return lqd log quantile density on lqdSup
#' 
#' @examples
#' x <- seq(0,2,length.out =512)
#' y <- rep(0.5,length.out =512)
#' y.lqd <- dens2lqd( dens=y, dSup = x) # should equate # log(2)
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2016} 
#' @export

dens2lqd = function(dens, dSup, N = length(dSup), lqdSup = NULL){

  if(is.null(lqdSup)){
      lqdSup = seq(0, 1, length.out = N)
  }

  # Check density requirements
  if(any(dens<=0)){
    stop('Please correct negative or zero probability density estimates.')
  }  
  if(abs( fdapace:::trapzRcpp(X = dSup, dens) - 1) > 1e-5){
    stop('Density does not integrate to 1 with tolerance of 1e-5 - please adjust.')
  }
  if(!all.equal( range(lqdSup),c(0,1) )){
    stop("Please check the support of the LQD domain's boundaries.")
  }
 
  # Get CDF  
  qtemp = fdapace:::cumtrapzRcpp(X = dSup, dens)
  # ind = duplicated(qtemp)
  # qtemp = unique(qtemp)
  # lqd_temp = -log(dens[!ind]);
  lqd_temp = -log(dens)
  lqd = spline(x = qtemp, y = lqd_temp, xout = lqdSup)$y
  return(lqd)
}