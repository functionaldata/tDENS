#' Function for converting log quantile densities to densities
#' 
#' @param  lqd log quantile density on lqSup
#' @param  lqSup support for LQ domain - must begin at 0 and end at 1
#' @param  dSup support for Density domain - max and min values are assumed to mark the boundary of the support.
#' @param  useSplines fit spline to the lqd when doing the numerical integration (default: TRUE)
#' 
#' @return dens density values on dSup
#' 
#' @examples
#' x <- seq(0,2,length.out =512)
#' y <- rep(0.5,length.out =512)
#' yOnLQ <- dens2lqd(dens=y, dSup = x) # should equate # -log(1/2)
#' yOnDens <- lqd2dens(dSup=x, lqd = yOnLQ)   
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2015} 
#' @export

lqd2dens = function(lqd, lqSup = seq(0, 1, length.out = length(lqd)), dSup, useSplines = TRUE){

  if(!all.equal( range(lqSup),c(0,1) )){
    stop("Please check the support of the LQ domain's boundaries.")
  }
  
  if( useSplines ){
    # Could fit spline if this yields more accurate numerical integration
    lqd_sp = splinefun(lqSup, lqd, method = 'natural')
    lqd_exp = function(t) exp(lqd_sp(t))
    # Get grid and function for density space
    dtemp = dSup[1] + c(0, cumsum(sapply(2:length(lqSup), function(i) integrate(lqd_exp, lqSup[i - 1], lqSup[i])$value)))
  } else {
    # Get grid and function for density space
    dtemp = dSup[1] + fdapace:::cumtrapzRcpp(lqSup, exp(lqd))
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
  dens = dens/fdapace:::trapzRcpp(X = dSup,Y = dens); # Normalize

  return(dens)
}

#MakeDENsample <- function(Qmatrix, alpha = 0, dSupUsedInLQ, dSup, useAlpha = FALSE ){
#  # A. Get the densities from the log-quantile-density projections
#  DENreg <- t(apply( Qmatrix, 1, function(u) lqd2dens( u, dSup = dSupUsedInLQ, useSplines = TRUE) ) ) 
#  
#  # B. Correct them according to their bump / Here I do not use alpha but rather the smallest value
#  DEN3 <- t(apply( DENreg,  1, function(u) 
#    DeregulariseByAlpha(y = u, x = dSupUsedInLQ, alpha = ifelse( useAlpha, alpha, min(u)) ) ) ) 
#  
#  # C. Pad the density with the appropariate 0.
#  DEN <- matrix(0, nrow = nrow(Qmatrix), ncol = length(dSup) )
#  DEN[,  match(dSupUsedInLQ, dSup)] <- DEN3
#  
#  return( list(DEN = DEN, origSupport = dSup) )
#}
#
