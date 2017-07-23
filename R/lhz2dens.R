#' Function for converting log quantile densities to densities
#' 
#' @param  lhz log hazard density on lhSup
#' @param  lhSup support for LH domain - must begin at 0 and end at 1 - delta
#' @param  dSup support for Density domain - max and min values are assumed to mark the boundary of the support.
#' @param  useSplines fit spline to the lqd when doing the numerical integration (default: TRUE)
#' @param  delta support trimming variable
#' 
#' @return dens density domain and values on dSup
#' 
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2015} 
#' @export

lhz2dens = function(lhz, lhSup, dSup, useSplines = TRUE, delta) {
  
  # Make subdomains based on delta
  endDensDomain = min(dSup) + (1 - delta) * (diff(range(dSup)))
  densDomain1 = dSup[dSup <= endDensDomain]
  densDomain2 = dSup[dSup > endDensDomain]
  
  # Make adjusted density domain the same length as log hazard domain
  dSup1_adj = seq(min(densDomain1), max(densDomain1), length.out = length(lhz))
  dSup2_adj = seq(min(densDomain2), max(densDomain2), length.out = length(lhz))
  
  # Calculate density on first subdomain
  dens_temp1 = exp(lhz - fdapace:::cumtrapzRcpp(X = dSup1_adj, Y = exp(lhz)))
  dens1 = spline(dSup1_adj, dens_temp1, xout = densDomain1)
  
  # Calculate density on second subdomain
  dens_temp2 = (max(dSup) - endDensDomain)^(-1) * exp(-fdapace:::trapzRcpp(X = dSup1_adj, Y = exp(lhz)))
  dens2 = rep(dens_temp2, length(densDomain2))
  
  # Get density and normalize
  out = list(x = dSup,
             y = c(dens1$y, dens2))
  out$y = out$y / fdapace:::trapzRcpp(out$x, out$y)
  
  return(out)
}

