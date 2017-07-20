#' Function for converting densities to log hazard functions 
#' 
#' @param dens density values on dSup - must be strictly positive and integrate to 1
#' @param dSup support (grid) for Density domain
#' @param lhSup support for LH domain - must begin at 0 and end at 1; default [0,1] with N-equidistant support points
#' @param N desired number of points on a [0,1] grid for quantile function; default length(dSup)
#' 
#' @return lhz log hazard density on lhSup
#' 
#' @examples
#' x <- seq(0,2,length.out =512)
#' y <- rep(0.5,length.out =512)
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2015} 
#' @export

dens2lhz = function(dens, dSup, N = length(dSup), lhSup = NULL, delta = NULL) {
  
  # Defalut log hazard support
  if(is.null(lhSup)) {
    lhSup = seq(0, 1, length.out = N)
  }

  # Check density requirements
  if(any(dens <= 0)) {
    stop('Please correct negative or zero probability density estimates.')
  }  
  if(abs(fdapace:::trapzRcpp(X = dSup, dens) - 1) > 1e-5) {
    stop('Density does not integrate to 1 with tolerance of 1e-5 - please adjust.')
  }
  if(!all.equal(range(lhSup), c(0, 1))) {
    stop("Please check the support of the LQ domain's boundaries.")
  }
  
  # Get CDF  
  cdf = fdapace:::cumtrapzRcpp(X = dSup, dens)
  if(is.null(delta)) {
    delta = 1 - max(seq(0, 1, length.out = N)[cdf < 1])
  }
  
  # Truncate using delta
  lhSup_use = lhSup[lhSup < 1 - delta]
  
  # Get log hazard and evaluate on lhSup
  lhz_temp = suppressWarnings(log(dens / (1 - cdf)))
  lhz_temp = replace(lhz_temp, is.infinite(lhz_temp), NA)
  lhz = spline(x = seq(0, 1, length.out = N),
               y = lhz_temp, xout = lhSup_use)$y

  # Output a list with domain, function, and delta values
  out = list(lhSup_use, lhz, delta)
  names(out) = c('x', 'y', 'delta')
  return(out)
}

