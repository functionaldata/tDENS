#' Convenience function for converting log quantile densities to densities                                                                                                     
#' 
#' See 'lqd2dens' and 'DeregulariseByAlpha' for more details.
#' This function transforms the log quantile densities in 'qmatrix' to density functions, optionally followed by deregularisation.
#' 
#' @param qmatrix Matrix holding the log quantile density values on [0,1]
#' @param lqdSup Support grid for input log quantile densities (default = seq(0, 1, length.out = ncol(qmatrix)))
#' @param dSup Support grid for output densities (default = seq(0, 1, length.out = ncol(qmatrix)))
#' @param useAlpha Logical indicator to deregularise the densities (default = FALSE)
#' @param alpha Scalar to deregularise the density - where possible, this will be the minimum value for the deregularised densities (default=0)

#' @return list with the 'DEN' transformed data, and 'dSup' that matches the input argument.
#' 
#' @examples  
#' 
#' x <- seq(0,1,length.out = 101)
#' # linear densities on (0, 1)
#' y <- t(sapply(seq(0.5, 1.5, length.out = 10), function(b) b + 2*(1 - b)*x)) 
#' 
#' # Get LQDs
#' 
#' y.lqd = MakeLQDsample(dmatrix = y, dSup = x)
#' matplot(y.lqd$lqdSup, t(y.lqd$LQD), ylab = 'LQD', type = 'l', lty = 1, col = 'black')
#' 
#' # Get Densities Back
#' 
#' y.dens = MakeDENsample(y.lqd$LQD, lqdSup = x, dSup = x) # should equate to y above
#' # These should look the same
#' matplot(y.dens$dSup, t(y.dens$DEN), ylab = 'Density', type = 'l', lty = 1, col = 'blue')
#' matplot(x, t(y), ylab = 'Original Density', type = 'l', lty = 1, col = 'red')
#' 
#' @seealso \code{\link{DeregulariseByAlpha},\link{lqd2dens}}
#' @references
#'   \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2016} 
#' @export


MakeDENsample <- function(qmatrix, lqdSup = seq(0, 1, length.out = ncol(qmatrix)), dSup = seq(0, 1, length.out = ncol(qmatrix)), useAlpha = FALSE,  alpha = 0){
  # A. Get the densities from the log-quantile-density projections
  DEN <- t(apply( qmatrix, 1, function(u) lqd2dens(u, lqdSup = lqdSup, dSup = dSup, useSplines = TRUE) ) ) 
  
  # B. Correct them according to alpha
  if(useAlpha){
    DENtmp <- t(apply( DEN,  1, function(u) 
      DeregulariseByAlpha(y = u, x = dSup, alpha = alpha ) ) ) 
  
    DEN = DENtmp
  }
  
  return( list(DEN = DEN, dSup = dSup) )
}

