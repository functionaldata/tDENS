#' Convenience function for converting densities to log-quantile densities
#' 
#' See 'dens2lqd' and 'RegulariseByAlpha' for more details.
#' This function first (transforms the densities in 'dmatrix' to log quantile density functions, optionally followed by regularisation.
#' 
#' @param dmatrix Matrix holding the density values on dSup - all rows must be strictly positive and integrate to 1
#' @param dSup Support (grid) for Density domain
#' @param lqdSup Support grid for lqd domain (default = seq(0, 1, length.out = length(dSup)))
#' @param useAlpha should regularisation be performed (default=FALSE) 
#' @param alpha Scalar to regularise the supports with (default=0.01)
#' 
#' @return list with 'LQD', a matrix of log quantile density functions, and 'lqdSup' that matches the input argument
#' 
#' @examples  
#' 
#' x <- seq(0,1,length.out = 101)
#' y <- t(sapply(seq(0.5, 1.5, length.out = 10), function(b) -log(b^2 + 4*(1-b)*x)/2)) # some log quantile densities on (0, 1)
#' 
#' # Get densities
#' 
#' y.dens = MakeDENsample(qmatrix = y, lqdSup = x, dSup = x)$DEN
#' matplot(x, t(y.dens), ylab = 'Density', type = 'l', lty = 1, col = 'black')
#' 
#' # Get LQDs Back
#' 
#' y.lqd = MakeLQDsample(y.dens, lqdSup = x, dSup = x)
#' # These should match
#' matplot(y.lqd$lqdSup, t(y.lqd$LQD), ylab = 'LQD', type = 'l', lty = 1, col = 'blue')
#' matplot(x, t(y), ylab = 'LQD', type = 'l', lty = 1, col = 'red')
#' 
#' @seealso \code{\link{RegulariseByAlpha},\link{dens2lqd}} 
#' @references
#'   \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2016} 
#' @export

MakeLQDsample <- function(dmatrix, dSup, lqdSup = seq(0, 1, length.out = length(dSup)), useAlpha = FALSE, alpha = 0.01){

  #pointsToKeep <-  which(0 < apply(dmatrix, 2, quantile, 0.995))
  #newSupport <- dInput[pointsToKeep]
  #dmatrix2 <- dmatrix[, pointsToKeep];
  #dmatrix3 <- dmatrix2 / apply(dmatrix2, 1, function(u) fdapace:::trapzRcpp(Y=  u, X = newSupport) )
  
  if(useAlpha){
    tmp <-  t(apply(dmatrix, 1, function(u) RegulariseByAlpha(u, x = dSup, alpha = alpha) ))
    dmatrix <- tmp
  }
  
  LQD <- t(apply( dmatrix, 1, dens2lqd, dSup = dSup, lqdSup = lqdSup) ) 
  return(list(LQD = LQD, lqdSup = lqdSup))
}

