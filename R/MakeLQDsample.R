#' Convenience function for converting densities to log-quantile functions
#' 
#' See 'dens2lqd' and 'RegulariseByAlpha' for more details.
#' This function will automatically trim out support points that the 99.5-th percentile (calculated using all available samples) is not positive.
#' It will then normalise the new densities to integraate to one, regularise by 'alpha' and calculate the LQ projection.
#' 
#' @param dmatrix Matrix holding the density values on dInput - assumed to be strictly positive and integrate to 1
#' @param dInput Support (grid) for Density domain
#' @param alpha Scalar to regularise the supports with (default=0)
#' 
#' @return list with the 'LQD' projected data, the 'alpha' used and the 'newSupport' that was defined after trimming the data.
#' 
#' @references
#'   \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2015} 
#' @export

MakeLQDsample <- function(dmatrix, alpha = 0, dInput ){
  # A. Find the points to keep and make sure they integrate to 1.
  pointsToKeep <-  which(0 < apply(dmatrix, 2, quantile, 0.995))
  newSupport <- dInput[pointsToKeep]
  dmatrix2 <- dmatrix[, pointsToKeep];
  dmatrix3 <- dmatrix2 / apply(dmatrix2, 1, function(u) fdapace:::trapzRcpp(Y=  u, X = newSupport) )
  
  # B. Regularise by alpha
  dmatrixReg <-  t(apply( dmatrix3, 1, function(u) RegulariseByAlpha(u, x = newSupport, alpha = alpha, deregularise = FALSE) ))
  
  # C. Get the log-quantile-density projections
  LQD <- t(apply( dmatrixReg, 1, dens2lqd, dSup = dInput[pointsToKeep]) ) 
  return( list(LQD = LQD, alpha = alpha, newSupport = newSupport) )
}

