#' Convenience function for converting log-quantile functions to densities                                                                                                     
#' 
#' See 'lqd2dens' and 'RegulariseByAlpha' for more details.
#' This function will automatically calculate the:sout support points that the 99.5-th percentile (calculated using all available samples) is not positive.
#' It will then normalise the new densities to integraate to one, regularise by 'alpha' and calculate the LQ projection.
#' 
#' @param qmatrix Matrix holding the log-quantiles values on [0,1] and assumed to lay strictly on dSupUsedInLQ
#' @param dSupUsedInLQ Support (grid) for Density domain that was used when calculating the LQ projections (maybe different from the original)
#' @param dSup Original support (grid) for Density domain of the data
#' @param useAlpha Logical indicator to regularise the support so the smallest value of each density is zero 

#' @return list with the 'DEN' projected data, and 'dSup' that was originally used for the data.
#' 
#' @references
#'   \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2015} 
#' @export


MakeDENsample <- function(qmatrix, alpha = 0, dSupUsedInLQ, dSup, useAlpha = FALSE ){
  # A. Get the densities from the log-quantile-density projections
  DENreg <- t(apply( qmatrix, 1, function(u) lqd2dens( u, dSup = dSupUsedInLQ, useSplines = TRUE) ) ) 
  
  # B. Correct them according to their bump / Here I do not use alpha but rather the smallest value
  DEN3 <- t(apply( DENreg,  1, function(u) 
    DeregulariseByAlpha(y = u, x = dSupUsedInLQ, alpha = ifelse( useAlpha, alpha, min(u)) ) ) ) 
  
  # C. Pad the density with the appropariate 0.
  DEN <- matrix(0, nrow = nrow(qmatrix), ncol = length(dSup) )
  DEN[,  match(dSupUsedInLQ, dSup)] <- DEN3
  
  return( list(DEN = DEN, origSupport = dSup) )
}

