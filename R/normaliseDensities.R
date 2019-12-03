#' Normalise Densities
#' 
#' Preprocessing function to ensure densities integrate to 1  
#'
#' @param dmatrix Matrix with rows representing distinct densities on dSup - all entries must be nonnegative
#' @param dSup Support (grid) for Density domain
#'
#' @details Uses trapezoidal integration to normalise the densities to have integral 1
#'  
#' @examples
#' ## Normalise collection of truncated normal densities
#' 
#' mu <- seq(-2, 2, by = 0.5)
#' dSup = seq(-3, 3, length.out = 101)
#' y <- t(sapply(mu, function(m) dnorm(x = dSup, mean = m)))
#' 
#' # Check integral constraint
#' apply(y, 2, function(yy) trapzRcpp(X = dSup, Y = yy))
#' 
#' # Normalise and check again
#' dens <- normaliseDensities(dmatrix = y, dSup = dSup)
#' apply(dens, 2, function(yy) trapzRcpp(X = dSup, Y = yy))
#' 
#' 
#' @return matrix 'dmatrix' consisting of rows of input of the same name that have been normalised to have integral 1
#' 
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2016} 

#' @export

normaliseDensities <-function(dmatrix, dSup = 1:ncol(dmatrix)){  
  
  if (any(dmatrix < 0)){
    stop("Density values must be nonnegative.")
  }
  
  dmatrix = t(apply(dmatrix, 1, function(d) d/trapzRcpp(X = dSup, Y = d)))
  
  return(dmatrix)
   
}