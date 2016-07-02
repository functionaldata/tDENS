#' Function for converting densities to quantile functions 
#' 
#' @param dens density values on dSup - must be strictly positive and integrate to 1
#' @param dSup support (grid) for Density domain
#' @param lqSup support for LQ domain - must begin at 0 and end at 1; default [0,1] with N-equidistant support points
#' @param N desired number of points on a [0,1] grid for quantile function; default length(dSup)
#' 
#' @return lqd log quantile density on lqSup
#' 
#' @examples
#' x <- seq(0,2,length.out =512)
#' y <- rep(0.5,length.out =512)
#' yOnLQ <- dens2lqd( dens=y, dSup = x) # should equate # -log(1/2)
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2015} 
#' @export

dens2lqd = function(dens, dSup, N = length(dSup), lqSup = NULL){

  if(is.null(lqSup)){
      lqSup = seq(0, 1, length.out = N)
  }

  # Check density requirements
  if(any(dens<=0)){
    stop('Please correct negative or zero probability density estimates.')
  }  
  if(abs( fdapace:::trapzRcpp(X = dSup, dens) - 1) > 1e-5){
    stop('Density does not integrate to 1 with tolerance of 1e-5 - please adjust.')
  }
  if(!all.equal( range(lqSup),c(0,1) )){
    stop("Please check the support of the LQ domain's boundaries.")
  }
 
  # Get CDF  
  qtemp = fdapace:::cumtrapzRcpp(X = dSup, dens)
  # ind = duplicated(qtemp)
  # qtemp = unique(qtemp)
  # lqd_temp = -log(dens[!ind]);
  lqd_temp = -log(dens)
  lqd = spline(x = qtemp, y = lqd_temp, xout = lqSup)$y
  return(lqd)
}

#MakeLQDsample <- function(dmatrix, alpha = 0, dInput ){
#  # A. Find the points to keep and make sure they integrate to 1.
#  pointsToKeep <-  which(0 < apply(dmatrix, 2, quantile, 0.995))
#  newSupport <- dInput[pointsToKeep]
#  dmatrix2 <- dmatrix[, pointsToKeep];
#  dmatrix3 <- dmatrix2 / apply(dmatrix2, 1, function(u) fdapace:::trapzRcpp(Y=  u, X = newSupport) )
#  
#  # B. Regularise by alpha
#  dmatrixReg <-  t(apply( dmatrix3, 1, function(u) RegulariseByAlpha(u, x = newSupport, alpha = alpha, deregularise = FALSE) ))
#  
#  # C. Get the log-quantile-density projections
#  LQD <- t(apply( dmatrixReg, 1, dens2lqd, dSup = dInput[pointsToKeep]) ) 
#  return( list(LQD = LQD, alpha = alpha, newSupport = newSupport) )
#}

