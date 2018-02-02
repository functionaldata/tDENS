#' Wasserstein Frechet Mean Computation
#' 
#' Function for computing the Wasserstein Frechet mean through quantile density averaging
#' 
#' @param dmatrix matrix of density values on dSup - must be strictly positive and each row must integrate to 1
#' @param dSup support (grid) for Density domain
#' @param qdSup support for LQ domain - must begin at 0 and end at 1; default [0,1] with N-equidistant support points
#' @param N desired number of points on a [0,1] grid for quantile density functions; default length(dSup)
#' @param useAlpha should regularisation be performed (default=FALSE) 
#' @param alpha Scalar to regularise the supports with (default=0.01)
#' 
#' @return wfmean the Wasserstein-Frechet mean density
#' 
#' @examples
#' x <- seq(0,1,length.out = 101)
#' # linear densities on (0, 1)
#' y <- t(sapply(seq(0.5, 1.5, length.out = 10), function(b) b + 2*(1 - b)*x)) 
#' wfmean = getWFmean(y, x)
#' 
#' # Plot WF mean with Euclidean Mean
#' matplot(x, t(y), ylab = 'Density', type = 'l', lty = 1, col = 'black')
#' lines(x, wfmean, lwd = 2, col = 'red')
#' lines(x, colMeans(y), lwd = 2, col = 'blue')
#' legend('topright', col = c('black', 'red', 'blue'), lwd = c(1, 2, 2), 
#'            legend = c('Densities', 'WF Mean', 'Euclidean Mean'))
#' 
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2016} 
#' @export
#' 

getWFmean = function(dmatrix, dSup, N = length(dSup), qdSup = seq(0, 1, length.out = N), useAlpha = FALSE, alpha = 0.01){
  
  if(useAlpha){
    tmp <-  t(apply(dmatrix, 1, function(u) RegulariseByAlpha(u, x = dSup, alpha = alpha) ))
    dmatrix <- tmp
  }
  
  if(any(apply(dmatrix, 1, function(d) abs(trapzRcpp(dSup, d) - 1) > 1e-5))){
    stop('Densities must all integrate to 1.')
  }
  
  qd = t(apply(dmatrix, 1, function(d) dens2qd(d, dSup = dSup, qdSup = qdSup)))
  
  wfmean = qd2dens(colMeans(qd), dSup = dSup, qdSup = qdSup)
  
  return(wfmean)
}
