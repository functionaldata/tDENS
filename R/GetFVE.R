#' Function for converting densities to quantile functions 
#' 
#' @param fpcaObj PACE output (FPCA on LQDs)
#' @param dmatrix matrix of densities measures on grid dout, rows correspond to individual densities
#' @param dSup support for Density domain - max and min values are assumed to mark the boundary of the support. 
#' @param alpha scalar to regularise with (default = 0.01)
#' 
#' @return FVEvector
#' 
#' @examples 
#' pts <- seq(0, 1, by=0.05) 
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2015} 
#' @export

GetFVE = function(fpcaObj, dmatrix, dSup, alpha=0.01){
  
  mu_dens = colMeans(dmatrix)
  pointsToKeep <-  which(0 < apply(dmatrix, 2, quantile, 0.995))
  
  vtot = mean(apply(dmatrix, 1, function(u) fdapace:::trapzRcpp(X = dSup, Y = (u - mu_dens)^2)))
  K = length(fpcaObj$lambda)
  FVEs = rep(0, K)
  # print(K)
  
  for(k in 1:K){
    print(k)
    fittedKQ <- fdapace:::fitted.FPCA(k=k, fpcaObj);
    fittedKD <- MakeDENsample(fittedKQ, alpha = alpha, dSupUsedInLQ = dSup[pointsToKeep], dSup, useAlpha = FALSE ) 
    vK <- mean(apply((dmatrix -  fittedKD$DEN)^2, 1, function(u) fdapace:::trapzRcpp(X = dSup, Y = u)))
    FVEs[k] <- (vtot - vK)/vtot 
  }
 
  return(FVEs)
}
