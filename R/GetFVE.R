#' Function for converting densities to quantile functions 
#' 
#' @param fpcaObj PACE output (FPCA on LQDs)
#' @param dmatrix matrix of original densities measures on grid dout, rows correspond to individual densities
#' @param dSup support for Density domain - max and min mark the boundary of the support. 
#' @param metric metric for measuring variance - 'L2' for Euclidean or 'W' for Wasserstein
#' @param useAlpha should deregularisation be performed? (default = FALSE)
#' @param alpha scalar to deregularise when transforming back to density space with (default = 0)
#' 
#' @return FVEvector
#' 
#' @examples 
#' 
#' data(Top50BabyNames)
#'
#' # Perform Transformation FPCA for male baby name densities
#' dSup = Top50BabyNames$x
#' X = FPCAdens(dmatrix = t(Top50BabyNames$dens$male), dSup = dSup, useAlpha = TRUE, 
#'                    optns = list(dataType = 'Dense', error = FALSE, methodSelectK = 8))
#'
#' # Compute FVE - must compare to regularized densities 
#' dens.reg = t(apply(Top50BabyNames$dens$male, 2, function(y) RegulariseByAlpha(dSup, y)))
#' 
#' fve.L2 = GetFVE(fpcaObj = X, dmatrix = dens.reg, dSup = dSup)
#' fve.W = GetFVE(fpcaObj = X, dmatrix = dens.reg, dSup = dSup, metric = 'W')
#' 
#' @seealso \code{\link{lqd2dens},\link{DeregulariseByAlpha}}
#' 
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2016} 
#' @export

GetFVE = function(fpcaObj, dmatrix, dSup, metric = 'L2', useAlpha = FALSE, alpha=0){
  
  if(!(metric %in% c('L2', 'W'))){
    stop('Unrecognized value for metric input.')
  }
  
  if(metric == 'L2'){
    mu_dens = colMeans(dmatrix)
    #pointsToKeep <-  which(0 < apply(dmatrix, 2, quantile, 0.995))
  
    vtot = mean(apply(dmatrix, 1, function(u) fdapace:::trapzRcpp(X = dSup, Y = (u - mu_dens)^2)))
    K = length(fpcaObj$lambda)
    FVEs = rep(0, K)
    # print(K)
  
    for(k in 1:K){
      #print(k)
      fittedKQ <- fdapace:::fitted.FPCA(K=k, fpcaObj);
      fittedKD <- MakeDENsample(fittedKQ, alpha = alpha, dSup = dSup, useAlpha = useAlpha) 
      vK <- mean(apply((dmatrix -  fittedKD$DEN)^2, 1, function(u) fdapace:::trapzRcpp(X = dSup, Y = u)))
      FVEs[k] <- (vtot - vK)/vtot 
    }
  }else {
    
    Qmatrix = t(apply(dmatrix, 1, function(d) dens2quantile(d, dSup)))
    mu_Q = colMeans(Qmatrix)
    
    vtot = mean(apply(Qmatrix, 1, function(u) fdapace:::trapzRcpp(X = dSup, Y = (u - mu_Q)^2)))
    K = length(fpcaObj$lambda)
    FVEs = rep(0, K)
    
    for(k in 1:K){
      fittedKQ <- fdapace:::fitted.FPCA(K=k, fpcaObj);
      fittedKD <- MakeDENsample(fittedKQ, alpha = alpha, dSup = dSup, useAlpha = useAlpha)
      fittedKQ2 <- t(apply(fittedKD$DEN, 1, function(d) dens2quantile(d, dSup)))
      vK <- mean(apply((Qmatrix -  fittedKQ2)^2, 1, function(u) fdapace:::trapzRcpp(X = dSup, Y = u)))
      FVEs[k] <- (vtot - vK)/vtot  
    }
    
    
  }
 
  return(FVEs)
}
