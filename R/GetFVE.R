#' Compute Metric-based Fraction of Variance Explained 
#' 
#' When FPCA is performed on the log quantile density functions, the fraction of variance explained by the first K components is computed based on 
#' the density reconstruction and chosen metric.
#' 
#' @param fpcaObj PACE output (FPCA on LQDs)
#' @param dmatrix matrix of original densities measures on grid dSup, rows correspond to individual densities
#' @param dSup support for Density domain - max and min mark the boundary of the support. 
#' @param metric metric for measuring variance - 'L2' for Euclidean or 'W' for Wasserstein
#' @param useAlpha should regularisation be performed to densities in dmatrix? This should be set to TRUE if densities were regularised prior to FPCA (default = FALSE)
#' @param alpha scalar to regularise before computing FVE.  If useAlpha = TRUE, this should match the value used to regularise prior to FPCA (default = 0.01)
#' 
#' @return FVEvector
#' 
#' @details The fraction of variance explained (FVE) by the first K principal components corresponding to the LQD functions is computed by taking the K-dimensional
#' LQD representations, transforming back to densities, and comparing the reconstruction to the original densities using the chosen metric.  If densities were regularised
#' prior to transformation and FPCA, the same regularisation parameters should be used here.
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
#' 
#' fveL2 = GetFVE(fpcaObj = X, dmatrix = t(Top50BabyNames$dens$male), dSup = dSup, useAlpha = TRUE)
#' fveW = GetFVE(fpcaObj = X, dmatrix = t(Top50BabyNames$dens$male), dSup = dSup, 
#'                    metric = 'W', useAlpha = TRUE)
#' 
#' @seealso \code{\link{RegulariseByAlpha},\link{lqd2quantile}}
#' 
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2016} 
#' @export

GetFVE = function(fpcaObj, dmatrix, dSup, metric = 'L2', useAlpha = FALSE, alpha=0.01){
  
  if(!(metric %in% c('L2', 'W'))){
    stop('Unrecognized value for metric input.')
  }
  
  if(useAlpha){
    tmp <-  t(apply(dmatrix, 1, function(u) RegulariseByAlpha(u, x = dSup, alpha = alpha) ))
    dmatrix <- tmp
  }
  
  if(metric == 'L2'){
    mu_dens = colMeans(dmatrix)
    #pointsToKeep <-  which(0 < apply(dmatrix, 2, quantile, 0.995))
  
    vtot = mean(apply(dmatrix, 1, function(u) trapzRcpp(X = dSup, Y = (u - mu_dens)^2)))
    K = length(fpcaObj$lambda)
    FVEs = rep(0, K)
    # print(K)
  
    for(k in 1:K){
      #print(k)
      fittedKQ <- fitted(K=k, fpcaObj);
      fittedKD <- MakeDENsample(fittedKQ, dSup = dSup) # Don't deregularise for FVE computation 
      vK <- mean(apply((dmatrix -  fittedKD$DEN)^2, 1, function(u) trapzRcpp(X = dSup, Y = u)))
      FVEs[k] <- (vtot - vK)/vtot 
    }
  }else {
    
    Qmatrix = t(apply(dmatrix, 1, function(d) dens2quantile(d, dSup)))
    mu_Q = colMeans(Qmatrix)
    
    vtot = mean(apply(Qmatrix, 1, function(u) trapzRcpp(X = dSup, Y = (u - mu_Q)^2)))
    K = length(fpcaObj$lambda)
    FVEs = rep(0, K)
    
    for(k in 1:K){
      fittedK.lqd <- fitted(K=k, fpcaObj);
      fittedKQ <- t(apply(fittedK.lqd, 1, function(y) lqd2quantile(lqd = y, lb = dSup[1])))
      vK <- mean(apply((Qmatrix -  fittedKQ)^2, 1, function(u) trapzRcpp(X = dSup, Y = u)))
      FVEs[k] <- (vtot - vK)/vtot  
    }
    
    
  }
 
  return(FVEs)
}
