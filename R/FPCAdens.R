#' FPCA for densities by log quantile density transformation
#' 
#' Perform FPCA on LQD-transformed densities 
#'
#' @param dmatrix Matrix holding the density values on dSup - all rows must be strictly positive and integrate to 1
#' @param dSup Support (grid) for Density domain
#' @param lqdSup Support grid for lqd domain (default = seq(0, 1, length.out = length(dSup)))
#' @param useAlpha should regularisation be performed (default=FALSE) 
#' @param alpha Scalar to regularise the supports with (default=0.01)
#' @param optns A list of options for FPCA.  See documentation for \code{FPCA}.
#'
#' @details Densities are transformed to log-quantile densities, followed by standard FPCA.  If \code{useAlpha = TRUE}, densities are regularized before transformation
#' 
#' @seealso \code{\link{RegulariseByAlpha},\link{lqd2dens},\link{MakeLQDsample},\link{FPCA}}
#' 
#' @examples
#' 
#' ## Densities for Top 50 Female Baby Names
#' data(Top50BabyNames)
#'
#' # Perform Transformation FPCA for male baby name densities
#' X = FPCAdens(dmatrix = t(Top50BabyNames$dens$female), dSup = Top50BabyNames$x, useAlpha = TRUE, 
#'                   optns = list(dataType = 'Dense', error = FALSE, methodSelectK = 2))
#' x = Top50BabyNames$x
#' 
#' # Plot Modes
#' 
#' Qvec = quantile(X$xiEst[,1], probs = c(0.1, 0.25, 0.75, 0.9))/sqrt(X$lambda[1])
#' CreateModeOfVarPlotLQ2D(X, k = 1, dSup = x, Qvec = Qvec, main = 'First Mode, Density space')
#' CreateModeOfVarPlotLQ2D(X, domain = 'Q', k = 1, dSup = x, Qvec = Qvec, 
#'                             main = 'First Mode, LQD space')
#'
#' Qvec = quantile(X$xiEst[,2], probs = c(0.1, 0.25, 0.75, 0.9))/sqrt(X$lambda[2])
#' CreateModeOfVarPlotLQ2D(X, k = 2, dSup = x, Qvec = Qvec, main = 'Second Mode, Density Space')
#' CreateModeOfVarPlotLQ2D(X, domain = 'Q', k = 2, dSup = x, Qvec = Qvec, 
#'                             main = 'Second Mode, LQD space')
#' 
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2016} 

#' @export

FPCAdens = function(dmatrix, dSup, lqdSup = seq(0, 1, length.out = length(dSup)), useAlpha = F, alpha = 0.01, optns = list(dataType = 'Dense', error = FALSE)){
  
  # Transform Densities
  
    lqd = MakeLQDsample(dmatrix=dmatrix,dSup=dSup,lqdSup=lqdSup,useAlpha=useAlpha,alpha=alpha)$LQD
    
  # Perform FPCA
    
    L = fdapace::MakeFPCAInputs(tVec = lqdSup, yVec = lqd)
    X = fdapace::FPCA(Ly = L$Ly, Lt = L$Lt, optns = optns)
    
    return(X)
    
}