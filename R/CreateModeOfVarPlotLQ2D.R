#' Transformation Mode of Variation Plot
#' 
#' Create the k-th transformation mode of variation plot.  
#'
#' @param fpcaObj An FPCA class object returned by FPCA() on the log quantile density functions. 
#' @param domain should the mode be plotted in LQD ('Q') or density space ('D', the default).
#' @param k The k-th mode of variation to plot (default k = 1)
#' @param dSup The common support of the original densities. Only relevant for \code{domain = 'D'}
#' @param Qvec Vector of values \eqn{Q}{Q} to be plotted.  If 0 is not included, it will be added (default is -2:2).  Only relevant for \code{domain = 'D'}
#' @param alpha (De)regularisation parameter (default is 0).  See details.
#' @param useAlpha logical - should deregularisation be performed?  Default:FALSE
#' @param ... Additional arguments for the 'plot' function.
#'
#' @details If \code{domain = 'D'} (the default), the a transformation mode of variation is plotted.  The red-line is 
#' \eqn{\psi^{-1}(\nu)}{psi^(-1)(nu)}, where \eqn{\nu}{nu} is the mean in LQD space and
#' \eqn{\psi}{psi} is the LQD transformation.  Other lines correspond to perturbations by 
#' adding multiples of the LQD eigenfunctions \eqn{\rho_k}{rho_k} (with eigenvalues \eqn{\tau_k}{tau_k})
#' : \eqn{\psi^{-1}(\nu  + Q \sqrt{\tau_k} \rho_k)}{psi^(-1)(nu + Q sqrt{tau_k} rho_k} 
#' for the values \eqn{Q}{Q} in \code{Qvec}.  If \code{alpha} is positive, will attempt to deregularise (see \code{DeregulariseByAlpha}). This
#' will throw an error if alpha is too large.
#'
#' If \code{domain = 'Q'}, ordinary modes of variation are plotted in LQD space (see documentation for \code{CreateModeOfVarPlot} in \code{fdapace}).
#' 
#' @seealso \code{\link{DeregulariseByAlpha}}
#' 
#' @examples
#' ## Densities for Top 50 Male Baby Names
#' data(Top50BabyNames)
#' x = Top50BabyNames$x
#'
#' # Perform Transformation FPCA for male baby name densities
#' X = FPCAdens(dmatrix = t(Top50BabyNames$dens$male), dSup = Top50BabyNames$x, useAlpha = TRUE, 
#'                   optns = list(dataType = 'Dense', error = FALSE, methodSelectK = 2))
#' 
#' # Plot Modes
#' 
#' Qvec = quantile(X$xiEst[,1], probs = c(0.1, 0.25, 0.75, 0.9))/sqrt(X$lambda[1])
#' CreateModeOfVarPlotLQ2D(X, k = 1, dSup = x, Qvec = Qvec, main = 'First Mode, Density Space')
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

CreateModeOfVarPlotLQ2D <-function(fpcaObj, domain = 'D', k = 1, dSup = NULL, Qvec = -2:2, alpha = 0, useAlpha = FALSE, ...){  
  
  if (!any(domain %in% c('D','Q'))){
    stop("Invalid domain option.")
  }
  
  
  if(domain == 'Q'){
    fdapace::CreateModeOfVarPlot(fpcaObj, k = k, ...)
    return(invisible(0))
  } else {
    
    if(is.null(dSup)){
      stop("Please provide a density domain for inverse transformation.")
    }
  
    if(!(0 %in% Qvec)){
      Qvec = sort(c(0, Qvec))
    }
    
    args1 <- list( main="Default Title", xlab='s', ylab='', lty=1, lwd= 2)  
    inargs <- list(...)
    args1[names(inargs)] <- inargs
    
    if(k > length(fpcaObj$lambda) ){
      stop("The FPCA object has fewer components than the requested value of k.")
    }  
    
    obsGrid = fpcaObj$obsGrid      
    s = fpcaObj$workGrid
    nu = fpcaObj$mu
    
    sigma = sqrt(fpcaObj$lambda[k])
    rho = fpcaObj$phi[,k]
    
    Qmatrix = Qvec %*% t(rho * sigma) + 
      matrix(rep(nu,length(Qvec)), nrow=length(Qvec), byrow = TRUE)
    
    DENreg <- t(apply( Qmatrix, 1, function(u) lqd2dens(u, dSup = dSup, useSplines = TRUE) ) )
    
    if(useAlpha){
      DENtmp <- t(apply( DENreg,  1, function(u)
        DeregulariseByAlpha(y = u, x = dSup, alpha = alpha) ) )
      
      DENreg <- t(apply( DENtmp, 1, function(u) 
        approx(y = u, x = dSup, xout = dSup)$y ) ) 
      DENreg[ is.na(DENreg) ] = 0;
    }
     
    thisColPalette = colorRampPalette(c("blue","red", "green"))(length(Qvec))
    
    do.call(matplot, c(list(type='l'), list(x=dSup), list(y=t(DENreg)), list(col= thisColPalette), args1))
    grid()    
    lines(x=dSup, y = DENreg[Qvec == 0,], col='red')
    if(length(Qvec) < 8 ){
      thisMode = Qvec
      legend("topright", col=thisColPalette, lty=1, border = FALSE, fill=FALSE, lwd= 2,
             legend = do.call(expression, sapply(1:length(Qvec), function(u){ 
               
               if(Qvec[u] == 0){
                 return(bquote(psi^-1*(nu)))
               } else if(Qvec[u] < 0){
                 return(bquote(psi^-1*(nu - .(round(abs(thisMode[u])*sqrt(sigma), 2)) * rho[ .(k)])))
               } else{
                 return(bquote(psi^-1*(nu + .(round(thisMode[u]*sqrt(sigma), 2)) * rho[ .(k)])))
               }
               })) )
    }
  }
}
