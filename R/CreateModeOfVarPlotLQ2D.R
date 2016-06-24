#' Functional Principal Component Analysis mode of variation plot
#' 
#' Create the k-th mode of variation plot around the mean. The red-line is
#' the functional mean, the grey shaded areas show the range of variations
#' around the mean: \eqn{ \pm Q \sqrt{\lambda_k} \phi_k}{+/- Q sqrt{lambda_k} phi_k}
#' for the dark grey area Q = 1, and for the light grey are Q = 2.
#'
#' @param fpcaObj An FPCA class object returned by FPCA(). 
#' @param k The k-th mode of variation to plot (default k = 1) 
#' @param domain character defining if we should plot on the Quantile ('Q') or the Density ('D') domain (default: 'Q')
#' @param dSup The support of the original density used by LQD(relevant only for density domain)
#' @param dSupPlot The support of the original density used for plotting (relevant only for density domain)
#' @param alpha regularisation parameter
#' @param numOfModes scalar number of principal modes to plot (relevant only for density domain, needs to be an odd number >1.) (default: 7)
#' @param ... Additional arguments for the 'plot' function.
#'
#' @examples
#' library(fdapace)
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$Ly, sampWiener$Lt)
#' CreateModeOfVarPlotLQ2D(res)
#' @export

CreateModeOfVarPlotLQ2D <-function(fpcaObj, dSup = NULL, k = 1, domain = 'Q', alpha = 0, numOfModes = 7, dSupPlot = NULL, ...){  
  
  
  thisMedian = median(1:numOfModes);
  if( !all.equal(thisMedian,  as.integer(thisMedian) ) ){
    stop( "The central tendency is not 0." )
  }
  
  
  if (!any(domain %in% c('D','Q'))){
    stop("Invalid domain option.")
  }
  
  if(domain == 'Q'){
    fdapace::CreateModeOfVarPlot(fpcaObj, k = 1, ...)
    return(invisible(0))
  } else { 
    if(is.null(dSup)){
      stop("Please provide a density domain over which to project the LQ domain.")
    }
    
    if(is.null(dSupPlot)){
      dSupPlot = dSup
    }
    
    args1 <- list( main="Default Title", xlab='s', ylab='', lty=1, lwd= 2)  
    inargs <- list(...)
    args1[names(inargs)] <- inargs
    
    if(k> length(fpcaObj$lambda) ){
      stop("You are asking to plot a mode of variation that is incomputable.")
    }  
    
    obsGrid = fpcaObj$obsGrid      
    s = fpcaObj$workGrid
    mu = fpcaObj$mu
    
    sigma = sqrt(fpcaObj$lambda[k])
    sigma1 = sqrt(fpcaObj$lambda[1])
    phi = fpcaObj$phi[,k]
    phi1 = fpcaObj$phi[,1]
    
    
    Qmatrix = (seq(-3 ,3 ,length.out = numOfModes) %*% t(fpcaObj$phi[,k] * sqrt(fpcaObj$lambda[k]))) + 
      matrix(rep(mu,numOfModes), nrow=numOfModes, byrow = TRUE)
    
    DENreg <- t(apply( Qmatrix, 1, function(u) lqd2dens( u, dSup = dSup, useSplines = TRUE) ) ) 
    
    # B. Correct them according to their bump / Here I do not use alpha but rather the smallest value
    useAlpha = FALSE
    DEN3 <- t(apply( DENreg,  1, function(u) 
      DeregulariseByAlpha(y = u, x = dSup, alpha =  min(u)) ) ) 
    
    DEN4 <- t(apply( DEN3, 1, function(u) 
      approx(y = u, x = dSup, xout = dSupPlot)$y ) ) 
    DEN4[ is.na(DEN4) ] = 0;
     
    thisColPalette = colorRampPalette(c("blue","red", "green"))(numOfModes)
    
    do.call(matplot, c(list(type='l'), list(x=dSupPlot), list(y=t(DEN4)), list(col= thisColPalette), args1))
    grid()    
    lines(x=dSupPlot, y = DEN4[median(1:nrow(DEN4)),], col='red')
    if(numOfModes < 10 ){
      thisMode = (seq(-3 ,3 ,length.out = numOfModes))
      legend("topright", col=thisColPalette, lty=1, border = FALSE, fill=FALSE, lwd= 2,
             legend = do.call(expression, sapply(1:numOfModes, function(u) return(bquote(mu +  ( .(thisMode[u]) * phi[ .(k)]))))) )
    }
    
  }
  
  
}
