#' Create density from raw data
#' 
#' Create kernel density estimate along the support of the raw data using the HADES method
#' 
#' @param y A vector of raw readings. 
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are 
#' \describe{
#' \item{userBwMu}{The bandwidth value for the smoothed mean function (using 'CV' or 'GCV'); positive numeric - default: determine automatically based on CV}
#' \item{nRegGrid}{The number of support points the KDE; numeric - default: 101}
#' \item{delta}{The size of the bin to be used; numeric - default: determine automatically as "max(c(diff(range(y))/1000, min(diff(sort(unique(y))))))"}
#' \item{kernel}{smoothing kernel choice, "rect", "gauss", "epan", "gausvar", "quar" - default: "gauss"}
#' \item{infSupport}{logical if we expect the distribution to have infinite support or not; logical - default: TRUE}
#' \item{outputGrid}{User defined output grid for the support of the KDE, it overrides nRegGrid; numeric - default: NULL}
#' }
#' 
#' @return A list containing the following fields:
#' \item{bw}{Variance for measure error.The bandwidth used by smoothing.}
#' \item{x}{A vector of length \emph{nGridReg} with the values of the KDE's support points.}
#' \item{y}{A vector of length \emph{nGridReg} with the values of the KDE at the support points.} 
#' 
#' @examples
#' par(mfrow=c(1,2))
#' makeComparisonPlotE <- function(N, mySeed = 123){
#'   set.seed(mySeed)
#'   asdf2  = (rexp(N, rate = 1.5))
#'   plot(density(asdf2, bw = "SJ"), main= "Exponential (rate=1.5)", 
#'     xlab = paste0(collapse = '', c( "N = ", as.character(N))) , ylim = c(0, 1.45))
#'   lines(density(asdf2), col='red')
#'   Unormal = CreateDensity(y = asdf2, optns = list(infSupport = FALSE))
#'   lines(col='green', x = Unormal$x, y = Unormal$y)
#'   lines(col='magenta' , Unormal$x, dexp(Unormal$x, rate = 1.5))
#'   abline(v = min(asdf2))
#'   abline(v = max(asdf2))
#'   legend(legend = c("SJ", "R-default", "HADES-like", "True PDF"),
#'     lwd= 2, col=c("black", "red", "green", "magenta"), bty = 'n', 'topright')
#' }
#' 
#' makeComparisonPlotE(100)
#' makeComparisonPlotE(2000)
#' 
#' @references
#' \cite{HG Mueller, JL Wang and WB Capra (1997). "From lifetables to hazard rates: The transformation approach." Biometrika 84, 881-892.}
#' @export



CreateDensity <- function(y, optns = list()){
  
  if(is.null(optns$kernel)){
    kernel = 'gauss'
  } else {
    kernel =  optns$kernel
  }
  
  if(is.null(optns$delta)){
    delta = max(c( diff(range(y))/1000, min(diff(sort(unique(y)))) ))
  } else {
    delta = optns$delta
  }
  
  if(is.null(optns$nRegGrid)){
    nRegGrid = 101
  } else {
    nRegGrid = optns$nRegGrid
  }
  
  if(is.null(optns$outputGrid)){
    outputGrid = NULL
  } else {
    outputGrid = optns$outputGrid
  }
  
  if(is.null(optns$infSupport)){
    infSupport = TRUE
  } else {
    infSupport = optns$infSupport
  }
  
  N = length(y)
  histgrid = seq( min(y)-delta*0.5, max(y)+delta*0.5, by = delta );
  if(  (max(y)+delta*0.5) > histgrid[length(histgrid)] ){
    histgrid[length(histgrid)] =  max(y)+delta*0.5;
  }
  M = length(histgrid)
  histObj =  hist(y, breaks = histgrid, plot = FALSE);
  yin = histObj$counts[1:M-1] / N / delta;
  xin = seq( min(y), max(y), by = delta);
  
  if( is.null(optns$userBwMu)){
    bw = fdapace:::CVLwls1D(y = yin, t = xin, kernel = kernel, npoly = 1, nder = 0, dataType = 'Dense', kFolds = 5)
  } else {
    bw = optns$userBwMu
  }
  
  densObj <- list()
  densObj$bw <- bw
  densObj$x <- outputGrid
  if( infSupport ){ 
    if(is.null(outputGrid)){
      densObj$x <- seq(min(y)-delta, max(y)+delta, length.out = nRegGrid)
    }
    # remove qpadding here
    qpadding = 10
    mu = fdapace::Lwls1D(bw = bw, kernel_type = kernel, win = rep(1,M+(-1+2*qpadding)), 
                         xin = c( min(y) - 11:2 *delta, xin, max(y) + 2:11 * delta), yin = c(rep(0,qpadding),yin,rep(0,qpadding)), xout = densObj$x)
  } else {
    if(is.null(outputGrid)){
      densObj$x <- seq(min(y), max(y), length.out = nRegGrid)
    }
    mu = fdapace::Lwls1D(bw = bw, kernel_type = kernel, win = rep(1,M-1), xin = xin, yin = yin, xout = densObj$x)
  }
  mu[mu<0] = 0;
  densObj$y = mu / fdapace:::trapzRcpp(densObj$x, mu);
  
  return(densObj)
  
}
