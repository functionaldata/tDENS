#' Function to deregularise densities to have (smaller) minimum value
#' 
#' If possible, deregularises the input density \code{y} to have minimum density value is \code{alpha}. See details.
#' 
#' @param x support of the density
#' @param y values of the density
#' @param alpha scalar to deregularise with (default = 0) - this will be the minimum value of the deregularised density, unless \code{min(y) < alpha}, in which case 
#' no deregularisation will be performed
#' 
#' @details 
#' 
#' If \code{min(y) <= alpha}, or \code{y} is the uniform distribution, no deregularisation is performed and \code{y} is returned.  If \code{min(y)*diff(range(x)) > 1}, 
#' the deregularisation is not possible and an error is thrown.  Otherwise, the deregularised density in an inverse manner to \code{RegulariseByAlpha}.
#' 
#' @return dens density values on x
#' 
#' @seealso \code{\link{RegulariseByAlpha}}
#' 
#' @examples
#' 
#'  x = seq(0,1,length.out=122)
#'  y = seq(0.1,1.9,length.out=122)
#'  z = DeregulariseByAlpha(x=x, y=y, alpha = 0)
#'  
#' @export
 
DeregulariseByAlpha <- function(x,y,alpha=0){

  if(abs(fdapace:::trapzRcpp(x, y) - 1) > 1e-5 || min(y) < 0){
    stop('y must be a density!')
  }
  
  if(min(y) <= alpha || all(diff(y) == 0)){
    y
  }else {
  
    totalBumpArea = diff(range(x))*min(y)
    if(totalBumpArea >= 1){
      stop('alpha deregularisation does not result in valid density.')
    }else {
      gam = (min(y) - alpha)/(1 - totalBumpArea)
      return((1 + gam*diff(range(x)))*y - gam)
    }
  }
}