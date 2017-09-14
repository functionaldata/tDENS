#' Function to regularise densities to have (larger) minimum value
#' 
#' If possible, regularises the input density \code{y} to have minimum density value is \code{alpha}. See details.
#' 
#' @param x support of the density
#' @param y values of the density
#' @param alpha scalar to regularise with (default = 0.01) - this will be the minimum value of the regularised density, unless \code{min(y) > alpha}, in which case 
#' no regularisation will be performed
#' 
#' @details 
#' 
#' If \code{min(y) >= alpha} or \code{y} is the uniform distribution, no regularisation is performed and \code{y} is returned.  If \code{alpha*diff(range(x)) > 1}, 
#' the regularisation is not possible and an error is thrown.  Otherwise, the regularised density is computed by adding an appropriate constant \code{gam} \code{y}, 
#' followed by renormalisation to have integral 1.
#' 
#' @return dens density values on dSup
#' 
#' @seealso \code{\link{DeregulariseByAlpha}}
#' 
#' @examples
#' 
#'  x = seq(0,1,length.out=122)
#'  y = seq(0,2,length.out=122)
#'  z = RegulariseByAlpha(x=x, y=y, alpha = 0.1)
#'  
#' @export
 
RegulariseByAlpha <- function(x,y,alpha=0.01){
  
  if(abs(fdapace:::trapzRcpp(x, y) - 1) > 1e-5 || min(y) < 0){
    stop('y must be a density!')
  }
  
  if(min(y) >= alpha || all(diff(y) == 0)){
    y
  }else {
  
    totalBumpArea = diff(range(x))*alpha
    if(totalBumpArea >= 1){
      stop('alpha regularisation does not result in valid density.')
    }else {
      gam = (alpha - min(y))/(1 - totalBumpArea)
      return((y + gam)/(1 + gam*diff(range(x))))
    }
  }
}