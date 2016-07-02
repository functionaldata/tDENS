#' Function to regularise non-strictly positive densities by a scalar alpha
#' 
#' Completely uniform distributions cannot be regularised by a level alpha.
#' 
#' @param x support of the density
#' @param y values of the density
#' @param alpha scalar to regularise with (default = 0.01)
#' @param deregularise logical to deregularise by alpha (instead of regularising) (default: FALSE )
#' 
#' @return dens density values on dSup
#' 
#' @examples
#' 
#'  x = seq(0,1,length.out =122)
#'  y = seq(0,2,length.out =122)
#'  z = RegulariseByAlpha(x=x, y=y, alpha = 0.1)
#'  
#' @export
 
RegulariseByAlpha <- function(x,y,alpha=0.01, deregularise = FALSE){ 
  totalBumpArea = diff(range(x)) * alpha;
  if( (1 - totalBumpArea) <= 0){
    stop('alpha regularisation does not result in valid density.')
  }
  ifelse(!deregularise, return( (y * (1 - totalBumpArea)) + alpha ), return( (y - alpha) / (1 - totalBumpArea) ))
}

# For convenience
DeregulariseByAlpha <- function(x,y,alpha){ 
  totalBumpArea = diff(range(x)) * alpha;
  if( (1 - totalBumpArea) <= 0){
    stop('alpha regularisation does not result in valid density.')
  }
  return( (y - alpha) / (1 - totalBumpArea) )
}
