#' Calculate integral via trapezoid integration
#' 
#' Calculate the value of an integra via trapezoid integration
#' 
#' @param y A vector of values at time f(x) 
#' @param x A vector of time-points 'x'; has to be sorted. 
#' 
#' @return A scalar with the associated value
#' 
#' @examples
#' 
#' n <- 1001
#' x <- seq(0, 3*pi, len = n)
#' y <- sin(x)
#' fastTrapz(x, y) # 1.999985
#'   
#' @export 
#' 

fastTrapz <- function(x, y){
  
  return( fdapace:::trapzRcpp( X = x, Y = y) )
  
}
