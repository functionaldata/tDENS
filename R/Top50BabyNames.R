#' Baby name popularity densities for 50 male and 50 female names in the USA
#' 
#' Baby name popularity densities, obtained by smoothing year-to-year popularity indices from 1950 to 2016, after normalization to have integral equal to 1.  
#' The top 50 names, in absolute popularity, are included for each gender.
#'
#' @name Top50BabyNames
#' @docType data
#' @format A list with two variables
#' \describe{
#' \item{x}{grid of years between 1950 and 2016, of length 67.}
#' \item{dens}{list of length two, corresponding to male (\code{dens$male}) and female(\code{dens$female}) names.  Each is a 67-by-50 matrix of density estimates, where each 
#' column corresponds to a unique baby name given by the corresponding column name.}
#' }  
#' @references
#' Data from the \code{R} package \code{babynames}, originally from the US Social Security Administration
NULL
