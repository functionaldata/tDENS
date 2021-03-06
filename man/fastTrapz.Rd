% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastTrapz.R
\name{fastTrapz}
\alias{fastTrapz}
\title{Calculate integral via trapezoid integration}
\usage{
fastTrapz(x, y)
}
\arguments{
\item{x}{A vector of time-points 'x'; has to be sorted.}

\item{y}{A vector of values at time f(x)}
}
\value{
A scalar with the associated value
}
\description{
Calculate the value of an integra via trapezoid integration
}
\examples{

n <- 1001
x <- seq(0, 3*pi, len = n)
y <- sin(x)
fastTrapz(x, y) # 1.999985
  
}
