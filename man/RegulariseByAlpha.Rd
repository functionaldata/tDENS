% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RegulariseByAlpha.R
\name{RegulariseByAlpha}
\alias{RegulariseByAlpha}
\title{Function to regularise non-strictly positive densities by a scalar alpha}
\usage{
RegulariseByAlpha(x, y, alpha = 0.01, deregularise = FALSE)
}
\arguments{
\item{x}{support of the density}

\item{y}{values of the density}

\item{alpha}{scalar to regularise with (default = 0.01)}

\item{deregularise}{logical to deregularise by alpha (instead of regularising) (default: FALSE )}
}
\value{
dens density values on dSup
}
\description{
Completely uniform distributions cannot be regularised by a level alpha.
}
\examples{

 x = seq(0,1,length.out =122)
 y = seq(0,2,length.out =122)
 z = RegulariseByAlpha(x=x, y=y, alpha = 0.1)
 
}
