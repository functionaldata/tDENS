% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MakeLQDsample.R
\name{MakeLQDsample}
\alias{MakeLQDsample}
\title{Convenience function for converting densities to log-quantile functions}
\usage{
MakeLQDsample(dmatrix, alpha = 0, dInput)
}
\arguments{
\item{dmatrix}{Matrix holding the density values on dInput - assumed to be strictly positive and integrate to 1}

\item{alpha}{Scalar to regularise the supports with (default=0)}

\item{dInput}{Support (grid) for Density domain}
}
\value{
list with the 'LQD' projected data, the 'alpha' used and the 'newSupport' that was defined after trimming the data.
}
\description{
See 'dens2lqd' and 'RegulariseByAlpha' for more details.
This function will automatically trim out support points that the 99.5-th percentile (calculated using all available samples) is not positive.
It will then normalise the new densities to integraate to one, regularise by 'alpha' and calculate the LQ projection.
}
\references{
\cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2015}
}
