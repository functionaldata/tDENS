% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateModeOfVarPlotLQ2D.R
\name{CreateModeOfVarPlotLQ2D}
\alias{CreateModeOfVarPlotLQ2D}
\title{Functional Principal Component Analysis mode of variation plot}
\usage{
CreateModeOfVarPlotLQ2D(fpcaObj, dSup = NULL, k = 1, domain = "Q",
  alpha = 0, numOfModes = 7, dSupPlot = NULL, ...)
}
\arguments{
\item{fpcaObj}{An FPCA class object returned by FPCA().}

\item{dSup}{The support of the original density used by LQD(relevant only for density domain)}

\item{k}{The k-th mode of variation to plot (default k = 1)}

\item{domain}{character defining if we should plot on the Quantile ('Q') or the Density ('D') domain (default: 'Q')}

\item{alpha}{regularisation parameter}

\item{numOfModes}{scalar number of principal modes to plot (relevant only for density domain, needs to be an odd number >1.) (default: 7)}

\item{dSupPlot}{The support of the original density used for plotting (relevant only for density domain)}

\item{...}{Additional arguments for the 'plot' function.}
}
\description{
Create the k-th mode of variation plot around the mean. The red-line is
the functional mean, the grey shaded areas show the range of variations
around the mean: \eqn{ \pm Q \sqrt{\lambda_k} \phi_k}{+/- Q sqrt{lambda_k} phi_k}
for the dark grey area Q = 1, and for the light grey are Q = 2.
}
\examples{
library(fdapace)
set.seed(1)
n <- 20
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- Sparsify(sampWiener, pts, 10)
res <- FPCA(sampWiener$Ly, sampWiener$Lt)
CreateModeOfVarPlotLQ2D(res)
}
