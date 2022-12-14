% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00_FEM.R
\name{sampling_from_BT}
\alias{sampling_from_BT}
\title{Sample from the approximate prior using BT splines}
\usage{
sampling_from_BT(x, a, k, region, boundary = TRUE, n = 1)
}
\arguments{
\item{x}{A vector that specifies where to evaluate the sample path.}

\item{a}{A positive scalar represents the periodicity parameter.}

\item{k}{A positive integer represents the number of knots used to define the BT basis. The number of
basis functions equals to 2 times k or 2 times (k-1) if boundary is TRUE.}

\item{region}{A vector of size 2 that denotes the upper and lower interval limit of the region of interest.}

\item{boundary}{A logical value to denote whether to consider the boundary conditions.}

\item{n}{The number of samples to draw.}
}
\value{
A matrix with each column denote a sample path f(x).
}
\description{
This function samples sample paths from the approximate prior using BT splines.
}
\examples{
sampling_from_BT(x = c(1,2,3), a = 1, k = 5, region = c(0,5))
}
