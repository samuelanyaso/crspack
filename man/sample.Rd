% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample.R
\docType{data}
\name{sample}
\alias{sample}
\title{Simulated cluster-correlated data}
\format{
A data frame with 3 variables:
\describe{
\item{cID}{Cluster ID}
\item{Y}{Cluster-correlated response variable}
\item{Z}{Unit-level continuous covariate}
...
}
}
\usage{
sample
}
\description{
A simple simulated data set with 3 clusters. The number of units within each cluster is
generated from Poisson with \eqn{\lambda=4}. The unit-level covariate, \emph{Z} is generated from Normal(\eqn{\mu=0, \sigma=0.25}).
The cluster-correlated responses \emph{Y} are obtained by \eqn{Y = 0.1*Z + \nu + \epsilon}, where \eqn{\nu} is a cluster-level
random effect generated from Normal(\eqn{\mu=0, \sigma=0.25}) and \eqn{\epsilon} is a random-error term generated
from Normal(\eqn{\mu=0, \sigma=0.25}).
}
\keyword{datasets}
