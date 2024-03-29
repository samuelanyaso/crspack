% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ASDstat.R
\name{ASDstat}
\alias{ASDstat}
\title{The ASD statistic}
\usage{
ASDstat(
  mdat,
  swapG = TRUE,
  iICGs = FALSE,
  useZ = "percentiles",
  npctiles = 10,
  Zvals = NULL
)
}
\arguments{
\item{mdat}{a numeric matrix with columns names 'cID', 'Y', and 'Z' that denotes the cluster ID, response, and unit-level covariate, respectively.}

\item{swapG}{a logical argument set to TRUE if we should swap the group indicator of a random observation in a cluster with incomplete ICG size. Default is \code{TRUE}.}

\item{iICGs}{a logical argument indicating if values of Z that creates incomplete ICG should be discarded. If \code{swapG=TRUE} then set to \code{FALSE}. Default is \code{FALSE}.}

\item{useZ}{a character argument that indicates if 'all', 'percentiles', or 'prespecified' values of \emph{Z} should be used to create data sets with a binary grouping variable. The DD statistic will be applied to these data sets. Default is 'percentiles'.}

\item{npctiles}{an integer indicating the number of equally spaced percentile points, if useZ='percentiles'. Default is 10.}

\item{Zvals}{a numeric vector indicating the pre-specified values of \emph{Z} if useZ='prespecified'.}
}
\value{
a list of three elements which are the ASD test statistic, a numeric vector of DD statistic evaluated at the Zs, and a numeric vector indicating the values/percentiles of \emph{Z} at which the DD stat were evaluated (returns an empty vector if useZ='all').
}
\description{
Obtains the test statistic proposed by Anyaso-Samuel and Datta (ASD). The statistic estimates the effect of a continuous unit-level covariate on a cluster-correlated response while accounting for informativeness.
}
\details{
Computes the ASD statistic by aggregating over the \insertCite{dutta2016;textual}{crspack} statistics (DD statistic, herein) obtained from the data sets with binary groups formed by selected values of the unit-level continuous covariate.
}
\examples{
## simple data simulation
m <- 3
mdat <- lapply(1:m, function(x){
  ni <- rpois(n=1, lambda=4)+2 # cluster size
  cID <- rep(x=x, times=ni) # cluster ids
  nu <- rnorm(n=1, mean=0, sd=0.25) # cluster-level random effect
  Z <- rnorm(n=ni, mean=1, sd=0.25) # unit-level covariate
  Y <- 0.1*Z + nu + rnorm(n=ni, mean=0, sd=0.25) # response
  data.frame(cID=cID, Y=Y, Z=Z)
})
mdat <- do.call(rbind, mdat)
mdat <- apply(as.matrix(mdat),2,as.numeric)

## Estimate the ASD statistic
ASDstat(mdat=mdat, swapG = TRUE, iICGs=FALSE, useZ='percentiles', npctiles=10, Zvals=NULL)
}
\references{
{
\insertRef{dutta2016}{crspack}
}
}
\author{
Samuel Anyaso-Samuel, Somnath Datta
}
