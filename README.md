# Rank-sum tests for cluster-correlated data

## Description 
This package contains functions to perform inference based on rank-sum statistics for cluster-correlated data with informativeness of the total cluster size, informativeness of a binary covariate distribution or informativeness of a subject-level covariate distribution.

## Installation
The functions require the R packages `Rcpp` and `RcppArmadillo` to be installed.

The `crspack` package can be installed using the `devtools` package:

    library(devtools) 
    devtools::install_github(repo="samuelanyaso/crspack") 

## Example
A simulated clustered-correlated data set `sample` is included in the package. The data set contains three variables, `cID`: cluster identification number, `Y`: the response and `Z`: the unit-level covariate. Based on the toy dataset, one can estimate the p-value for the test statistic proposed by Anyaso-Samuel and Datta as follows:

    library(crspack) 
    head(sample)
    ASDpv(mdat=sample, swapG = TRUE, iICGs=FALSE, useZ='percentiles',
      npctiles=10, Zvals=NULL, K = 10, parallel=TRUE, cores=5)




