##' Compute p-value for the ASD statistic
##'
##' Uses permutation tests to obtain the p-value of the ASD statistic. The procedure proceeds by permuting the subject-level covariate within each cluster, the ASD-statistic is computed for each permuted sample. The p-value is the proportion of the permuted statistic greater than the observed statistic.
##' @title Compute p-value for the ASD statistic
##' @param dw a numeric matrix with columns names 'cID', 'Y', and 'Z' that denotes the cluster ID, response, and subject-level covariate, respectively.
##' @param K an integer indicating the number of permutations. The default is 10000.
##' @return a list of three elements which are the observed test statistic, a numeric vector of the K permuted statistics, and the computed p-value.
##' @author Samuel Anyaso-Samuel, Somnath Datta
##' @import stats
##' @export
##' @useDynLib crspack
##' @references{
##'   \insertRef{dutta2016}{crspack}
##' }
##' @importFrom Rdpack reprompt
##' @importFrom Rcpp sourceCpp
##' @examples
##' ## simple data simulation
##' m <- 3
##' mdat <- lapply(1:m, function(x){
##'   ni <- rpois(n=1, lambda=10)+2 # cluster size
##'   cID <- rep(x=x, times=ni) # cluster ids
##'   nu <- rnorm(n=1, mean=0, sd=0.25) # cluster-level random effect
##'   Z <- rnorm(n=ni, mean=1, sd=0.25) # subject-level covariate
##'   Y <- 0.1*Z + nu + rnorm(n=ni, mean=0, sd=0.25) # response
##'   data.frame(cID=cID, Y=Y, Z=Z)
##' })
##' mdat <- do.call(rbind, mdat)
##' mdat <- apply(as.matrix(mdat),2,as.numeric)
##' 
##' ## Estimate the state occupation probabilities
##' ASDpv(dw=mdat, K=10)
ASDpv <- function(dw, K = 10000) {
    ASDpvC(dw, K)
}