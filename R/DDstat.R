##' Computes the Dutta & Datta statistic.
##'
##' Estimates the test statistic, expected value and variance of \insertCite{dutta2016;textual}{crspack}.
##' @title The Dutta & Datta statistic
##' @param mdat a numeric matrix with columns names 'cID', 'Y', and 'G' that denotes the cluster ID, response, and group indicator, respectively.
##' @param iICGs a logical argument indicating if the data has incomplete ICG structure. 
##' @param get_var a logical argument indicating whether the variance of the estimate should be obtained. Default is `TRUE`.
##' @return a list of the DD test statistic, expected value, variance (optional), and Z-score (optional).
##' @author Samuel Anyaso-Samuel, Somnath Datta
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
##'   ni <- rpois(n=1, lambda=4)+2 # cluster size
##'   cID <- rep(x=x, times=ni) # cluster ids
##'   nu <- rnorm(n=1, mean=0, sd=0.25) # cluster-level random effect
##'   G <- rbinom(n=ni, size=1, prob=0.5) # grouping variable
##'   err <- rnorm(n=ni, mean=0, sd=0.25)
##'   err[which(G == 1)] <- rnorm(n=sum(G == 1), mean=0.15, sd=0.25)
##'   Y <- 0.5 + nu + err # response
##'   data.frame(cID=cID, Y=Y, G=G)
##' })
##' mdat <- do.call(rbind, mdat)
##' mdat <- apply(as.matrix(mdat),2,as.numeric)
##' 
##' ## Estimate the DD statistic
##' DDstatV2(mdat=mdat, iICGs=TRUE, get_var=TRUE)
DDstatV2 <- function(mdat, iICGs = TRUE, get_var = TRUE) {

    # figure out the data
    arg2 <- c("cID", "Y", "G")
    arg2chk <- charmatch(arg2, colnames(mdat))
    if (any(is.na(arg2chk))) {
        stop("mdat should contain 'cID', 'Y', and 'G'")
    }
    mdat <- mdat[, arg2]

    # create the data matrix
    mdat <- apply(as.matrix(mdat), 2, as.numeric)

    ######################################### begin computation for the DD statistic
    res <- DDstatV2C(dw = mdat, iICGs = iICGs, get_var = get_var)

    return(res)
}
