##' Obtains the ASD statistic.
##'
##' Computes the ASD statistic by aggregating over the \insertCite{dutta2016;textual}{crspack} statistics (DD statistic, herein) obtained from the data sets with binary groups formed by selected values of the subject-level continuous covariate.
##' @title The ASD statistic
##' @param mdat a numeric matrix with columns names 'cID', 'Y', and 'Z' that denotes the cluster ID, response, and subject-level covariate, respectively.
##' @param swapG a logical argument set to TRUE if we should swap the group indicator of a random observation in a cluster with incomplete ICG size. Default is `TRUE`.
##' @param iICGs a logical argument indicating if values of Z that creates incomplete ICG should be discarded. If `swapG=TRUE` then set to `FALSE`. Default is `FALSE`.
##' @param useZ a character argument that indicates if 'all', 'percentiles', or 'prespecified' values of *Z* should be used to create data sets with a binary grouping variable. The DD statistic will be applied to these data sets. Default is 'percentiles'.
##' @param npctiles an integer indicating the number of equally spaced percentile points, if useZ='percentiles'. Default is 10.
##' @param Zvals a numeric vector indicating the pre-specified values of *Z* if useZ='prespecified'.
##' @return a list of three elements which are the ASD test statistic, a numeric vector of DD statistic evaluated at the Zs, and a numeric vector indicating the values/percentiles of *Z* at which the DD stat were evaluated (returns an empty vector if useZ='all').
##' @author Samuel Anyaso-Samuel, Somnath Datta
##' @import doParallel stats
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
##'   Z <- rnorm(n=ni, mean=1, sd=0.25) # subject-level covariate
##'   Y <- 0.1*Z + nu + rnorm(n=ni, mean=0, sd=0.25) # response
##'   data.frame(cID=cID, Y=Y, Z=Z)
##' })
##' mdat <- do.call(rbind, mdat)
##' mdat <- apply(as.matrix(mdat),2,as.numeric)
##'
##' ## Estimate the ASD statistic
##' ASDstatV3R(mdat=mdat, swapG = TRUE, iICGs=FALSE, useZ='percentiles', npctiles=10, Zvals=NULL)
ASDstatV3R <- function(mdat, swapG = TRUE, iICGs = FALSE, useZ = "percentiles", npctiles = 10, Zvals = NULL) {

    # figure out useZ
    arg0 <- c("all", "percentiles", "prespecified")
    arg0chk <- charmatch(useZ, arg0)
    if (is.na(arg0chk)) {
        stop("useZ should be either 'all', 'percentiles' or 'prespecified' ")
    }

    # figure out the data
    arg2 <- c("cID", "Y", "Z")
    arg2chk <- charmatch(arg2, colnames(mdat))
    if (any(is.na(arg2chk))) {
        stop("mdat should contain 'cID', 'Y', and 'Z'")
    }
    mdat <- mdat[, arg2]

    # create the data matrix
    mdat <- apply(as.matrix(mdat), 2, as.numeric)
    mdat <- mdat[order(mdat[, "cID"]), ]  # sort data by cID

    ######################################### begin computation for the ASD statistic
    Z <- mdat[, "Z"]

    # values of Z at which DD stat should be evaluated
    if (useZ == "all") {

        # use all values of Z - if this option is used, parallelize this function and unparallelize the ASDpvR function.
        Zvals <- Z
    } else if (useZ == "percentiles") {

        # use the values of Z at certain percentiles
        Zvals <- as.numeric(quantile(Z, sort(c(0.5, seq(0, 1, length.out = npctiles)))))
    } else if (useZ == "prespecified") {

        # use pre-specified values of Z
        Zvals <- Zvals
    }

    ## computes the ASD statistic
    res0 <- ASDstatV3C(dw = mdat, Z = Zvals, iICGs = iICGs, swapG = swapG)

    if (useZ == "all")
        Zvals = NULL

    res <- list(vstat = res0$v_stat, Vstar = res0$Vstar, Zvals = Zvals)
    return(res)
}
