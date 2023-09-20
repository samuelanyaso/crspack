##' Compute p-value for the test statistic proposed by Anyaso-Samuel and Datta (ASD). The p-value for testing the effect of a continuous unit-level covariate on a cluster-correlated response while accounting for informativeness.
##'
##' Uses permutation tests to obtain the p-value of the ASD statistic. The procedure proceeds by permuting the unit-level covariate within each cluster, the ASD-statistic is computed for each permuted sample. The p-value is the proportion of the permuted statistic greater than the observed statistic.
##' @title p-value for the ASD statistic
##' @param mdat a numeric matrix with columns names 'cID', 'Y', and 'Z' that denotes the cluster ID, response, and unit-level covariate, respectively.
##' @param swapG a logical argument set to TRUE if we should swap the group indicator of a random observation in a cluster with incomplete ICG size. Default is `TRUE`.
##' @param iICGs a logical argument indicating if values of Z that creates incomplete ICG should be discarded. If `swapG=TRUE` then set to `FALSE`. Default is `FALSE`.
##' @param useZ a character argument that indicates if 'all', 'percentiles', or 'prespecified' values of *Z* should be used to create data sets with a binary grouping variable. The DD statistic will be applied to these data sets. Default is 'percentiles'.
##' @param npctiles an integer indicating the number of equally spaced percentile points, if useZ='percentiles'. Default is 10.
##' @param Zvals a numeric vector indicating the pre-specified values of *Z* if useZ='prespecified'.
##' @param K an integer indicating the number of permutations. The default is 10000.
##' @param parallel a logical argument indicating whether to parallelize the permuatation procedure.
##' @param cores an integer indicating the number of cores to use for parallelization. For instance see `parallel::detectCores()`.
##' @return a list of three elements which are the observed test statistic, a numeric vector of the *K* permuted statistics, and the computed p-value.
##' @author Samuel Anyaso-Samuel, Somnath Datta
##' @import doParallel parallel foreach stats
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
##'   Z <- rnorm(n=ni, mean=1, sd=0.25) # unit-level covariate
##'   Y <- 0.1*Z + nu + rnorm(n=ni, mean=0, sd=0.25) # response
##'   data.frame(cID=cID, Y=Y, Z=Z)
##' })
##' mdat <- do.call(rbind, mdat)
##' mdat <- apply(as.matrix(mdat),2,as.numeric)
##'
##' ## Estimate the p-value of the ASD statistic
##' ASDpv(mdat=mdat, swapG = TRUE, iICGs=FALSE, useZ='percentiles',
##' npctiles=10, Zvals=NULL, K = 10, parallel=TRUE, cores=5)
ASDpv <- function(mdat, swapG = TRUE, iICGs = FALSE, useZ = "percentiles", npctiles = 10, Zvals = NULL, K = 10000, parallel = TRUE, cores = 5) {

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

    ## compute the observed test statistic
    Z <- mdat[, "Z"]

    # values of Z at which DD stat should be evaluated
    if (useZ == "all") {

        # use all values of Z - if this option is used, parallelize this function and unparallelize the ASDpvR function.
        Zvals <- Z
    } else if (useZ == "percentiles") {

        # use the values of Z at certain percentiles
        Zvals <- as.numeric(quantile(Z, seq(0, 1, length.out = npctiles)))
    } else if (useZ == "prespecified") {

        # use pre-specified values of Z
        Zvals <- Zvals
    }
    obstat <- ASDstat(mdat = mdat, swapG = swapG, iICGs = iICGs, useZ = "prespecified", Zvals = Zvals)
    obstat <- obstat$vstat

    if (parallel == TRUE) {

        ## set up parallel backend

        ## for performance reasons, CRAN limits the num of cores available to packages to 2
        chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
        if (nzchar(chk) && chk == "TRUE") {
            # use 2 cores in CRAN/Travis/AppVeyor
            cores <- 2L
        }

        # initiate cluster
        cores <- as.integer(cores)
        cl <- parallel::makeCluster(cores)
        registerDoParallel(cl)

        ## compute the ASDstat based on permuted samples
        permstat <- foreach(k = 1:K, .combine = "c", .errorhandling = "remove", .export = c("ASDstat"), .packages = c("doParallel")) %dopar%
            {

                ## create the permuted sample
                pdat <- mdat

                # permute the Z within each cluster
                pdat[, "Z"] <- unlist(lapply(unique(mdat[, "cID"]), function(x) {
                  sample(mdat[which(mdat[, "cID"] == x), "Z"], replace = FALSE)
                }))

                # compute the ASDstat for the permuted data
                pstat <- ASDstat(mdat = pdat, swapG = swapG, iICGs = iICGs, useZ = "prespecified", Zvals = Zvals)

                # return the values for both stats
                pstat$vstat
            }

        # stop cluster
        parallel::stopCluster(cl)

        # compute the p-value for the ASD statistic
        permstat <- permstat[!is.na(permstat)]
        pvalue <- sum(permstat >= obstat)/length(permstat)

        res <- list(obstat = obstat, permstat = permstat, pvalue = pvalue)
    } else {
        ## Note that the ASD statistic embedded in this function is based on the data sets created by all Z.
        res <- ASDpv3C(dw = mdat, K = K, Zvals = Zvals, swapG = swapG, iICGs = iICGs)
    }
    return(res)
}
