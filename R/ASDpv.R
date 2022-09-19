##' Compute p-value for the ASD statistic.
##'
##' Uses permutation tests to obtain the p-value of the ASD statistic. The procedure proceeds by permuting the subject-level covariate within each cluster, the ASD-statistic is computed for each permuted sample. The p-value is the proportion of the permuted statistic greater than the observed statistic.
##' @title p-value for the ASD statistic
##' @param mdat a numeric matrix with columns names 'cID', 'Y', and 'Z' that denotes the cluster ID, response, and subject-level covariate, respectively.
##' @param iICGs a logical argument indicating if values of *Z* that creates incomplete ICG should be discarded. Default is `TRUE`.
##' @param useZ a character argument that indicates if 'all', 'percentiles', or 'prespecified' values of Z should be used to create data sets with a binary grouping variable. The \insertCite{dutta2016;textual}{crspack} statistic will be applied to these data sets. Default is 'percentiles'.
##' @param npctiles an integer indicating the number of equally spaced percentile points, if useZ='percentiles'. Default is 10.
##' @param Zvals a numeric vector indicating the pre-specified values of *Z* if useZ='prespecified'.
##' @param method a character argument indicating whether the DD statistic should be 'scaled' or not ('nonscaled') by its standard deviation. Default is 'nonscaled'.
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
##'   Z <- rnorm(n=ni, mean=1, sd=0.25) # subject-level covariate
##'   Y <- 0.1*Z + nu + rnorm(n=ni, mean=0, sd=0.25) # response
##'   data.frame(cID=cID, Y=Y, Z=Z)
##' })
##' mdat <- do.call(rbind, mdat)
##' mdat <- apply(as.matrix(mdat),2,as.numeric)
##' 
##' ## Estimate the p-value of the ASD statistic
##' ASDpv(mdat=mdat, iICGs=TRUE, useZ='percentiles',
##' npctiles=10, Zvals=NULL, method='nonscaled',
##' K = 10, parallel=TRUE, cores=5)
ASDpv <- function(mdat, iICGs = TRUE, useZ = "percentiles", npctiles = 10, Zvals = NULL, method = "nonscaled", K = 10000, parallel = TRUE,
    cores = 5) {

    # figure out useZ
    arg0 <- c("all", "percentiles", "prespecified")
    arg0chk <- charmatch(useZ, arg0)
    if (is.na(arg0chk)) {
        stop("useZ should be either 'all', 'percentiles' or 'prespecified' ")
    }

    # figure out the method
    arg1 <- c("nonscaled", "scaled")
    arg1chk <- charmatch(method, arg1)
    if (is.na(arg1chk)) {
        stop("method should be either 'nonscaled' or 'scaled'")
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

    if (parallel == TRUE) {

        ## compute the observed test statistic
        obstat <- ASDstatV3R(mdat = mdat, iICGs = iICGs, useZ = useZ, npctiles = npctiles, method = method)
        Zvals <- obstat$Zvals
        obstat <- obstat$vstat

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
        permstat <- foreach(k = 1:K, .combine = "c", .errorhandling = "remove", .export = c("ASDstatV3R"), .packages = c("doParallel")) %dopar%
            {

                ## create the permuted sample
                pdat <- mdat

                # permute the Z within each cluster
                pdat[, "Z"] <- unlist(lapply(unique(mdat[, "cID"]), function(x) {
                  sample(mdat[which(mdat[, "cID"] == x), "Z"], replace = FALSE)
                }))

                # compute the ASDstat for the permuted data
                if (useZ != "all") {
                  pstat <- ASDstatV3R(mdat = pdat, iICGs = iICGs, useZ = "prespecified", Zvals = Zvals, method = method)
                } else {
                  pstat <- ASDstatV3R(mdat = pdat, iICGs = iICGs, useZ = "all", method = method)
                }

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
        res <- ASDpvC(dw = mdat, K = K)
    }
    return(res)
}
