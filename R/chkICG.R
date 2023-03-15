##' Incomplete ICG structure
##'
##' Checks if data has incomplete ICG structure.
##' @title Incomplete ICG structure
##' @param mdat a numeric matrix with columns names 'cID', 'Y', and 'G' that denotes the cluster ID, response, and group indicator, respectively.
##' @param swapG a logical argument set to TRUE if we should swap the group indicator of a random observation in a cluster with incomplete ICG size. Default is `TRUE`.
##' @return a list of three elements which are - boolean indicating if data has incomplete ICG structure, a matrix of cluster id and indicators of whether such cluster has incomplete ICG, and the resulting dataframe.
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
##' ## Checks if data has incomplete ICG structure
##' chkICG(mdat=mdat, swapG=TRUE)
chkICG <- function(mdat, swapG = TRUE) {

    # figure out the data
    arg2 <- c("cID", "Y", "G")
    arg2chk <- charmatch(arg2, colnames(mdat))
    if (any(is.na(arg2chk))) {
        stop("mdat should contain 'cID', 'Y', and 'G'")
    }
    mdat <- mdat[, arg2]

    # create the data matrix
    mdat <- apply(as.matrix(mdat), 2, as.numeric)

    ######################################### Checks if data has incomplete ICG structure
    res <- foo7v2C(dw = mdat, swapG)

    return(res)
}
