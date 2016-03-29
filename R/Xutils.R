####============================================================
##  matchmz
##
##  Find the closest match based on mz value.
##  input: a data.frame or matrix with columns: mz, intensity.
##         if intensity is missing it just returns the closest entry,
##         otherwise, if more than one are within the range, the one
##         with the largest intensity.
##  mzdev: max allowed deviation.
##  ppm:
##  returns: a data.frame with values pos, mass, dif (position/idx,
##           mass and difference in masses)
##  discussion: distance to mz or to mzmin/mzmax? If mzmin and mzmax
##              are present I could use as distance the mean of the
##              distances to mzmin and mzmax.
####------------------------------------------------------------
.matchmz <- function(x, mz, mzdev=0.001, ppm=5){
    ## Support that x is only a vector of mz values.
    if(is.null(ncol(x))){
        x <- matrix(x, ncol=1)
        colnames(x) <- "mz"
    }
    if(!any(colnames(x) == "mz"))
        stop("'x' does not have the required column named 'mz'!")
    ## For now, drop columns with name mzmin and mzmax... eventually add them later.
    x <- x[, !(colnames(x) %in% c("mzmin", "mzmax")), drop=FALSE]
    ## ## Need an "average" mz value, if it's not there.
    ## if(all(c("mzmin", "mzmax") %in% colnames(x)) & !(any(colnames(x) == "mz")))
    ##     x <- cbind(x, mz=rowMeans(x[, c("mzmin", "mzmax")]))

    res <- lapply(mz, FUN=function(z){
        return(.singlematchmz(x=x, mz=z, mzdev=mzdev, ppm=ppm))
    })
    res <- do.call(rbind, res)
    return(res)
}
## The ranges thing doesn't work yet!!!
.singlematchmz <- function(x, mz, mzdev=0.001, ppm=5){
    if(all(c("mzmin", "mzmax") %in% colnames(x))){
        diffs <- mean(c(abs(x[, "mzmin"] - mz), abs(x[, "mzmax"] - mz)))
    }else{
        diffs <- abs(x[, "mz"] - mz)
    }
    ## Set those that are too far away to NA
    diffs[diffs > (mzdev + x[, "mz"]/(1000000 * ppm))] <- NA
    if(all(is.na(diffs))){
        ## Just return NA.
        return(c(pos=NA, dif=NA))
    }
    ## Get the index of those that are within the allowed deviation.
    validDiffs <- which(!is.na(diffs))
    ## Now it's just a matter of re-ordering them.
    if(any(colnames(x) == "intensity")){
        idx <- order(x[validDiffs, "intensity"], decreasing=TRUE)
    }else{
        idx <- order(diffs[validDiffs])
    }
    return(c(pos=validDiffs[idx[1]], dif=diffs[validDiffs][idx[1]]))
}


## ## Returns always a matrix with 3 columns: rt, mz, intensity. Is in principal the
## ## same thing as implemented in the rawMat method.
## .getWhatever <- function(xraw, mzrange=NULL, rtrange=NULL){
##     if(is.null(mzrange))
##         mzrange <- range(xraw@env$mz)
##     if(is.null(rtrange))
##         rtrange <- range(xraw@scantime)
##     ## Define which of the time points are within the time range.
##     timeRange <- range(which(xraw@scantime >= rtrange[1] & xraw@scantime <= rtrange[2]))

##     ## Get the indices of values within the mz and the rt range.
##     massIdx <- which(xraw@env$mz >= mzrange[1] & xraw@env$mz <= mzrange[2])
##     subsetIdx <- intersect()
##     ##scanRange <-
## }

