####============================================================
##  Methods for MSslice
##
####------------------------------------------------------------

####============================================================
##  show
##
####------------------------------------------------------------
setMethod("show", "MSslice", function(object){
    cat(class(object), " object:\n", sep="")
    cat("| Number of MSdata objects: ", length(object@data), "\n", sep="")
    if(length(object@mzrange) == 2)
        cat("| m/z range: ", object@mzrange[1], " - ", object@mzrange[2], "\n", sep="")
    if(length(object@rtrange) == 2)
        cat("| RT range: ", object@rtrange[1], " - ", object@rtrange[2], "\n", sep="")
})

####============================================================
##  rtrange
##
##  Getter the rtrange slot.
####------------------------------------------------------------
setMethod("rtrange", "MSslice", function(object){
    return(object@rtrange[1:2])
})

####============================================================
##  mzrange
##
##  Getter the mzrange slot.
####------------------------------------------------------------
setMethod("mzrange", "MSslice", function(object){
    return(object@mzrange[1:2])
})

####============================================================
##  intrange
##
##  Get the intensity range.
####------------------------------------------------------------
setMethod("intrange", "MSslice", function(object){
    ints <- range(unlist(lapply(msData(object), intrange)))
    return(ints)
})

####============================================================
##  msData
##
##  Get the MSdata object.
####------------------------------------------------------------
setMethod("msData", "MSslice", function(object, ...){
    return(object@data)
})

####============================================================
##  length
##
##  Get the number of MSdata objects
####------------------------------------------------------------
setMethod("length", "MSslice", function(x){
    return(length(x@data))
})

####============================================================
##  getChromatogram
##
##  Extracts the chromatogram and returns a matrix with the values.
####------------------------------------------------------------
setMethod("getChromatogram", "MSslice",
          function(object, FUN=max, bins=NULL, nbin=NULL, binSize=NULL, ...){
              ## Extract a list with a chromatogram-matrix for each sample.
              chrL <- .getChromList(object, FUN=FUN, bins=bins, nbin=nbin,
                                    binSize=binSize)
              return(.list2mat(chrL))
          })
## That one extracts the chromatogram and returns a list with matrices,
## one for each MSdata.
.getChromList <- function(x, FUN=max, bins=NULL, nbin=NULL, binSize=NULL){
    if(is.null(bins)){
        if(!is.null(nbin) | !is.null(binSize)){
            bins <- .getBins(rtrange(x), nbin=nbin, binSize=binSize)
            nbin <- NULL
            binSize <- NULL
        }
    }
    resList <- lapply(msData(x), FUN=function(z, theFun){
        return(getChromatogram(z, FUN=theFun, bins=bins, nbin=nbin,
                               binSize=binSize))
    }, theFun=FUN)
    return(resList)
}
## This function converts a list into a matrix.
.list2mat <- function(x){
    ## Get the list of unique time points
    unt <- sort(unique(unlist(lapply(x, function(z)z[,1]))))
    vals <- lapply(x, function(z){
        tmp <- rep(NA, length(unt))
        tmp[match(z[, 1], unt)] <- z[, 2]
        return(tmp)
    })
    Res <- do.call(cbind, vals)
    rownames(Res) <- unt
    return(Res)
}

####============================================================
##  plotChromatogram
##
##  Basically, plotting a chromatogram for each of the samples.
####------------------------------------------------------------
setMethod("plotChromatogram", "MSslice",
          function(object, FUN=max, bins=NULL, nbin=NULL, binSize=NULL,
                   add=FALSE, main=paste(format(mzrange(object), 2), collapse="-"),
                   xlab="Retention time", ylab="Intensity",
                   ...){
              ## First get all of the chromatograms.
              chrs <- lapply(msData(object), function(z, theFun, ...){
                  return(getChromatogram(z, FUN=theFun, ...))
              }, theFun=FUN, )
              ## Get the range of intensity values, range of rts.
              xlim <- rtrange(object)
              ## Plot the empty plot.
              ## Check if we've got col or lty specified.
              theDots <- list(...)
              ## plot the individual chromatograms; for loop.

              ##
              dat <- getChromatogram(object, FUN=FUN, bins=bins, nbin=nbin,
                                     binSize=binSize)
              if(add){
                  points(dat[, 1], dat[, 2], ...)
              }else{
                  plot(dat[, 1], dat[, 2], main=main, xlab=xlab, ylab=ylab, ...)
              }
          })
