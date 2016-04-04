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
                   main=paste(format(mzrange(object), 2), collapse="-"),
                   xlab="Retention time", ylab="Intensity", col=1, lty=1,
                   ...){
              ## col and lty check
              if(length(col) > 1){
                  if(length(col) != length(object)){
                      warning("Length of 'col' does not match length of 'object';",
                              " using only the first value.")
                      col <- rep(col[1], length(object))
                  }
              }else{
                  col <- rep(col, length(object))
              }
              if(length(lty) > 1){
                  if(length(lty) != length(object)){
                      warning("Length of 'lty' does not match length of 'object';",
                              " using only the first value.")
                      lty <- rep(lty[1], length(object))
                  }
              }else{
                  lty <- rep(lty, length(object))
              }
              chrM <- getChromatogram(object, FUN=FUN, bins=bins, nbin=NULL, binSize=NULL)
              ## Do some rounding here???
              xVals <- round(as.numeric(rownames(chrM)), digits=2)
              xlim <- range(xVals)
              ylim <- range(chrM, na.rm=TRUE)
              ## Plot the empty plot.
              plot(3, 3, pch=NA, main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
              ## plot the individual chromatograms; for loop.
              for(i in 1:ncol(chrM)){
                  points(xVals, chrM[, i], col=col[i], lty=lty[i], ...)
              }
          })

####============================================================
##  getSpectrum
##
##  The same as getChromatogram, but for the spectrum.
####------------------------------------------------------------
setMethod("getSpectrum", "MSslice",
          function(object, FUN=max, bins=NULL, nbin=NULL, binSize=NULL, ...){
              ## Extract a list with a chromatogram-matrix for each sample.
              spcL <- .getSpecList(object, FUN=FUN, bins=bins, nbin=nbin,
                                   binSize=binSize)
              return(.list2mat(spcL))
          })
## That one extracts the chromatogram and returns a list with matrices,
## one for each MSdata.
.getSpecList <- function(x, FUN=max, bins=NULL, nbin=NULL, binSize=NULL){
    if(is.null(bins)){
        if(!is.null(nbin) | !is.null(binSize)){
            bins <- .getBins(mzrange(x), nbin=nbin, binSize=binSize)
            nbin <- NULL
            binSize <- NULL
        }
    }
    resList <- lapply(msData(x), FUN=function(z, theFun){
        return(getSpectrum(z, FUN=theFun, bins=bins, nbin=nbin,
                           binSize=binSize))
    }, theFun=FUN)
    return(resList)
}

####============================================================
##  plotSpectrum
##
##  Basically, plotting a spectrum for each of the samples.
####------------------------------------------------------------
setMethod("plotSpectrum", "MSslice",
          function(object, FUN=max, bins=NULL, nbin=NULL, binSize=NULL,
                   main=paste(format(rtrange(object), 2), collapse="-"),
                   xlab="M/Z", ylab="Intensity", col=1, lty=1,
                   ...){
              ## col and lty check
              if(length(col) > 1){
                  if(length(col) != length(object)){
                      warning("Length of 'col' does not match length of 'object';",
                              " using only the first value.")
                      col <- rep(col[1], length(object))
                  }
              }else{
                  col <- rep(col, length(object))
              }
              if(length(lty) > 1){
                  if(length(lty) != length(object)){
                      warning("Length of 'lty' does not match length of 'object';",
                              " using only the first value.")
                      lty <- rep(lty[1], length(object))
                  }
              }else{
                  lty <- rep(lty, length(object))
              }
              spcM <- getSpectrum(object, FUN=FUN, bins=bins, nbin=NULL, binSize=NULL)
              xVals <- as.numeric(rownames(spcM))
              xlim <- range(xVals)
              ylim <- range(spcM, na.rm=TRUE)
              ## Plot the empty plot.
              plot(3, 3, pch=NA, main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
              ## plot the individual chromatograms; for loop.
              for(i in 1:ncol(spcM)){
                  points(xVals, spcM[, i], col=col[i], lty=lty[i], ...)
              }
          })
