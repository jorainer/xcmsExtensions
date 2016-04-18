####============================================================
##  Methods for MSdata
##
####------------------------------------------------------------

####============================================================
##  show
##
####------------------------------------------------------------
setMethod("show", "MSdata", function(object){
    cat(class(object), " object:\n", sep="")
    cat("| MS level: ", object@mslevel, "\n", sep="")
    if(length(object@mzrange) == 2)
        cat("| m/z range: ", object@mzrange[1], " - ", object@mzrange[2], "\n", sep="")
    if(length(object@rtrange) == 2)
        cat("| RT range: ", object@rtrange[1], " - ", object@rtrange[2], "\n", sep="")
    if(length(object@intrange) == 2)
        cat("| Intensity range: ", object@intrange[1], " - ", object@intrange[2], "\n", sep="")
    cat("| Number of data points: ", length(object@mz), "\n", sep="")
})

####============================================================
##  rtime
##
##  Getter/setter for the rtime slot.
####------------------------------------------------------------
setMethod("rtime", "MSdata", function(object, return.type="numeric"){
    return.type <- match.arg(return.type, c("numeric", "Rle"))
    if(return.type == "Rle")
        return(object@rtime)
    return(as.numeric(object@rtime))
})
setReplaceMethod("rtime", "MSdata", function(object, value){
    object@rtrange <- range(value)
    object@rtime <- Rle(value)
    validObject(object)
    return(object)
})

####============================================================
##  intensity
##
##  Getter for the intensity slot.
####------------------------------------------------------------
setMethod("intensity", "MSdata", function(object){
    return(object@intensity)
})
setReplaceMethod("intensity", "MSdata", function(object, value){
    object@intensity <- value
    object@intrange <- range(value)
    validObject(object)
    return(object)
})

####============================================================
##  mz
##
##  Getter for the mz slot.
####------------------------------------------------------------
setMethod("mz", "MSdata", function(object){
    return(object@mz)
})
setReplaceMethod("mz", "MSdata", function(object, value){
    object@mz <- value
    object@mzrange <- range(value)
    validObject(object)
    return(object)
})

####============================================================
##  rtrange
##
##  Getter the rtrange slot.
####------------------------------------------------------------
setMethod("rtrange", "MSdata", function(object){
    return(object@rtrange[1:2])
})

####============================================================
##  mzrange
##
##  Getter the rtrange slot.
####------------------------------------------------------------
setMethod("mzrange", "MSdata", function(object){
    return(object@mzrange[1:2])
})

####============================================================
##  intrange
##
##  Getter the rtrange slot.
####------------------------------------------------------------
setMethod("intrange", "MSdata", function(object){
    if(length(object@intensity) == 0){
        return(c(NA, NA))
    }
    return(range(intensity(object)))
})


####============================================================
##  msSlice
##
##  Create an MSslice object from the MSdata.
####------------------------------------------------------------
setMethod("msSlice", "MSdata",
          function(object, call=match.call(), ...){
              ## if(!is.null(mzrange) | !is.null(rtrange)){
              ##     warning("Ignoring arguments 'mzrange' and 'rtrange' if ")
              ## }
              res <- new("MSslice", data=list(object), mzrange=mzrange(object),
                         rtrange=rtrange(object), call=call)
              validObject(res)
              return(res)
          })

####============================================================
##  as.matrix
##
##  Transform the data into a matrix.
####------------------------------------------------------------
setMethod("as.matrix", "MSdata",
          function(x, ...){
              return(cbind(rtime=as.numeric(rtime(x)),
                           mz=mz(x),
                           intensity=as.numeric(intensity(x))))
          })


####============================================================
##  .aggregateWithBins
##
##  Aggregate values in x using method FUN depending on the bins.
##  Argument bins might be returned by the .getBins function.
.aggregateWithBins <- function(x, bins, FUN=sum){
    return(unlist(lapply(split(x, findInterval(x, bins, all.inside=TRUE)), FUN=FUN),
                  recursive=FALSE, use.names=FALSE))
}
.aggregateWithIntervals <- function(x, ints, FUN=sum){
    return(unlist(lapply(split(x, ints), FUN=FUN),
                  recursive=FALSE, use.names=FALSE))
}

####============================================================
##  .getBins
##
##  Define the interval/breaks to be used for the binning of x.
##  x can also be a range.
##  Arguments nbin and binSize are mutually exclusive.
.getBins <- function(x, nbin, binSize){
    if(missing(nbin))
        nbin <- NULL
    if(missing(binSize))
        binSize <- NULL
    if(!is.null(nbin) & !is.null(binSize))
        stop("Arguments 'nbin' and 'binSize' are mutually exclusive!")
    rangs <- range(x, na.rm=TRUE)
    if(!is.null(nbin)){
        return(seq(rangs[1], rangs[2], length.out=(nbin+1)))
    }else{
        tmp <- seq(rangs[1], rangs[2], by=binSize)
        return(c(tmp, tmp[length(tmp)] + binSize))
    }
}
####============================================================
##  .binMid
##
##  Get (roughly) the midpoint of equal spaced bins.
.binMid <- function(bins){
    tmp <- bins + diff(bins[1:2])/2
    return(tmp[-length(tmp)])
}

####============================================================
##  .getChrom
##
##  get the chromatogram. If bins is not specified, the function /just/
##  aggregates values with the same retention time using the function FUN.
##  x... MSdata object.
##  FUN... the function to aggregate intensity values.
##  bins... binning variable, supposed to be generated by the .getBins method.
##  nbin, binSize... if one of the two is not NULL, the bins will be internally
##                   calculated using the .getBins function.
##  The function returns a two-column matrix with the retention time and the
##  intensity (columns rtime and intensity).
####------------------------------------------------------------
.getChrom <- function(x, FUN=sum, bins=NULL, nbin=NULL, binSize=NULL){
    if(!is.null(nbin) & !is.null(binSize))
        stop("Arguments 'nbin' and 'binSize' are mutually exclusive")
    if(is.null(bins) & (!is.null(nbin) | !is.null(binSize))){
        ## Calculate the bins internally.
        bins <- .getBins(rtrange(x), nbin=nbin, binSize=binSize)
    }

    ## Do the aggregation on bins if it's not empty
    if(!is.null(bins)){
        ## Get the rt intervals for which we do have data.
        rtIntervals <- findInterval(rtimeOrdered(x), bins, all.inside=TRUE)
        aggInts <- .aggregateWithIntervals(intensityOrderedByRtime(x), rtIntervals,
                                           FUN=FUN)
        return(cbind(rtime=.binMid(bins)[unique(rtIntervals)], intensity=aggInts))
    }else{
        ## Check if there would be any need to aggregate the data.
        if(length(x@rtime) == length(unique(x@rtime))){
            ## No need
            return(cbind(rtime=rtimeOrdered(x),
                         intensity=intensityOrderedByRtime(x)))
        }else{
            ## Aggregate.
            aggVals <- unlist(lapply(split(intensityOrderedByRtime(x),
                                           rtimeOrdered(x)), FUN=FUN),
                              use.names=FALSE)
            return(cbind(rtime=unique(rtimeOrdered(x)), intensity=aggVals))
        }
    }
}

####============================================================
##  getChromatogram
##
##  Method to extract the chromatogram of an MSdata object. Data will be eventually
##  binned in rt dimension and intensities aggregated in m/z dimension.
##  The method returns a matrix with two columns, rtime and intensity.
setMethod("getChromatogram", "MSdata",
          function(object, FUN=max, bins=NULL, nbin=NULL, binSize=NULL, ...){
              ## Input arg checking... at least to some degree.
              if(!is.null(bins)){
                  if(!is.numeric(bins) | length(bins) < 2)
                      stop("'bins' should be a numeric vector of length > 2",
                           " specifying the bins in which the data should be binned.")
              }
              if(!is.null(nbin)){
                  if(!is.numeric(nbin) | length(nbin) > 1)
                      stop("'nbin' should be a numeric vector of length 1!")
              }
              if(!is.null(binSize)){
                  if(!is.numeric(binSize) | length(binSize) > 1)
                      stop("'binSize' should be a numeric vector of length 1!")
              }
              return(.getChrom(object,FUN=FUN, bins=bins, nbin=nbin,
                               binSize=binSize))
          })

####============================================================
##  plotChromatogram
##
####------------------------------------------------------------
setMethod("plotChromatogram", "MSdata",
          function(object, FUN=max, bins=NULL, nbin=NULL, binSize=NULL,
                   add=FALSE, main=paste(format(mzrange(object), 2), collapse="-"),
                   xlab="Retention time", ylab="Intensity",
                   ...){
              dat <- getChromatogram(object, FUN=FUN, bins=bins, nbin=nbin,
                                     binSize=binSize)
              if(add){
                  points(dat[, 1], dat[, 2], ...)
              }else{
                  plot(dat[, 1], dat[, 2], main=main, xlab=xlab, ylab=ylab, ...)
              }
          })
####============================================================
##  chromatogram
##
##  The same...
####------------------------------------------------------------
setMethod("chromatogram", "MSdata",
          function(object, FUN=max, bins=NULL, nbin=NULL, binSize=NULL,
                   add=FALSE, main=paste(format(mzrange(object), 2), collapse="-"),
                   xlab="Retention time", ylab="Intensity", ...){
              plotChromatogram(object=object, FUN=FUN, bins=bins, nbin=nbin,
                               binSize=binSize, add=add, main=main, xlab=xlab,
                               ylab=ylab, ...)
          })


####============================================================
##  getSpectrum
##
## Get the spectrum from an MSdata object. Data will be eventually binned
## in m/z dimension and intensities aggregated within these bins and along
## the rt dimension.
####------------------------------------------------------------
setMethod("getSpectrum", "MSdata",
          function(object, FUN=max, bins=NULL, nbin=NULL, binSize=NULL, ...){
              ## Input arg checking... at least to some degree.
              if(!is.null(bins)){
                  if(!is.numeric(bins) | length(bins) < 2)
                      stop("'bins' should be a numeric vector of length > 2",
                           " specifying the bins in which the data should be binned.")
              }
              if(!is.null(nbin)){
                  if(!is.numeric(nbin) | length(nbin) > 1)
                      stop("'nbin' should be a numeric vector of length 1!")
              }
              if(!is.null(binSize)){
                  if(!is.numeric(binSize) | length(binSize) > 1)
                      stop("'binSize' should be a numeric vector of length 1!")
              }
              return(.getSpec(object,FUN=FUN, bins=bins, nbin=nbin,
                              binSize=binSize))

          })
####============================================================
##  .getSpec
##
##  Basically the same as .getChrom, just along the other dimension.
##  We're extracting all values ordered by Mz in this function, i.e.
##  mzOrdered and intensityOrderedByMz.
####------------------------------------------------------------
.getSpec <- function(x, FUN=sum, bins=NULL, nbin=NULL, binSize=NULL){
    if(!is.null(nbin) & !is.null(binSize))
        stop("Arguments 'nbin' and 'binSize' are mutually exclusive")
    if(is.null(bins) & (!is.null(nbin) | !is.null(binSize))){
        ## Calculate the bins internally.
        bins <- .getBins(mzrange(x), nbin=nbin, binSize=binSize)
    }
    ## Do the aggregation on bins if it's not empty
    if(!is.null(bins)){
        ## Get the rt intervals for which we do have data.
        mzIntervals <- findInterval(mzOrdered(x), bins, all.inside=TRUE)
        aggInts <- .aggregateWithIntervals(intensityOrderedByMz(x), mzIntervals,
                                           FUN=FUN)
        return(cbind(mz=.binMid(bins)[unique(mzIntervals)], intensity=aggInts))
    }else{
        ## Check if there would be any need to aggregate the data.
        if(length(x@mz) == length(unique(x@mz))){
            ## No need
            return(cbind(mz=mzOrdered(x), intensity=intensityOrderedByMz(x)))
        }else{
            ## Aggregate.
            aggVals <- unlist(lapply(split(intensityOrderedByMz(x),
                                           mzOrdered(x)), FUN=FUN),
                              use.names=FALSE)
            return(cbind(mz=unique(mzOrdered(x)), intensity=aggVals))
        }
    }
}
####============================================================
##  plotSpectrum
##
####------------------------------------------------------------
setMethod("plotSpectrum", "MSdata",
          function(object, FUN=max, bins=NULL, nbin=NULL, binSize=NULL,
                   add=FALSE, main=paste(format(rtrange(object), 2), collapse="-"),
                   xlab="M/Z", ylab="Intensity",
                   ...){
              dat <- getSpectrum(object, FUN=FUN, bins=bins, nbin=nbin,
                                 binSize=binSize)
              if(add){
                  points(dat[, 1], dat[, 2], ...)
              }else{
                  plot(dat[, 1], dat[, 2], main=main, xlab=xlab, ylab=ylab, ...)
              }
          })


####============================================================
##  binMz
##
##  bin a MSdata object in M/Z dimension.
####------------------------------------------------------------
setMethod("binMz", "MSdata",
          function(object, FUN=max, bins=NULL, nbin=NULL, binSize=NULL){
              if(!is.null(nbin) & !is.null(binSize))
                  stop("Arguments 'nbin' and 'binSize' are mutually exclusive")
              ## Input arg checking... at least to some degree.
              if(!is.null(bins)){
                  if(!is.numeric(bins) | length(bins) < 2)
                      stop("'bins' should be a numeric vector of length > 2",
                           " specifying the bins in which the data should be binned.")
              }
              if(!is.null(nbin)){
                  if(!is.numeric(nbin) | length(nbin) > 1)
                      stop("'nbin' should be a numeric vector of length 1!")
              }
              if(!is.null(binSize)){
                  if(!is.numeric(binSize) | length(binSize) > 1)
                      stop("'binSize' should be a numeric vector of length 1!")
              }
              if(is.null(nbin) & is.null(bins) & is.null(binSize)){
                  ## Well, just return the object.
                  return(object)
              }else{
                  return(.binMz(object, FUN=FUN, bins=bins, nbin=nbin,
                                binSize=binSize))
              }
          })

####============================================================
##  binRtime
##
##  bin a MSdata object in retention time dimension
####------------------------------------------------------------
setMethod("binRtime", "MSdata",
          function(object, FUN=max, bins=NULL, nbin=NULL, binSize=NULL){
              if(!is.null(nbin) & !is.null(binSize))
                  stop("Arguments 'nbin' and 'binSize' are mutually exclusive.")
              ## Input arg checking... at least to some degree.
              if(!is.null(bins)){
                  if(!is.numeric(bins) | length(bins) < 2)
                      stop("'bins' should be a numeric vector of length > 2",
                           " specifying the bins in which the data should be binned.")
              }
              if(!is.null(nbin)){
                  if(!is.numeric(nbin) | length(nbin) > 1)
                      stop("'nbin' should be a numeric vector of length 1!")
              }
              if(!is.null(binSize)){
                  if(!is.numeric(binSize) | length(binSize) > 1)
                      stop("'binSize' should be a numeric vector of length 1!")
              }
              if(is.null(nbin) & is.null(bins) & is.null(binSize)){
                  ## Well, just return the object.
                  return(object)
              }else{
                  return(.binRtime(object, FUN=FUN, bins=bins, nbin=nbin,
                                   binSize=binSize))
              }
          })

####============================================================
##  binMzRtime
##
##  bin first in M/Z, then in retention time.
####------------------------------------------------------------
setMethod("binMzRtime", "MSdata",
          function(object, FUN=max, mzNbin=NULL, mzBinSize=NULL,
                   rtNbin=NULL, rtBinSize=NULL){
              ## Argument checking.
              if(!is.null(mzNbin) & !is.null(mzBinSize))
                  stop("Arguments 'mzNbin' and 'mzBinSize' are mutually exclusive.")
              if(!is.null(rtNbin) & !is.null(rtBinSize))
                  stop("Arguments 'rtNbin' and 'rtBinSize' are mutually exclusive.")
              if(!is.null(mzNbin)){
                  if(!is.numeric(mzNbin) | length(mzNbin) > 1)
                      stop("'mzNbin' should be a numeric vector of length 1!")
              }
              if(!is.null(mzBinSize)){
                  if(!is.numeric(mzBinSize) | length(mzBinSize) > 1)
                      stop("'mzBinSize' should be a numeric vector of length 1!")
              }
              if(!is.null(rtNbin)){
                  if(!is.numeric(rtNbin) | length(rtNbin) > 1)
                      stop("'rtNbin' should be a numeric vector of length 1!")
              }
              if(!is.null(rtBinSize)){
                  if(!is.numeric(rtBinSize) | length(rtBinSize) > 1)
                      stop("'rtBinSize' should be a numeric vector of length 1!")
              }
              ## Do the stuff.
              return(.binRtime(.binMz(object, nbin=mzNbin,
                                      binSize=mzBinSize, FUN=FUN, sort=FALSE),
                               , nbin=rtNbin, binSize=rtBinSize, FUN=FUN))
          })


####============================================================
##  mapMatrix
##
##  Returns a sparse matrix with rows being M/Z, columns being
##  retention time.
####------------------------------------------------------------
setMethod("mapMatrix", "MSdata",
          function(object){
              return(.msData2mapSparseMatrix(object))
          })

####============================================================
##  subset
##
##  Subset the MSdata object by mzrange and/or rtrange
####------------------------------------------------------------
setMethod("subset", "MSdata", function(x, mzrange=NULL, rtrange=NULL){
    if(is.null(mzrange) & is.null(rtrange)){
        return(x)
    }
    returnMeMz <- -9
    returnMeRt <- -9
    ## Checking mzrange input
    if(!is.null(mzrange)){
        if(!is.numeric(mzrange))
            stop("'mzrange' has to be a numeric vector of length 2.")
        if(length(mzrange) != 2)
            stop("'mzrange' has to be a numeric vector of length 2.")
        mzrange <- sort(mzrange)
        returnMeMz <- which(mz(x) >= mzrange[1] & mz(x) <= mzrange[2])
        ## Return an empty object if we've got no matching data.
        if(length(returnMeMz) == 0){
            return(MSdata(rtime=numeric(), mz=numeric(), intensity=numeric()))
        }
    }
    ## Checking rtrange input
    if(!is.null(rtrange)){
        if(!is.numeric(rtrange))
            stop("'rtrange' has to be a numeric vector of length 2.")
        if(length(rtrange) != 2)
            stop("'rtrange' has to be a numeric vector of length 2.")
        rtrange <- sort(rtrange)
        returnMeRt <- which(rtime(x) >= rtrange[1] & rtime(x) <= rtrange[2])
        ## Return an empty object if we've got no matching data.
        if(length(returnMeRt) == 0){
            return(MSdata(rtime=numeric(), mz=numeric(), intensity=numeric()))
        }
    }
    ## If we got so far we can /really/ subset the data.
    if(returnMeMz[1] > 0){
        ## OK we've got a subset for MZ
        if(returnMeRt[1] > 0){
            ## Got also retention time subset: make an intersect.
            returnMeMz <- intersect(returnMeMz, returnMeRt)
            if(length(returnMeMz) == 0)
                return(MSdata(rtime=numeric(), mz=numeric(), intensity=numeric()))
        }
        return(MSdata(rtime=rtime(x)[returnMeMz],
                      mz=mz(x)[returnMeMz],
                      intensity=intensity(x)[returnMeMz]))
    }else{
        ## Only rtime subset.
        return(MSdata(rtime=rtime(x)[returnMeRt],
                      mz=mz(x)[returnMeRt],
                      intensity=intensity(x)[returnMeRt]))
    }
})

####============================================================
####============================================================
##  "private" methods
##
####------------------------------------------------------------

####============================================================
##  rtimeOrdered
##
##  return the retention time as an ordered numeric.
####------------------------------------------------------------
setMethod("rtimeOrdered", "MSdata",
          function(object){
              return(sort(rtime(object)))
          })
####============================================================
##  mzOrdered
##
##  Return the mz values ordered increasingly.
####------------------------------------------------------------
setMethod("mzOrdered", "MSdata",
          function(object){
              return(sort(mz(object)))
          })
####============================================================
##  intensityOrderedByMz
##
##  Return the intensity values ordered by mz. That way, the mz
##  values returned by mzOrdered are in-sync (i.e. ordered in the
##  same way) with the intensity values.
####------------------------------------------------------------
setMethod("intensityOrderedByMz", "MSdata",
          function(object){
              idx <- order(mz(object))
              return(intensity(object)[idx])
          })
####============================================================
##  intensityOrderedByRtime
##
##  Return the intensity values ordered by retention time. That way,
##  the mz values returned by mzOrdered are in-sync (i.e. ordered in
##  the same way) with the intensity values.
####------------------------------------------------------------
setMethod("intensityOrderedByRtime", "MSdata",
          function(object){
              idx <- order(rtime(object))
              return(intensity(object)[idx])
          })

####============================================================
##  .msData2mapMatrix
##
##  Transforms the data in an MSdata object into a two-dimensional
##  "map" matrix with rows representing mz values, columns retention
##  times. The matrix will be a i x j matrix with i being the length
##  of unique mz values and j the length of unique retention times.
##  x: MSdata object.
####------------------------------------------------------------
.msData2mapMatrix <- function(x){
    ## Eventually do some sort of binning...
    ## x is supposed to be ordered by rt, but it shouldn't make
    ## a difference in the end...
    tmpL <- split(data.frame(mz=mz(x), intensity=intensity(x)),
                  f=rtime(x))
    mzUnique <- sort(unique(mz(x)))
    mzL <- length(mzUnique)
    tmp <- rep(NA, mzL)  ## Can re-use this. R will copy this variable.
    vals <- unlist(lapply(tmpL, function(z){
        ## Note: tmp is a copy, so any changes to tmp will not be propagated or be
        ## visible outside of the lapply.
        tmp[match(z$mz, mzUnique)] <- z$intensity
        return(tmp)
    }), use.names=FALSE)
    retMat <- matrix(vals, nrow=length(mzUnique))
    rownames(retMat) <- mzUnique
    colnames(retMat) <- unique(rtimeOrdered(x))
    return(retMat)
}

####============================================================
##  .msData2mapSparseMatrix
##
##  Transforms the data of an MSdata object into a 2-dimensional
##  (map) matrix in SparseMatrix format.
.msData2mapSparseMatrix <- function(x){
    ## the j (column) index will be each rt value repeated unique
    ## mz value times.
    mzUnique <- unique(mzOrdered(x))
    rtUnique <- unique(rtimeOrdered(x))
    iIdx <- unlist(lapply(split(mz(x), rtime(x)),
                          function(z){
                              return(match(z, mzUnique))
                          }
                          ), use.names=FALSE)
    ## jIdx <- unlist(lapply(split(rtime(x), mz(x)),
    ##                       function(z){
    ##                           return(match(z, rtUnique))
    ##                       }
    ##                       ), use.names=FALSE)
    jIdx <- match(rtime(x), rtUnique)
    spM <- sparseMatrix(i=iIdx, j=jIdx, x=intensity(x), dimnames=list(mzUnique, rtUnique))
    return(spM)
}

####============================================================
##  .binMz
##
##  bin the mz range. this means intensities falling within a mz bin
##  THAT HAVE THE SAME RETENTION TIME are aggregated.
##  x is supposed to be a matrix returned by as.matrix(MSdata)
##  NOTE: we're assuming all the argument checking has been performed
##        already, this means, only one of bins, nbin or binSize is
##        defined.
##  sort: if sort is FALSE we're not sorting by rt and mz, otherwise we do.
####------------------------------------------------------------
.binMz <- function(x, bins=NULL, nbin=NULL, binSize=NULL, FUN=max, sort=TRUE){
    ## That's definitely not the ideal way...
    ## theDf <- data.frame(rtime=rtime(x), mz=mz(x), intensity=intensity(x))
    ## binning on mz, splitting by the bins, and within each sub-data.frame split
    ## again but by retention time and aggregate the intensities.

    ## Alternative: define the splitting factor as a combination of bin and rt.
    if(is.null(bins)){
        bins <- .getBins(mzrange(x), nbin=nbin, binSize=binSize)
    }
    mzInts <- findInterval(mz(x), bins, all.inside=TRUE)
    ## Define the factors to split the intensities.
    splitF <- factor(paste0(mzInts, ":", rtime(x)))
    ## Extract the bin idx and the retention time from the splitF
    idMapMat <- matrix(as.numeric(do.call(rbind, strsplit(levels(splitF), split=":"))), ncol=2)
    ## The ordering... have to fix that, order by rtime, then by mz.
    ## The data is supposedly already sorted by mz, so just ordering by rt.
    if(sort){
        idxRt <- order(idMapMat[, 2], idMapMat[, 1])
        idMapMat <- idMapMat[idxRt, ]
    }else{
        idxRt <- 1L:nrow(idMapMat)
    }
    return(MSdata(mz=.binMid(bins)[idMapMat[, 1]],
                  rtime=idMapMat[, 2],
                  intensity=.aggregateWithIntervals(intensity(x), splitF, FUN=FUN)[idxRt]
                  ))
}
## Same as above, but does not any sorting/ordering of data by rt then by mz.
.binMzUnsorted <- function(x, bins=NULL, nbin=NULL, binSize=NULL, FUN=max){
    ## That's definitely not the ideal way...
    ## theDf <- data.frame(rtime=rtime(x), mz=mz(x), intensity=intensity(x))
    ## binning on mz, splitting by the bins, and within each sub-data.frame split
    ## again but by retention time and aggregate the intensities.

    ## Alternative: define the splitting factor as a combination of bin and rt.
    if(is.null(bins)){
        bins <- .getBins(mzrange(x), nbin=nbin, binSize=binSize)
    }
    mzInts <- findInterval(mz(x), bins, all.inside=TRUE)
    ## Define the factors to split the intensities.
    splitF <- factor(paste0(rtime(x), ":", mzInts))
    ## Extract the bin idx and the retention time from the splitF
    idMapMat <- do.call(rbind, strsplit(levels(splitF), split=":"))
    return(MSdata(mz=.binMid(bins)[as.numeric(idMapMat[, 2])],
                  rtime=as.numeric(idMapMat[, 1]),
                  intensity=.aggregateWithIntervals(intensity(x), splitF, FUN=FUN)
                  ))
}
## The same method, just with some slower code... to check that it works.
.binMzSlow <- function(x, bins=NULL, nbin=NULL, binSize=NULL, FUN=max){
    if(is.null(bins)){
        bins <- .getBins(mzrange(x), nbin=nbin, binSize=binSize)
    }
    mzInts <- findInterval(mz(x), bins, all.inside=TRUE)
    ## Make a data.frame of the data.
    datDfL <- split(data.frame(rtime=rtime(x), mz=mz(x), intensity=intensity(x)),
                    f=mzInts)
    ## Now, within each element, split by retention time and aggregate the intensities
    resL <- lapply(datDfL, function(z){
        tmp <- lapply(split(z[, c("rtime", "intensity")], f=z$rtime), function(zx, theFun){
            return(c(zx[1, "rtime"], theFun(zx[, "intensity"])))
        }, theFun=FUN)
        return(do.call(rbind, tmp))
    })
    ## Now creating the "final" data.
    newMz <- .binMid(bins)[rep(as.numeric(names(resL)), unlist(lapply(resL, FUN=nrow)))]
    otherDat <- do.call(rbind, resL)
    return(MSdata(mz=newMz,
                  rtime=otherDat[, 1],
                  intensity=otherDat[, 2]))
}


####============================================================
##  .binRtime
##
##  same as .binMz, just for retention time.
####------------------------------------------------------------
.binRtime <- function(x, bins=NULL, nbin=NULL, binSize=NULL, FUN=max, sort=TRUE){
    if(is.null(bins)){
        bins <- .getBins(rtime(x), nbin=nbin, binSize=binSize)
    }
    rtInts <- findInterval(rtime(x), bins, all.inside=TRUE)
    ## Define the factors to split the intensities.
    splitF <- factor(paste0(rtInts, ":", mz(x)))
    ## Extract the bin idx and the retention time from the splitF
    idMapMat <- matrix(as.numeric(do.call(rbind, strsplit(levels(splitF), split=":"))), ncol=2)
    if(sort){
        idxRt <- order(idMapMat[, 1], idMapMat[, 2])
        idMapMat <- idMapMat[idxRt, ]
    }else{
        idxRt <- 1L:nrow(idMapMat)
    }
    return(MSdata(rtime=.binMid(bins)[idMapMat[, 1]],
                  mz=idMapMat[, 2],
                  intensity=.aggregateWithIntervals(intensity(x), splitF, FUN=FUN)[idxRt]
                  ))
}
## The same method, just with some slower code... to check that it works.
.binRtimeSlow <- function(x, bins=NULL, nbin=NULL, binSize=NULL, FUN=max){
    if(is.null(bins)){
        bins <- .getBins(rtime(x), nbin=nbin, binSize=binSize)
    }
    rtInts <- findInterval(rtime(x), bins, all.inside=TRUE)
    ## Make a data.frame of the data.
    datDfL <- split(data.frame(rtime=rtime(x), mz=mz(x), intensity=intensity(x)),
                    f=rtInts)
    ## Now, within each element, split by mz and aggregate the intensities
    resL <- lapply(datDfL, function(z){
        tmp <- lapply(split(z[, c("mz", "intensity")], f=z$mz), function(zx, theFun){
            return(c(zx[1, "mz"], theFun(zx[, "intensity"])))
        }, theFun=FUN)
        return(do.call(rbind, tmp))
    })
    ## Now creating the "final" data.
    newRt <- .binMid(bins)[rep(as.numeric(names(resL)), unlist(lapply(resL, FUN=nrow)))]
    otherDat <- do.call(rbind, resL)
    return(MSdata(rtime=newRt,
                  mz=otherDat[, 1],
                  intensity=otherDat[, 2]))
}

