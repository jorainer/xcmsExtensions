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
##  Create an msSlice object from the MSdata.
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
