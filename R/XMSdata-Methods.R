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
    if(is.null(bins) & (!is.null(nbin) | !is.null(binSize)){
        ## Calculate the bins internally.
        bins <- .getBins()
    }

    ## Do the aggregation on bins if it's not empty
    if(!is.null(bins)){
    }else{
        ## Check if there would be any need to aggregate the data.
        if(length(msd@rtime) == length(unique(msd@rtime))){
            ## No need
        }else{
            ## Aggregate.
        }
    }
}
