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
    return(object@rtrange)
})

####============================================================
##  mzrange
##
##  Getter the rtrange slot.
####------------------------------------------------------------
setMethod("mzrange", "MSdata", function(object){
    return(object@mzrange)
})

####============================================================
##  intrange
##
##  Getter the rtrange slot.
####------------------------------------------------------------
setMethod("intrange", "MSdata", function(object){
    if(length(object@intensity) == 0){
        return(c(0, 0))
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
