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
