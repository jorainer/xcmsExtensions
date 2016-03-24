####============================================================
##  MSsliceList
##
####------------------------------------------------------------
MSsliceList <- function(...){
    slices <- list(...)
    if(length(slices) > 0){
        slices <- as.list(unlist(slices))
        obj <- new("MSsliceList", slices=slices)
    }else{
        obj <- new("MSsliceList")
    }
    validObject(obj)
    return(obj)
}
.validateMSsliceList <- function(object){
    ## Ensure that, if not empty, all elements in slices are
    ## MSslice objects!
    sl <- slices(object)
    if(length(sl) > 0){
        if(!all(unlist(lapply(sl, function(z){
            return(is(z, "MSslice"))
        }))))
            return("Only MSslice objects allowed in slot @slices!")
    }
    return(TRUE)
}
setValidity("MSsliceList", .validateMSsliceList)
setMethod("initialize", "MSsliceList", function(.Object, ...){
    OK <- .validateMSsliceList(.Object)
    if(is(OK, "character"))
        stop(OK)
    callNextMethod(.Object, ...)
})


####============================================================
##  MSslice
##
####------------------------------------------------------------
MSslice <- function(...){
    data <- list(...)
    if(length(data) > 0){
        data <- as.list(unlist(data))
        obj <- msSlice(data)
    }else{
        obj <- new("MSslice")
    }
    validObject(obj)
    return(obj)
}
setMethod("msSlice", "list", function(object, call=match.call(), ...){
    ## Have to evaluate that all elements are MSdata objects.
    if(length(object) > 0){
        if(!all(unlist(lapply(object, function(z){
            return(is(z, "MSdata"))
        }))))
            stop("Argument 'object' has to be a list of MSdata objects!")
        ## Define the rtrange and mzrange
        mzrs <- lapply(object, mzrange)
        rtrs <- lapply(object, rtrange)
        rtrange <- range(unlist(rtrs))
        mzrange <- range(unlist(mzrs))
    }else{
        rtrange <- numeric()
        mzrange <- numeric()
    }
    res <- new("MSslice", data=object, rtrange=rtrange, mzrange=mzrange,
               call=call)
    validObject(res)
    return(res)
})
.validateMSslice <- function(object){
    if(!all(unlist(lapply(object@data, function(z){
        return(is(z, "MSdata"))
    }))))
        return("Only MSdata objects allowed in slot @data!")
    return(TRUE)
}
setValidity("MSslice", .validateMSslice)
setMethod("initialize", "MSslice", function(.Object, ...){
    OK <- .validateMSslice(.Object)
    if(is(OK, "character"))
        stop(OK)
    callNextMethod(.Object, ...)
})

####============================================================
##  MSdata
##
####------------------------------------------------------------
MSdata <- function(mz, rtime, intensity, mslevel=1){
    if(missing(mz) | missing(rtime) | missing(intensity))
        stop("All three of 'mz', 'rtime' and 'intensity' have to be specified!")
    if(length(rtime) == 0){
        rtr <- c(0, 0)
    }else{
        rtr <- range(rtime)
    }
    if(length(mz) == 0){
        mzr <- c(0, 0)
    }else{
        mzr <- range(mz)
    }
    if(length(intensity) == 0){
        intr <- c(0, 0)
    }else{
        intr <- range(intensity)
        }
    res <- new("MSdata", mz=as.numeric(mz), rtime=Rle(rtime),
               intensity=as.integer(intensity), mslevel=mslevel, rtrange=rtr,
               mzrange=mzr, intrange=intr)
    return(res)
}
.validateMSdata <- function(object){
    ## mz, rtime and intensity all have to have the same length!
    if(length(unique(c(length(object@mz), length(object@rtime),
                       length(object@intensity))))!=1)
        return(paste0("'mz', 'rtime' and 'intensity' have different lengths!"))
    return(TRUE)
}
setValidity("MSdata", .validateMSdata)
setMethod("initialize", "MSdata", function(.Object, ...){
    OK <- .validateMSdata(.Object)
    if(is(OK, "character"))
        stop(OK)
    callNextMethod(.Object, ...)
})
