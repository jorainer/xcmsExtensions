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
MSslice <- function(x, ...){
    if(missing(x))
        return(new("MSslice"))
    if(is(x, "MSdata"))
        x <- list(x)
    if(!is(x, "list"))
        stop("Argument 'x' has to be a list of MSdata objects!")
##    data <- list(...)
    if(length(x) > 0){
        ## x <- as.list(unlist(data))
        obj <- msSlice(x, ...)
    }else{
        obj <- new("MSslice")
    }
    validObject(obj)
    return(obj)
}
setMethod("msSlice", "list", function(object, ...){
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
    ## res <- new("MSslice", data=object, rtrange=rtrange, mzrange=mzrange,
    ##            call=call)
    res <- new("MSslice", assayData=object, rtrange=rtrange, mzrange=mzrange, ...)
    validObject(res)
    return(res)
})
.validateMSslice <- function(object){
    ## data <- object@data
    data <- assayData(object)
    if(length(data) > 0){
        if(!all(unlist(lapply(data, function(z){
            return(is(z, "MSdata"))
        }))))
            return("Only MSdata objects allowed in slot addayData!")
    }
    ## Check the phenoData
    pd <- phenoData(object)
    ## nrow pd has to match length of assayData!
    if(nrow(pd) > 0){
        if(nrow(pd) != length(data))
            return("The number of rows of the pheno data does not match the number of MSdata objects!")
    }

    ## ## Check the names slot
    ## if(!is.null(names(object))){
    ##     if(length(object@data) != length(object@names))
    ##         return("The number of names does not match the number of MSdata objects!")
    ## }
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
        rtr <- numeric()
    }else{
        rtr <- range(rtime)
    }
    if(length(mz) == 0){
        mzr <- numeric()
    }else{
        mzr <- range(mz)
    }
    if(length(intensity) == 0){
        intr <- numeric()
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


####============================================================
##  SimpleCompoundDb
##
##  SimpleCompoundDb constructor.
####------------------------------------------------------------
SimpleCompoundDb <- function(x){
    ## x is supposed to be the file name of the SQLite database!
    if(missing(x))
        stop("No SQLite database file provided!")
    lite <- dbDriver("SQLite")
    con <- dbConnect(lite, dbname=x, flags=SQLITE_RO)
    db <- new("SimpleCompoundDb", con=con)
    ## Get the tables and store that.
    theTab <- .doListTables(db)
    db@tables <- theTab
    return(db)
}
.validateSimpleCompoundDb <- function(object){
    if(!is.null(object@con)){
        con <- object@con
        reqTab <- c("compound_basic", "metadata")
        gotTab <- dbListTables(con)
        if(!all(reqTab %in% gotTab)){
            errStr <- paste0("The SQLite database does not provide the required",
                             " tables ", paste(sQuote(reqTab), collapse=", "), "!")
            return(errStr)
        }
    }
    return(TRUE)
}
setValidity("SimpleCompoundDb", .validateSimpleCompoundDb)
setMethod("initialize", "SimpleCompoundDb", function(.Object, ...){
    OK <- .validateSimpleCompoundDb(.Object)
    if(is(OK, "character"))
        stop(OK)
    callNextMethod(.Object, ...)
})


####============================================================
##  CompoundidFilter
##
####------------------------------------------------------------
CompoundidFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("CompoundidFilter", condition=condition, value=as.character(value)))
}

####============================================================
##  MassrangeFilter
##
####------------------------------------------------------------
MassrangeFilter <- function(value, condition="[]", column="mass"){
    obj <- new("MassrangeFilter", value=value, condition=condition, column=column)
    return(obj)
}
## setValidity("MassrangeFilter", .validateBasicrangeFilter)
## setMethod("initialize", "MassrangeFilter", function(.Object, ...){
##     OK <- .validateBasicrangeFilter(.Object)
##     if(is(OK, "character"))
##         stop(OK)
##     return(.Object)
##     ##callNextMethod(.Object, ...)
## })
