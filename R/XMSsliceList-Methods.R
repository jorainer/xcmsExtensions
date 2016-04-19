####============================================================
##  Methods for MSsliceList
##
####------------------------------------------------------------

####============================================================
##  show
##
####------------------------------------------------------------
setMethod("show", "MSsliceList", function(object){
    cat(class(object), " object with ", length(object), " slices:\n\n", sep="")
    if(length(object) > 4){
        idxs <- c(1, 2, length(object)-1, length(object))
        for(i in 1:2){
            cat("[[", i,"]]\n", sep="")
            show(slices(object)[[i]])
            cat("\n")
        }
        cat("   ...   \n\n")
        for(i in (length(object)-1):length(object)){
            cat("[[", i,"]]\n", sep="")
            show(slices(object)[[i]])
            cat("\n")
        }
    }else{
        show(slices(object))
    }
})


####============================================================
##  slices
##
##  Returns the content of the slices slot.
####------------------------------------------------------------
setMethod("slices", "MSsliceList", function(object){
    return(object@slices)
})
setReplaceMethod("slices", "MSsliceList", function(object, value){
    if(!is(value, "list")){
        value <- list(value)
    }
    object@slices <- value
    validObject(object)
    return(object)
})

####============================================================
##  length
##
####------------------------------------------------------------
setMethod("length", "MSsliceList", function(x){
    return(length(slices(x)))
})

####============================================================
##  rtranges
##
##  Get the rtrange of each MSslice object.
####------------------------------------------------------------
setMethod("rtranges", "MSsliceList", function(object){
    return(lapply(slices(object), rtrange))
})

####============================================================
##  rtrange
##
##  Get the range of rtranges of each MSslice object.
####------------------------------------------------------------
setMethod("rtrange", "MSsliceList", function(object){
    return(range(rtranges(object), na.rm=TRUE))
})

####============================================================
##  mzranges
##
##  Get the rtrange of each MSslice object.
####------------------------------------------------------------
setMethod("mzranges", "MSsliceList", function(object){
    return(lapply(slices(object), mzrange))
})

####============================================================
##  mzrange
##
##  Get the range of rtranges of each MSslice object.
####------------------------------------------------------------
setMethod("mzrange", "MSsliceList", function(object){
    return(range(mzranges(object), na.rm=TRUE))
})

####============================================================
##  intranges
##
##  Get the rtrange of each MSslice object.
####------------------------------------------------------------
setMethod("intranges", "MSsliceList", function(object){
    return(lapply(slices(object), intrange))
})

####============================================================
##  intrange
##
##  Get the range of rtranges of each MSslice object.
####------------------------------------------------------------
setMethod("intrange", "MSsliceList", function(object){
    return(range(intranges(object), na.rm=TRUE))
})


####============================================================
##  [
##
##  Single bracket subsetting.
####------------------------------------------------------------
## [ subsetting.
.bracketSubset <- function(x, i, j, ..., drop){
    if(!missing(j))
        stop("Subsetting by columns ('j') is not supported.")
    ## Check if index i is outside of what we've got.
    haveEls <- length(slices(x))
    if(haveEls == 0){
        warning("Can not subset an empty object.")
        return(x)
    }
    if(missing(i))
        i <- 1:haveEls
    i <- .checkElementIndices(i, haveEls, NULL)
    slices(x) <- slices(x)[i]
    return(x)
}
.checkElementIndices <- function(i, elNo, elNames){
    if(is.character(i)){
        if(is.null(elNames))
            stop("Can not subset by names as the object has no names!")
        ## match is to elNames.
        if(!all(i %in% elNames)){
            stop("Element(s) with name(s) ", paste(i[!(i %in% elNames)], collapse=", "),
                 " not found!")
        }
        return(match(i, elNames))
    }
    if(is.logical(i)){
        if(length(i) != elNo)
            stop("Length of 'i' does not match number of elements in 'x'!",
                 " If 'i' is a logical vector its length has to match that of 'x'.")
        i <- which(i)
    }
    ## Check if we've got elements outside...
    outs <- sum(!(i %in% 1:elNo))
    if(outs > 0)
        stop(outs, " indices are larger than the number of available slices!")
    return(i)
}
setMethod("[", "MSsliceList", .bracketSubset)

####============================================================
##  [[
##
##  Extract a single MSslice object.
####------------------------------------------------------------
setMethod("[[", "MSsliceList", function(x, i, j, ...){
    haveEls <- length(slices(x))
    if(!missing(j))
        stop("Subsetting by 'j' ([[ i, j ]] is not supported!")
    if(missing(i)){
        return(x)
    }
    if(length(i) > 1)
        stop("Can only extract a single element with [[,  but ", length(i),
             " indices were submitted.")
    i <- .checkElementIndices(i, haveEls, NULL)
    return(slices(x)[[i]])
})
