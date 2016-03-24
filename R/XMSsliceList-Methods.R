####============================================================
##  Methods for MSsliceList
##
####------------------------------------------------------------

####============================================================
##  show
##
####------------------------------------------------------------
setMethod("show", "MSsliceList", function(object){
    cat(class(object), " object:\n", sep="")
    cat("| Number of slices: ", length(object), "\n", sep="")
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

## To implement:
## + rtranges, mzranges, intranges.


