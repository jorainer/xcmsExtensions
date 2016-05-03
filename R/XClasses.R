setClassUnion("MatrixOrNull", c("matrix", "NULL", "missing"))
setClassUnion("MatrixOrNumericOrNull", c("matrix", "numeric", "NULL"))
setClassUnion("NumericOrNull", c("numeric", "integer","NULL", "missing"))
setClassUnion("NumericOrMissing", c("numeric", "integer","missing"))
setClassUnion("LogicalOrMissing", c("logical","missing"))
setClassUnion("AssayData", c("list", "environment"))
####============================================================
##  MSslice
##
##  Class representing a slice from the MS data (mz/rt/intentsity).
##  The class allows to store data from several files, but defined
##  by the same region.
####------------------------------------------------------------
setClass("MSsliceOld",
         representation(data="list",
                        call="call",
                        mzrange="numeric",
                        rtrange="numeric",
                        names="character"),
         prototype(
             data=list(),
             call=call("new"),
             mzrange=numeric(),
             rtrange=numeric(),
             names=character()
         )
         )

## Test what happens if we extend the eSet; it's not build for this type
## of data, but, well, so what.
## Methods implemented are:
## pData
## experimentData
## assayData
## all other "benefits" form the eSet method.
setClass("MSslice",
         representation=representation(
             experimentData="MIAxE",
             mzrange="numeric",
             rtrange="numeric",
             phenoData="AnnotatedDataFrame",
             assayData="AssayData"
         ),
         contains="Versioned",
         prototype=prototype(
             new("VersionedBiobase",
                 versions=c(MSslice="0.0.1")),
             experimentData=new("MIAPE"),
             assayData=list(),
             annotation="No feature annotation.",
             mzrange=numeric(),
             rtrange=numeric()))


####============================================================
##  MSsliceList
##
##  A collection of slices defined by different ranges.
####------------------------------------------------------------
setClass("MSsliceList",
         representation(slices="list"),
         prototype(slices=list()))

####============================================================
##  MSdata
##
##  Class containing the acutal (raw data) from a MS scan, i.e.
##  the intensity for each mz/retention time tuple.
####------------------------------------------------------------
setClass("MSdata",
         representation(mslevel="numeric",
                        mz="numeric",
                        intensity="integer",
                        rtime="Rle",
                        rtrange="numeric",
                        mzrange="numeric",
                        intrange="numeric"),
         prototype(
             mslevel=1,
             mz=numeric(),
             intensity=integer(),
             rtime=Rle(),
             rtrange=numeric(),
             mzrange=numeric(),
             intrange=numeric())
         )


####============================================================
##  SimpleCompoundDb
##
##  Amazing and crazy metabolite/compound database. By design this
##  should in the end hold some basic data from a lot of different
##  dabases.
####------------------------------------------------------------
setClass("SimpleCompoundDb",
         representation(con="DBIConnection", tables="list", .properties="list"),
         prototype=list(con=NULL, .properties=list(), tables=list()))


####============================================================
##  CompoundidFilter
##
####------------------------------------------------------------
setClass("CompoundidFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
            )
        )

####============================================================
##  BasicrangeFilter
##
##  That's a generic range filter. Just important that the validity method
##  checks a) that value is numeric, b) that condition is () (] [] [) and
##  c) that value has length 2.
.validBasicrangeFilter <- function(object){
    vals <- object@value
    if(length(vals) != 2)
        return(paste0("Length of 'value' has to be 2!"))
    suppressWarnings(
        vals <- as.numeric(vals)
    )
    if(any(is.na(vals)))
        return(paste0("'value' should only contain numeric values!"))
    ## Check condition.
    if(!(object@condition %in% c("()", "[)", "[]", "(]")))
        return(paste0("'condition' can only take values '()', '(]', '[]' or '[)'."))
    return(TRUE)
}
setClass("BasicrangeFilter", contains="BasicFilter",
         prototype=list(
             condition="[]",
             value=c("-Inf", "Inf"),
             .valueIsCharacter=FALSE
         ),
         validity=.validBasicrangeFilter
         )
##  initialize method.
## We're overwriting the method to avoid calling the initialize method of the BasicFilter,
## that wouldn't allow the conditions we need here. Thus we're NOT calling the
## callNextMethod as we would usually do to fill in all slots, but have to do this
## manually here. All classes extending this class do then also call this method.
setMethod("initialize", "BasicrangeFilter", function(.Object, value, condition){
    ## .Object <- callNextMethod()  ## Can not call that as it would call the validate of
    ## BasicFilter.
    if(!missing(value))
        .Object@value <- as.character(value)
    if(!missing(condition))
        .Object@condition <- as.character(condition)
    OK <- .validBasicrangeFilter(.Object)
    if(is(OK, "character"))
        stop(OK)
    ## Oh, man, doing the initialization here; otherwise I would have to update
    ## the ensembldb implementation!
    return(.Object)
    ##callNextMethod(.Object, ...)
})


####============================================================
##  MassrangeFilter
##
##  Extend BasicrangeFilter, adding also a slot column and over-
##  writing the corresponding methods.
##  The validity method we're using here ensures that the superclass is OK;
##  we can't call validObject, since that starts with the BasicFilter!
####------------------------------------------------------------
.validMassrangeFilter <- function(object){
    ## Call the method from the BasicrangeFilter
    return(getValidity(getClassDef("BasicrangeFilter"))(object))
}
setClass("MassrangeFilter", contains="BasicrangeFilter",
         representation(column="character"),
         prototype=list(
             column="mass"
         ),
         validity=.validMassrangeFilter)
setMethod("initialize", "MassrangeFilter", function(.Object, value, condition,
                                                    column){
    .Object <- callNextMethod(.Object, value, condition)
    if(!missing(column))
        .Object@column <- column
    return(.Object)
})
