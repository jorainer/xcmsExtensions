setClassUnion("MatrixOrNull", c("matrix", "NULL", "missing"))
setClassUnion("MatrixOrNumericOrNull", c("matrix", "numeric", "NULL"))
setClassUnion("NumericOrNull", c("numeric", "integer","NULL", "missing"))
setClassUnion("NumericOrMissing", c("numeric", "integer","missing"))
setClassUnion("LogicalOrMissing", c("logical","missing"))

####============================================================
##  MSslice
##
##  Class representing a slice from the MS data (mz/rt/intentsity).
##  The class allows to store data from several files, but defined
##  by the same region.
####------------------------------------------------------------
setClass("MSslice",
         representation(data="list",
                        call="call",
                        mzrange="numeric",
                        rtrange="numeric"),
         prototype(
             data=list(),
             call=call("new"),
             mzrange=numeric(),
             rtrange=numeric()
         )
         )

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
