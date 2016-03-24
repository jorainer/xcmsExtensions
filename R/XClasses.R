setClassUnion("MatrixOrNull", c("matrix", "NULL", "missing"))
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

