setGeneric("closestGroup", function(x, ...)
    standardGeneric("closestGroup"))

setGeneric("getChromatogram", function(object, ...)
    standardGeneric("getChromatogram"))
setGeneric("getSpectrum", function(object, ...)
    standardGeneric("getSpectrum"))
setGeneric("plotChromatogram", function(object, ...)
    standardGeneric("plotChromatogram"))
setGeneric("plotSpectrum", function(object, ...)
    standardGeneric("plotSpectrum"))

setGeneric("getData", function(x, mzrange=NULL, rtrange=NULL, intrange=NULL, ...)
    standardGeneric("getData"))

setGeneric("scantimes", function(x)
    standardGeneric("scantimes"))

setGeneric("rtime<-", function(object, value)
    standardGeneric("rtime<-"))

setGeneric("intensity<-", function(object, value)
    standardGeneric("intensity<-"))

setGeneric("intrange", function(object, ...)
    standardGeneric("intrange"))

setGeneric("msData", function(object, ...)
    standardGeneric("msData"))

setGeneric("msSlice", function(object, ...)
    standardGeneric("msSlice"))

## For MSsliceList
setGeneric("slices", function(object, ...)
    standardGeneric("slices"))
setGeneric("slices<-", function(object, value)
    standardGeneric("slices<-"))
setGeneric("mzranges", function(object, ...)
    standardGeneric("mzranges"))
setGeneric("rtranges", function(object, ...)
    standardGeneric("rtranges"))
setGeneric("intranges", function(object, ...)
    standardGeneric("intranges"))

####============================================================
##  internal methods.
setGeneric("rtimeOrdered", function(object, ...)
    standardGeneric("rtimeOrdered"))
setGeneric("mzOrdered", function(object, ...)
    standardGeneric("mzOrdered"))
setGeneric("intensityOrderedByMz", function(object, ...)
    standardGeneric("intensityOrderedByMz"))
setGeneric("intensityOrderedByRtime", function(object, ...)
    standardGeneric("intensityOrderedByRtime"))
