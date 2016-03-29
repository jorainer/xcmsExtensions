setGeneric("closestGroup", function(x, ...)
    standardGeneric("closestGroup"))

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
