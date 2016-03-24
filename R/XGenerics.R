setGeneric("closestGroup", function(x, mz=NULL, rt=NULL, ...)
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

setGeneric("msSlice", function(object, mzrange=NULL, rtrange=NULL, ...)
    standardGeneric("msSlice"))

setGeneric("slices", function(object, ...)
    standardGeneric("slices"))
setGeneric("slices<-", function(object, value)
    standardGeneric("slices<-"))
