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

setGeneric("binMz", function(object, ...)
    standardGeneric("binMz"))
setGeneric("binRtime", function(object, ...)
    standardGeneric("binRtime"))
setGeneric("binMzRtime", function(object, ...)
    standardGeneric("binMzRtime"))

setGeneric("mapMatrix", function(object, ...)
    standardGeneric("mapMatrix"))

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

setGeneric("peakGroupSummary", function(object, ...)
    standardGeneric("peakGroupSummary"))

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

## ## Remove that stuff for Bioc 3.3
setGeneric("condition<-", function(x, value)
    standardGeneric("condition<-"))
setGeneric("value", function(x, db, ...)
    standardGeneric("value"))
setGeneric("value<-", function(x, value)
    standardGeneric("value<-"))


## Database stuff:
## That one could be imported from ensembldb
if(!isGeneric("listTables")){
    setGeneric("listTables", function(x, ...)
        standardGeneric("listTables"))
}
setGeneric("compounds", function(x, columns, filter=list(), ...)
    standardGeneric("compounds"))
setGeneric("cleanColumns", function(x, columns, ...)
    standardGeneric("cleanColumns"))
setGeneric("cleanTables", function(x, tables, ...)
    standardGeneric("cleanTables"))
setGeneric("cleanFilter", function(x, filter, ...)
    standardGeneric("cleanFilter"))
setGeneric("buildQuery", function(x, ...)
    standardGeneric("buildQuery"))
setGeneric("prefixColumns", function(x, columns, clean=TRUE, with.tables)
    standardGeneric("prefixColumns"))
setGeneric("buildJoinQuery", function(x, ...)
    standardGeneric("buildJoinQuery"))
setGeneric("buildFilterQuery", function(x, ...)
    standardGeneric("buildFilterQuery"))
setGeneric("sortTablesByDegree", function(x, tables, ...)
    standardGeneric("sortTablesByDegree"))
setGeneric("addRequiredJoinTables", function(x, tables)
    standardGeneric("addRequiredJoinTables"))
setGeneric("getWhat", function(x, ...)
    standardGeneric("getWhat"))


## Utils and alike
setGeneric("grow", function(x, by, ...)
    standardGeneric("grow"))
setGeneric("mzmatch", function(x, mz, mzdev=0, ppm=10, ...)
    standardGeneric("mzmatch"))
setGeneric("mass2adductmz", function(x, ionAdduct=NULL, ...)
    standardGeneric("mass2adductmz"))
setGeneric("adductmz2mass", function(x, ionAdduct=NULL, ...)
    standardGeneric("adductmz2mass"))
setGeneric("supportedIonAdducts", function(x, charge)
    standardGeneric("supportedIonAdducts"))

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
