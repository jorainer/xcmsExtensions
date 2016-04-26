###
### Load any db objects whenever the package is loaded.
###

.onLoad <- function(libname, pkgname)
{
    options(useFancyQuotes=FALSE)
    ns <- asNamespace(pkgname)
    path <- system.file("extdata", package=pkgname, lib.loc=libname)
    files <- dir(path, ".sqlite$")
    for(i in seq_len(length(files))){
        db <- SimpleCompoundDb(system.file("extdata", files[[i]], package=pkgname,
                                           lib.loc=libname))
        objname <- sub(".sqlite$","",files[[i]])
        assign(objname, db, envir=ns)
        namespaceExport(ns, objname)
    }
    ## Read the Fiehn lab ESI-adducts.
    ESIADDUCTS <- read.table(system.file("txt/ESI-MS-adducts.txt", package=pkgname, lib.loc=libname),
                             header=TRUE, as.is=TRUE, sep="\t")
    assign("ESIADDUCTS", ESIADDUCTS, envir=ns)
}

