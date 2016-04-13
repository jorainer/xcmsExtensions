detach("package:xcmsExtensions", unload=TRUE)
library(xcmsExtensions)

test_parseXml <- function(){
    xf <- system.file("extdata/HMDB00001.xml", package="xcmsExtensions")
    res <- xcmsExtensions:::.simpleParseHmdbXml(xf)
    checkEquals(res[1, "accession"], "HMDB00001")
    checkEquals(res[1, "avg_molecular_weight"], 169.1811)
}

test_parseXmlFolder <- function(){
    xf <- system.file("extdata/", package="xcmsExtensions")
    res <- xcmsExtensions:::.simpleParseHmdbXmlFolder(xf)
    checkEquals(nrow(res), 5)
    checkEquals(res$accession, c("HMDB00001", "HMDB00002", "HMDB00005", "HMDB01346", "HMDB01347"))
}


