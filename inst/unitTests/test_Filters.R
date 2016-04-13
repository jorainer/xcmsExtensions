####============================================================
##  Testing the filters...
##
####------------------------------------------------------------
detach("package:xcmsExtensions", unload=TRUE)
library(xcmsExtensions)

test_CompoundidFilter <- function(){
    cf <- CompoundidFilter("abcd")
    checkEquals(value(cf), "abcd")
    condition(cf)
    checkException(condition(cf) <- ">")

    checkEquals(column(cf), NULL)
    checkEquals(column(cf, scDb), "compound_basic.accession")
}

test_checkFilters <- function(){
    ## Check what happens if we provide wrong filters.
    cf <- CompoundidFilter("abcd")
    gf <- ensembldb::GeneidFilter("atereadf")
    res <- xcmsExtensions:::cleanFilter(x=scDb, filter=list(cf, gf))
    checkTrue(length(res) == 1)
}
