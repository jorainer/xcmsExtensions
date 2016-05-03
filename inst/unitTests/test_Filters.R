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

####============================================================
##  MassrangeFilter
##
test_MassrangeFilter <- function(){
    ## Just testing basic functionality/contrsuctor.
    ## xcmsExtensions:::BasicrangeFilter()
    mrf <- MassrangeFilter()  ## Just to check that we don't get any error.
    checkException(mrf <- MassrangeFilter(value="d"))
    checkException(mrf <- new("MassrangeFilter", value="d"))
    checkException(mrf <- MassrangeFilter(value=c(3, "a")))
    checkException(mrf <- MassrangeFilter(condition="="))
    ## OK, now speaking
    mrf <- MassrangeFilter(value=c(3, 4), condition="(]", column="something")
    checkEquals(mrf@condition, "(]")
    checkEquals(mrf@column, "something")

    ## Check the column.
    checkEquals(column(mrf), "something")
    checkException(column(mrf, scDb))

    ## Now, check the validate thing; what happens if we manually change slot
    ## values.
    mrf@condition <- "="
    checkException(value(mrf))
    checkException(value(mrf) <- c(2, 3))
    mrf <- MassrangeFilter(value=c(3, 4), condition="(]", column="something")

    ## Check column with SimpleCompoundDb; should throw an error.
    checkException(column(mrf, scDb))

    ## Check the where call:
    ## without db
    checkEquals(where(mrf), "something > 3 and something <= 4")
    ## with db, no tables
    checkException(where(mrf, scDb))
    mrf@column <- "monoisotopic_molecular_weight"
    checkEquals(where(mrf, scDb), "compound_basic.monoisotopic_molecular_weight > 3 and compound_basic.monoisotopic_molecular_weight <= 4")
    ## with db and tables
    checkEquals(where(mrf, scDb, "compound_basic"), "compound_basic.monoisotopic_molecular_weight > 3 and compound_basic.monoisotopic_molecular_weight <= 4")
    ## Eventually check also [], (], (), [)
    condition(mrf) <- "[]"
    checkEquals(where(mrf),
                "monoisotopic_molecular_weight >= 3 and monoisotopic_molecular_weight <= 4")
    condition(mrf) <- "(]"
    checkEquals(where(mrf),
                "monoisotopic_molecular_weight > 3 and monoisotopic_molecular_weight <= 4")
    condition(mrf) <- "()"
    checkEquals(where(mrf),
                "monoisotopic_molecular_weight > 3 and monoisotopic_molecular_weight < 4")
    condition(mrf) <- "[)"
    checkEquals(where(mrf),
                "monoisotopic_molecular_weight >= 3 and monoisotopic_molecular_weight < 4")

    ## Check that the default column matches to monoiso...
    mrf2 <- MassrangeFilter()
    checkEquals(column(mrf2), "mass")
    checkEquals(column(mrf2, scDb), "compound_basic.monoisotopic_molecular_weight")
}


notrun_test <- function(){
    ## Starting from scratch; no initialize implemented for Basic
    brf <- new("BasicrangeFilter")
    brf <- new("BasicrangeFilter", value="4")
    brf <- new("BasicrangeFilter", value=c("4", "b"))
    brf <- new("BasicrangeFilter", value=c("4", "8"))
    brf <- new("BasicrangeFilter", value=c("4", "8"), condition="!=")
    brf <- new("BasicrangeFilter", condition="!=")
    brf <- new("BasicrangeFilter", value=c(4, 7), condition="(]")
 }

