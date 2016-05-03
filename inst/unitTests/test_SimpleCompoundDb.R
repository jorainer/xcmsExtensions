####============================================================
##  Testing some of the (basic) functionality of the SimpleCompoundDb.
##
####------------------------------------------------------------
detach("package:xcmsExtensions", unload=TRUE)
library(xcmsExtensions)

test_SimpleCompoundDb <- function(){
    ## Just a silly check
    checkException(SimpleCompoundDb())

    ## Just running the methods to see whether we would get any error.
    tmp <- columns(scDb)
    tmp <- dbconn(scDb)
    tmp <- as.data.frame(scDb)
    tmp <- listTables(scDb)
}

test_mzmatch_db <- function(){
    realMz <- c(169.2870, 169.5650, 346.4605)
    queryMz <- realMz - (floor(realMz) / 1000000) * 10
    comps <- c(300.1898, 298.1508, 491.2000, 169.13481, 169.1348, queryMz)

    Res <- mzmatch(comps, scDb, column="avg_molecular_weight")

    ## Compare that to the mzmatch on integer, integer.
    masses <- compounds(scDb, columns=c("accession", "avg_molecular_weight"))
    Res2 <- mzmatch(comps, masses$avg_molecular_weight)
    ## Get the accessions for those
    Res2 <- lapply(Res2, function(z){
        if(!any(is.na(z[, 1]))){
            return(data.frame(idx=masses[z[, 1], "accession"],
                              deltaMz=z[, 2], stringsAsFactors=FALSE)
                   )
        }else{
            return(data.frame(idx=NA, deltaMz=NA))
        }
    })
    tmp1 <- do.call(rbind, Res)
    tmp2 <- do.call(rbind, Res2)
    tmp2 <- cbind(tmp2, adduct=rep("M", nrow(tmp2)), stringsAsFactors=FALSE)
    rownames(tmp1) <- NULL
    rownames(tmp2) <- NULL
    checkEquals(tmp1, tmp2)

    ## Error checking
    checkException(mzmatch(comps, scDb, ionAdduct="sdfdkjf"))
}

test_mzmatch_matrix <- function(){
    ## Real life example: 5, 7, 9
    mzs <- cbind(mzmed=c(324.1584, 327.1989, 329.2000),
                 mzmin=c(324.1238, 327.1683, 329.1970),
                 mzmax=c(324.1632, 327.2000, 329.2000))
    ## Do the M+H search based on mzmed:
    mzmedRes <- mzmatch(mzs[, "mzmed"], scDb, ppm=10, ionAdduct="M+H")
    mzmedMat <- mzmatch(mzs[, c("mzmin", "mzmax")], scDb, ppm=10, ionAdduct="M+H")
    ## We require that all compounds identified by mzmed are also found in the matrix version
    checkTrue(all(do.call(rbind, mzmedRes)$idx %in% do.call(rbind, mzmedMat)$idx))

    ## Check it
    for(i in 1:3){
        ## Now get the masses for the first...
        cmps <- compounds(scDb, filter=CompoundidFilter(mzmedMat[[i]]$idx),
                          columns="monoisotopic_molecular_weight")
        ## Convert masses
        massMin <- adductmz2mass(mzs[i, "mzmin"], ionAdduct="M+H")[[1]]
        massMax <- adductmz2mass(mzs[i, "mzmax"], ionAdduct="M+H")[[1]]
        ## Check if all the masses of the identified compounds are within that range.
        massMin <- massMin - (massMin * 10/1000000)
        massMax <- massMax + (massMax * 10/1000000)
        checkTrue(all(cmps$monoisotopic_molecular_weight > massMin &
                      cmps$monoisotopic_molecular_weight < massMax))
    }

    ## And now compare that the the matrix,numeric version.
    compmasses <- compounds(scDb, columns=c("accession", "monoisotopic_molecular_weight"))
    for(i in 1:3){
        minmass <- adductmz2mass(mzs[i, "mzmin"], ionAdduct="M+H")[[1]]
        maxmass <- adductmz2mass(mzs[i, "mzmax"], ionAdduct="M+H")[[1]]
        massmat <- matrix(c(minmass, maxmass), nro=1)
        SingleRes <- mzmatch(massmat, mz=compmasses$monoisotopic_molecular_weight, ppm=10)[[1]]
        SingleRes <- data.frame(id=compmasses[SingleRes[, "idx"], "accession"], SingleRes,
                                stringsAsFactors=FALSE)
        dbRes <- mzmedMat[[i]]
        dbRes <- dbRes[order(dbRes$deltaMz), ]
        checkEquals(dbRes$idx, SingleRes$id)
        ## Also the distance?
        checkEquals(dbRes$deltaMz, SingleRes$deltaMz)
    }
}


test_mzmatch_db_new <- function(){
    ## This uses now the ion adducts.
    realMz <- c(169.2870, 169.5650, 346.4605)

    queryMz <- realMz - (floor(realMz) / 1000000) * 10
    comps <- c(300.1898, 298.1508, 491.2000, 169.13481, 169.1348, queryMz)

    ###########
    ## The "front-end" methods.
    Res <- mzmatch(comps, scDb, ppm=10, ionAdduct="M+H")


    ###########
    ## The internal functions.
    Res <- xcmsExtensions:::.mzmatchCompoundDbSQL(comps, scDb)

    ## Test the new one.
    Res2 <- xcmsExtensions:::.mzmatchMassCompoundDbSql(comps, mz=scDb, ionAdduct=NULL)
    Res <- do.call(rbind, Res)
    Res2 <- do.call(rbind, Res2)
    rownames(Res) <- NULL
    rownames(Res2) <- NULL
    checkEquals(Res, Res2[, 1:2])

    ## Test the other new one.
    Res3 <- xcmsExtensions:::.mzmatchMassPpmCompoundDbSql(comps, mz=scDb, ionAdduct=NULL)
    Res3 <- do.call(rbind, Res3)
    rownames(Res3) <- NULL
    checkEquals(Res, Res3[, 1:2])

    ## The full version with ppm on the MZ
    Res4 <- xcmsExtensions:::.mzmatchMassCompoundDbSql(comps, mz=scDb,
                                                       ionAdduct=supportedIonAdducts())
    ## and ppm on the mass
    Res5 <- xcmsExtensions:::.mzmatchMassPpmCompoundDbSql(comps, mz=scDb,
                                                          ionAdduct=supportedIonAdducts())
    tmp1 <- Res4[[1]]
    rownames(tmp1) <- NULL
    tmp2 <- Res5[[1]]
    rownames(tmp2) <- NULL
    checkEquals(tmp1, tmp2)
}

notrun_mzmatch_performance <- function(){
    ## Just testing the performance of x SQL queries against one SQL query
    ## and doing the rest in R...

    realMz <- c(169.2870, 169.5650, 346.4605)
    ## Should get them with a 10 ppm:
    queryMz <- realMz - (floor(realMz) / 1000000) * 10
    sqlRes <- xcmsExtensions:::.mzmatchCompoundDbSQL(queryMz, scDb)
    sqlRes2 <- xcmsExtensions:::.mzmatchCompoundDbSQL(realMz, scDb)
    checkEquals(sqlRes, sqlRes2)

    comps <- c(300.1898, 298.1508, 491.2000, 169.13481, 169.1348, queryMz)

    bigComps <- rep(comps, 1000)
    system.time(
        sqlRes <- xcmsExtensions:::.mzmatchCompoundDbSQL(bigComps, scDb)
    )
    ## Takes 4 secs for 8000 compounds; 7 seconds including the distance calc.

    system.time(
        rRes <- xcmsExtensions:::.mzmatchCompoundDbPlain(bigComps, scDb)
    )
    ## Incredible! takes 7.8 secs!!! with accession retrieval: 8.7
}

notrun_test_as.data.frame <- function(){
    full <- as.data.frame(scDb)
    other <- RSQLite::dbGetQuery(dbconn(scDb), "select * from compound_basic order by accession")
    checkEquals(full, other)
}

test_compounds <- function(){
    cf <- CompoundidFilter(c("HMDB00010", "HMDB00002", "HMDB00011"))
    res <- compounds(scDb, filter=cf)
    checkEquals(res$accession, sort(value(cf)))

    ## Just selected columns
    res <- compounds(scDb, filter=cf, columns=c("name", "inchi"))
    checkEquals(res$accession, sort(value(cf)))
    checkEquals(colnames(res), c("accession", "name", "inchi"))

    ## Optional arguments
    res <- compounds(scDb, filter=cf, columns=c("name", "inchi"), return.all.columns=FALSE)
    checkEquals(colnames(res), c("name", "inchi"))
}


test_cleanColumns <- function(){
    res <- xcmsExtensions:::cleanColumns(scDb, c("accession", "gene_id", "bla"))
    checkEquals(res, "accession")
}

test_prefixColumns <- function(){
    ## with and without clean.
    res <- xcmsExtensions:::prefixColumns(scDb, columns=c("accession", "gene_id", "name"))
    checkEquals(res[[1]], c("compound_basic.accession", "compound_basic.name"))

    res <- xcmsExtensions:::prefixColumns(scDb, columns=c("value", "gene_id", "name"),
                                           clean=FALSE)
    checkEquals(names(res), "metadata")
    checkEquals(res[[1]], c("metadata.name", "metadata.value"))
}

test_cleanTables <- function(){
    res <- xcmsExtensions:::cleanTables(scDb, tables=c("agfkg", "asdfdfd"))
    checkEquals(res, NULL)
    res <- xcmsExtensions:::cleanTables(scDb, tables=c("metadata", "compound_basic"))
    checkEquals(res, c("metadata", "compound_basic"))
}

test_sortTablesByDegree <- function(){
    res <- xcmsExtensions:::sortTablesByDegree(scDb, tables=c("adsds", "metadata", "compound_basic"))
    checkEquals(res, c("compound_basic", "metadata"))
}

test_addRequiredJoinTables <- function(){
    res <- xcmsExtensions:::addRequiredJoinTables(scDb, "metadata")
    checkEquals(res, "metadata")
    res <- xcmsExtensions:::addRequiredJoinTables(scDb, "asfkdf")
}

test_buildFilterQuery <- function(){
    res <- xcmsExtensions:::buildFilterQuery(scDb)
    cf <- CompoundidFilter("adffdf")
    gf <- ensembldb::GeneidFilter("asdasdfd")
    res <- xcmsExtensions:::buildFilterQuery(scDb, filter=cf)
    checkEquals(res, " where compound_basic.accession = 'adffdf'")

    res <- xcmsExtensions:::buildFilterQuery(scDb, filter=list(cf, gf))
    checkEquals(res, " where compound_basic.accession = 'adffdf'")
    res <- xcmsExtensions:::buildFilterQuery(scDb, filter=list(cf, cf))
    checkEquals(res, " where compound_basic.accession = 'adffdf' and compound_basic.accession = 'adffdf'")
}

test_buildJoinQuery <- function(){
    res <- xcmsExtensions:::buildJoinQuery(scDb, "name")
    checkEquals(res, "compound_basic")
    res <- xcmsExtensions:::buildJoinQuery(scDb, c("name", "adff"))
    checkEquals(res, "compound_basic")

    res <- xcmsExtensions:::buildJoinQuery(scDb, c("metadata", "asdds"))
    checkEquals(res, NULL)

}

test_buildQuery <- function(){
    res <- xcmsExtensions:::buildQuery(scDb, columns=c("accession", "name"))
    checkEquals(res, "select distinct compound_basic.accession,compound_basic.name from compound_basic")
    res <- xcmsExtensions:::buildQuery(scDb, columns=c("accession", "name"), order.by="smiles")
    checkEquals(res, "select distinct compound_basic.accession,compound_basic.name,compound_basic.smiles from compound_basic order by compound_basic.smiles asc")
    res <- xcmsExtensions:::buildQuery(scDb, columns=c("accession", "name"), order.by=c("smiles,dfadskfd"))
    checkEquals(res, "select distinct compound_basic.accession,compound_basic.name,compound_basic.smiles from compound_basic order by compound_basic.smiles asc")

    ## And with a filter.
    cf <- CompoundidFilter("abc")
    res <- xcmsExtensions:::buildQuery(scDb, columns=c("accession", "name"), filter=cf)
    checkEquals(res, "select distinct compound_basic.accession,compound_basic.name from compound_basic where compound_basic.accession = 'abc'")
}

test_getWhat <- function(){
    ## Get all the data
    cf <- CompoundidFilter("HMDB00002")
    res <- xcmsExtensions:::getWhat(scDb, filter=cf)
    checkTrue(nrow(res) == 1)
    res <- xcmsExtensions:::getWhat(scDb, filter=cf, columns=c("name", "inchi"))
    checkEquals(colnames(res), c("accession", "name", "inchi"))
    res <- xcmsExtensions:::getWhat(scDb, filter=cf, columns=c("name", "inchi"), return.all.columns=FALSE)
    checkEquals(colnames(res), c("name", "inchi"))
}

test_compounds_MassrangeFilter <- function(){
    mrf <- MassrangeFilter(c(300, 310))
    cmps <- compounds(scDb, filter=mrf, columns=c("accession", "avg_molecular_weight",
                                                  "monoisotopic_molecular_weight"))
    checkTrue(all(cmps$monoisotopic_molecular_weight >= 300 &
                  cmps$monoisotopic_molecular_weight <= 310))
    condition(mrf) <- "()"
    cmps <- compounds(scDb, filter=mrf, columns=c("accession", "avg_molecular_weight",
                                                  "monoisotopic_molecular_weight"))
    checkTrue(all(cmps$monoisotopic_molecular_weight > 300 &
                  cmps$monoisotopic_molecular_weight < 310))
    ## Changing the column to avg
    ## mrf@column <- "avg_molecular_weight"

    ## Combine filters.
    cmps <- compounds(scDb, filter=list(mrf, CompoundidFilter("HMDB60116")),
                      columns=c("accession", "avg_molecular_weight",
                                "monoisotopic_molecular_weight"))
    value(mrf) <- c(304, 310)
    cmps <- compounds(scDb, filter=list(mrf, CompoundidFilter("HMDB60116")),
                      columns=c("accession", "avg_molecular_weight",
                                "monoisotopic_molecular_weight"))
    checkTrue(nrow(cmps)==0)
}
