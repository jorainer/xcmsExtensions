####============================================================
##  Methods for the SimpleCompoundDb
####------------------------------------------------------------

####============================================================
##  show
##
##  show method for SimpleCompoundDb
####------------------------------------------------------------
setMethod("show", "SimpleCompoundDb", function(object){
    cat(class(object), " object:\n", sep="")
    metad <- dbGetQuery(dbconn(object), "select * from metadata;")
    cat("| metadata: \n")
    for(i in 1:nrow(metad)){
        cat("| o", metad[i, "name"], ": ", metad[i, "value"],"\n")
    }
    cat("| Data summary:\n")
    sources <- dbGetQuery(dbconn(object), "select distinct source from compound_basic;")[, "source"]
    for(i in 1:length(sources)){
        cat("| o ", sources[i], ": ",
            dbGetQuery(dbconn(object), paste0("select count(*) from compound_basic where source=",
                                              sQuote(sources[i])))[1, 1], " compounds.\n")
    }
})


####============================================================
##  dbconn
##
##  getter method for the database connection.
####------------------------------------------------------------
setMethod("dbconn", "SimpleCompoundDb", function(x){
    return(x@con)
})

## ####============================================================
## ##  keys
## ##
## ##  returns all ATC codes stored in the database.
## ####------------------------------------------------------------
## setMethod("keys", "AtcDb", function(x, level, ...){
##     wQuery <- .levelConditionQuery(level)
##     if(length(wQuery) > 0)
##         wQuery <- paste0("where ", wQuery)
##     Q <- paste0("select key from atc ", wQuery, " order by key;")
##     res <- dbGetQuery(dbconn(x), Q)
##     return(res[, "key"])
## })

####============================================================
##  columns
##
##  returns all columns available  in the database.
####------------------------------------------------------------
setMethod("columns", "SimpleCompoundDb", function(x){
    Tabs <- listTables(x)
    return(unique(sort(unlist(Tabs[names(Tabs) != "metadata"], use.names=FALSE))))
})

####============================================================
##  listTables
##
##  lists the database tables
####------------------------------------------------------------
setMethod("listTables", "SimpleCompoundDb", function(x, ...){
    if(length(x@tables) > 0)
        return(x@tables)
    ## Else: query:
    return(.doListTables)
})
.doListTables <- function(x){
    con <- x@con
    Tables <- dbListTables(con)
    theTables <- vector("list", length(Tables))
    for(i in 1:length(Tables)){
        theTables[[i]] <- colnames(dbGetQuery(con,
                                              paste0("select * from ",
                                                     Tables[i],
                                                     " limit 1;")))
    }
    names(theTables) <- Tables
    return(theTables)
}

####============================================================
##  as.data.frame
##
##
####------------------------------------------------------------
setMethod("as.data.frame", "SimpleCompoundDb", function(x, ...){
    columns <- cleanColumns(x, columns(x))
    Res <- getWhat(x, columns=columns, start.table="compound_basic", join="left join",
                   order.by="accession")

    ## cols <- unlist(.checkColumns(x), use.names=FALSE)
    ## Res <- dbGetQuery(dbconn(x),
    ##                   paste0("select ",
    ##                          paste0(cols, collapse=","),
    ##                          " from compound_basic order by accession;"))
    ## ## setting rownames
    ## ## rownames(Res) <- Res$key
    return(Res)
})


## ####============================================================
## ##  .checkColumns
## ##
## ##  check the columns and define the tables from which to return values
## ##  based on them. Returns a list (!) of column names, with the names of
## ##  the list being the table names. Note that the column names have also
## ##  been "sanitized" for the database table (i.e. <table>.<column>).
## ####------------------------------------------------------------
## .checkColumns <- function(x, columns=NULL){
##     haveTables <- listTables(x)
##     tableDf <- cbind(table=rep(names(haveTables), unlist(lapply(haveTables, length))),
##                      column=unlist(haveTables, use.names=FALSE))
##     ## Exclude the metadata table. Eventually we would need to re-order the df...
##     tableDf <- tableDf[tableDf[, 1] != "metadata", ]
##     allowedCols <- unique(tableDf[, 2])
##     if(is.null(columns))
##         columns <- allowedCols
##     notAllowed <- columns[!(columns %in% allowedCols)]
##     columns <- columns[columns %in% allowedCols]
##     if(length(columns) == 0)
##         stop("None of the submitted columns are present in the database!")
##     if(length(notAllowed) > 0){
##         warning("Columns ", paste0(notAllowed, ", "), " not present in ",
##                 "the database. These have been removed.")
##     }
##     ## Use match to select the first occurence in case we have two columns with
##     ## the same name.
##     idx <- match(columns, tableDf[, 2])
##     tableDf <- tableDf[idx, , drop=FALSE]
##     columns <- paste0(tableDf[, 1], ".", tableDf[, 2])
##     return(split(columns, f=tableDf[, 1]))
## }



####============================================================
##  mzmatch
##
##  Crazy cool mzmatch method that looks compounds up in the
##  HMDB database!
####------------------------------------------------------------
setMethod("mzmatch", signature(x="numeric", mz="SimpleCompoundDb"),
          function(x, mz, mzdev=0, ppm=10, column="monoisotopic_molecular_weight",
                   ionAdduct=NULL){
              column <- match.arg(column, c("avg_molecular_weight",
                                            "monoisotopic_molecular_weight"))
              ## Check if we would have the ion adduct available.
              adds <- ionAdduct
              if(!is.null(adds)){
                  notPresent <- !(adds %in% supportedIonAdducts())
                  if(all(notPresent))
                      stop("None of the provided ion adduct names are supported!",
                           " Use the 'supportedIonAdducts()' method to list all supported names.")
                  if(any(notPresent)){
                      warning("Adducts with name ", paste(adds[notPresent], collapes=", "), " are not supported",
                              " and have thus been removed.")
                      ionAdduct <- ionAdduct[!notPresent]
                  }
              }
              ## return(.mzmatchCompoundDbSQL(x=x, mz=mz, mzdev=mzdev, ppm=ppm,
              ##                              column=column))
              return(.mzmatchMassPpmCompoundDbSql(x=x, mz=mz, mzdev=mzdev, ppm=ppm,
                                                  column=column, ionAdduct=ionAdduct))
          })
## Here we support calculation and comparison of all masses for M/Z being all possible
## ion adducts for a compound.
## We are first calculating the minmz and maxmz given the mzdev and ppm, then we get
## all possible masses for these (considering ionName) and feeding that to the SQL query.
## So, here the ppm is on the M/Z, not on the mass.
.mzmatchMassCompoundDbSql <- function(x, mz, mzdev=0, ppm=10,
                                      column="monoisotopic_molecular_weight",
                                      ionAdduct=NULL){
    fixi <- 1e-9
    mzd <- x * ppm/1000000 + mzdev + fixi
    minmasses <- adductmz2mass(x-mzd, ionAdduct=ionAdduct)
    maxmasses <- adductmz2mass(x+mzd, ionAdduct=ionAdduct)
    con <- dbconn(mz)
    df <- data.frame(mzidx=rep(1:length(mzd), length(minmasses)),
                     adduct=rep(names(minmasses), each=length(x)),
                     minmass=unlist(minmasses, use.names=FALSE),
                     maxmass=unlist(maxmasses, use.names=FALSE),
                     stringsAsFactors=FALSE
                     )
    res <- lapply(split(df, f=1:nrow(df)), function(z, col){
        quer <- paste0("select accession, ", col, " from compound_basic where ", col,
                       " between ", z$minmass, " and ", z$maxmass, ";")
        tmp <- dbGetQuery(con, quer)
        if(nrow(tmp) == 0){
            return(data.frame(idx=NA, deltaMz=NA, mzidx=z$mzidx, adduct=z$adduct,
                              stringsAsFactors=FALSE))
        }
        ## Calculate distance.
        return(data.frame(idx=tmp$accession,
                          deltaMz=abs(tmp[, col] - mean(as.numeric(z[1, c(3, 4)]))),
                          mzidx=rep(z$mzidx, nrow(tmp)),
                          adduct=rep(z$adduct, nrow(tmp)),
                          stringsAsFactors=FALSE))
    },
    col=column
    )
    ## rbind that and split it again by mzidx; dropping the names later.
    res <- do.call(rbind, c(res, make.row.names=FALSE))
    rownames(res) <- NULL
    ## Split it by index of input.
    res <- split(res[, c(1, 2, 4)], f=res$mzidx)
    names(res) <- NULL
    return(res)
}
## Same as .mzmatchMassCompoundDbSql, but FIRST the masses of all adducts is calculated
## and then the query is performed using the mass +- its PPM, so PPM is calculated on the
## MASS.
.mzmatchMassPpmCompoundDbSql <- function(x, mz, mzdev=0, ppm=10,
                                      column="monoisotopic_molecular_weight",
                                      ionAdduct=NULL){
    fixi <- 1e-9
    ppm <- ppm/1000000
    ## Calculate the ion Adducts.
    masses <- adductmz2mass(x, ionAdduct=ionAdduct)
    con <- dbconn(mz)
    df <- data.frame(mzidx=rep(1:length(x), length(masses)),
                     adduct=rep(names(masses), each=length(x)),
                     mass=unlist(masses, use.names=FALSE),
                     stringsAsFactors=FALSE)
    ## Split it by row and perform the query.
    res <- lapply(split(df, f=1:nrow(df)), function(z, col){
        addP <- z$mass * ppm + mzdev + fixi
        quer <- paste0("select accession, ", col, " from compound_basic where ", col,
                       " between ", z$mass - addP, " and ", z$mass + addP, ";")
        tmp <- dbGetQuery(con, quer)
        if(nrow(tmp) == 0){
            return(data.frame(idx=NA, deltaMz=NA, mzidx=z$mzidx, adduct=z$adduct,
                              stringsAsFactors=FALSE))
        }
        ## Calculate distance.
        return(data.frame(idx=tmp$accession,
                          deltaMz=abs(tmp[, col] - as.numeric(z[1, 3])),
                          mzidx=rep(z$mzidx, nrow(tmp)),
                          adduct=rep(z$adduct, nrow(tmp)),
                          stringsAsFactors=FALSE))
    }, col=column)
    ## rbind that and split it again by mzidx; dropping the names later.
    res <- do.call(rbind, c(res, make.row.names=FALSE))
    rownames(res) <- NULL
    ## Split it by index of input.
    res <- split(res[, c(1, 2, 4)], f=res$mzidx)
    names(res) <- NULL
    return(res)
}

## This method aims to use SQL commands to find the match.
.mzmatchCompoundDbSQL <- function(x, mz, mzdev=0, ppm=10, column="monoisotopic_molecular_weight"){
    ## hm, could use SQL between: where xx between LB and UB; (between is >= and <=).
    ## Create a list of query strings and apply them.

    ## OK, what do I have to do:
    ## calculate the min and the max mz based on the mzdev and ppm. convert that to the mass assuming
    ## mz is an ion adduct.
    ## Use this as input in the search.

    fixi <- 1e-9
    con <- dbconn(mz)
    res <- lapply(x, function(z, mzd, ppm, col){
        addP <- (z/1000000) * ppm + mzdev + fixi
        query <- paste0("select accession, ", col, " from compound_basic where ", col,
                        " between ", z - addP ," and ", z + addP, ";")
        tmp <- dbGetQuery(con, query)
        if(nrow(tmp) == 0)
            return(data.frame(idx=NA, deltaMz=NA))
        ## Calculate the distance...
        return(data.frame(idx=tmp$accession,
                          deltaMz=abs(tmp[, col] - z),
                          stringsAsFactors=FALSE))
    }, mzd=mzdev, ppm=ppm, col=column)
    return(res)
}

## This method first fetches all the data and performs the search
## in R. It is however slower than the one above!
.mzmatchCompoundDbPlain <- function(x, mz, mzdev=0, ppm=10, column="avg_molecular_weight"){
    ## Just get all of the data.
    mz <- as.data.frame(mz)
    res <- mzmatch(x, mz=mz[, column], mzdev=mzdev, ppm=ppm)
    ## Need to get the accession for each.
    return(lapply(res, function(z){
        if(any(is.na(z[, 1])))
            return(z)
        return(data.frame(idx=mz[z[, 1], "accession"],
                          deltaMz=z[, 2], stringsAsFactors=FALSE))
    }))
}

####============================================================
##  mzmatch
##
##  mzmatch supporting a matrix with mzmin and mzmax as input. Masses in the database
##  are then matched against this, i.e. matched if the mass is within the range (+-ppm)
####------------------------------------------------------------
setMethod("mzmatch", signature(x="matrix", mz="SimpleCompoundDb"),
          function(x, mz, mzdev=0, ppm=10, column="monoisotopic_molecular_weight",
                   ionAdduct=NULL){
              column <- match.arg(column, c(
                                            "monoisotopic_molecular_weight",
                                            "avg_molecular_weight"))
              ## Check input matrix.
              if(ncol(x) > 2){
                  ## Require mzmin and mzmax.
                  if(!all(c("mzmin", "mzmax") %in% colnames(x)))
                      stop("If x has more than two columns, it has to have columns named ",
                           "'mzmin' and 'mzmax'!")
                  x <- x[, c("mzmin", "mzmax")]
              }
              ## Check that min is always <= max
              if(any((x[, 1] - x[, 2]) > 0))
                  stop("The first column ('mzmin') is expected to have smaller",
                       " values than the second ('mzmax').")
              ## Check if we would have the ion adduct available.
              adds <- ionAdduct
              if(!is.null(adds)){
                  notPresent <- !(adds %in% supportedIonAdducts())
                  if(all(notPresent))
                      stop("None of the provided ion adduct names are supported!",
                           " Use the 'supportedIonAdducts()' method to list all supported names.")
                  if(any(notPresent)){
                      warning("Adducts with name ", paste(adds[notPresent], collapes=", "), " are not supported",
                              " and have thus been removed.")
                      ionAdduct <- ionAdduct[!notPresent]
                  }
              }
              ## return(.mzmatchCompoundDbSQL(x=x, mz=mz, mzdev=mzdev, ppm=ppm,
              ##                              column=column))
              return(.mzmatchMassPpmMatrixCompoundDbSql(x=x, mz=mz, mzdev=mzdev, ppm=ppm,
                                                  column=column, ionAdduct=ionAdduct))
          })
## That's essentially the same than the mzmatchMassPpmCompoundDbSql function,
## just assuming x is a matrix with mzmin and mzmax. These are first converted
## to mass, and the ppm is then applied to the mass search.
.mzmatchMassPpmMatrixCompoundDbSql <- function(x, mz, mzdev=0, ppm=10,
                                      column="monoisotopic_molecular_weight",
                                      ionAdduct=NULL){
    ## settings:
    ## cat("mzdev: ", mzdev, "\nppm: ", ppm, "\ncolumn: ", column, "\n")
    fixi <- 1e-9
    ppm <- ppm/1000000
    if(ncol(x) != 2)
        stop("'x' is supposed to be a two-column matrix!")
    minmasses <- adductmz2mass(x[, 1], ionAdduct=ionAdduct)
    maxmasses <- adductmz2mass(x[, 2], ionAdduct=ionAdduct)
    con <- dbconn(mz)
    df <- data.frame(mzidx=rep(1:nrow(x), length(minmasses)),
                     adduct=rep(names(minmasses), each=nrow(x)),
                     minmass=unlist(minmasses, use.names=FALSE),
                     maxmass=unlist(maxmasses, use.names=FALSE),
                     stringsAsFactors=FALSE
                     )
    res <- lapply(split(df, f=1:nrow(df)), function(z, col){
        quer <- paste0("select accession, ", col, " from compound_basic where ", col,
                       " between ", z$minmass - (z$minmass * ppm + mzdev + fixi),
                       " and ", z$maxmass + (z$maxmass * ppm + mzdev + fixi), ";")
        tmp <- dbGetQuery(con, quer)
        if(nrow(tmp) == 0){
            return(data.frame(idx=NA, deltaMz=NA, mzidx=z$mzidx, adduct=z$adduct,
                              stringsAsFactors=FALSE))
        }
        ## Calculate distance: minimal distance of the mass to the peak defined by min and max mass,
        ## i.e. the distance is the min of the distance mass to min mass and mass to max mass.
        ## ## massmat <- matrix(rep(as.numeric(z[1, c(3, 4)]), each=nrow(tmp)), ncol=2)
        ## massmat <- matrix(as.numeric(z[1, c(3, 4)]), nrow=nrow(tmp), ncol=2, byrow=TRUE)
        ## massDist <- abs(rowMin(tmp[, col] - massmat))
        ## Note: better calculate the distance to the middle position of the range.
        massDist <- abs(tmp[, col] - mean(as.numeric(z[1, c(3, 4)])))
        return(data.frame(idx=tmp$accession,
                          deltaMz=massDist,
                          mzidx=rep(z$mzidx, nrow(tmp)),
                          adduct=rep(z$adduct, nrow(tmp)),
                          stringsAsFactors=FALSE))
    },
    col=column
    )
    ## rbind that and split it again by mzidx; dropping the names later.
    res <- do.call(rbind, c(res, make.row.names=FALSE))
    rownames(res) <- NULL
    ## Split it by index of input.
    res <- split(res[, c(1, 2, 4)], f=res$mzidx)
    names(res) <- NULL
    return(res)
}
## ## Alternative version
## .mzmatchMassPpmMatrixCompoundDbSql2 <- function(x, mz, mzdev=0, ppm=10,
##                                                 column="monoisotopic_molecular_weight",
##                                                 ionAdduct=NULL){
##     ## Alternative version of the function (eventually faster?):
##     ## o Calculate the mean mz, and an mzdev ()
##     mzMean <- rowMeans(x)
##     ## o The difference to the minz and maxz
##     mzDevs <- abs(mzMean - x[, 1]) + mzdev
##     fixi <- 1e-9
##     ppm <- ppm/1000000
##     if(ncol(x) != 2)
##         stop("'x' is supposed to be a two-column matrix!")
##     ## o Just calculate the mass once.
##     masses <- adductmz2mass(mzMean, ionAdduct=ionAdduct)
##     ## masses is now a list, names being adducts with elements being the calculated
##     ## mass for each of the input mean M/Z
##     con <- dbconn(mz)
##     df <- data.frame(mzidx=rep(1:nrow(x), length(minmasses)),
##                      adduct=rep(names(masses), each=nrow(x)),
##                      meanmass=unlist(masses, use.names=FALSE),
##                      mzdev=rep(mzDevs, length(minmasses)),
##                      stringsAsFactors=FALSE
##                      )
##     res <- lapply(split(df, f=1:nrow(df)), function(z, col){
##         mzd <- z$meanmass * ppm + z$mzdev + fixi
##         quer <- paste0("select accession, ", col, " from compound_basic where ", col,
##                        " between ", z$meamass - mzd,
##                        " and ", z$meanmass + mzd, ";")
##         tmp <- dbGetQuery(con, quer)
##         if(nrow(tmp) == 0){
##             return(data.frame(idx=NA, deltaMz=NA, mzidx=z$mzidx, adduct=z$adduct,
##                               stringsAsFactors=FALSE))
##         }
##         ## Calculate distance.
##         stop("Distance calculation has to be re-checked!")
##         return(data.frame(idx=tmp$accession,
##                           deltaMz=min(abs(tmp[, col] - c(as.numeric(z[1, 3]) - as.numeric(z[1, 4]),
##                                                          as.numeric(z[1, 3]) + as.numeric(z[1, 4])))),
##                           mzidx=rep(z$mzidx, nrow(tmp)),
##                           adduct=rep(z$adduct, nrow(tmp)),
##                           stringsAsFactors=FALSE))
##     },
##     col=column
##     )
##     ## rbind that and split it again by mzidx; dropping the names later.
##     res <- do.call(rbind, c(res, make.row.names=FALSE))
##     rownames(res) <- NULL
##     ## Split it by index of input.
##     res <- split(res[, c(1, 2, 4)], f=res$mzidx)
##     names(res) <- NULL
##     return(res)
## }


####============================================================
##  compounds
##
##  The main accessor method to retrieve compounds from the database.
####------------------------------------------------------------
setMethod("compounds", "SimpleCompoundDb",
          function(x, columns, filter=list(), order.by="accession", ...){
              startTable <- "compound_basic"
              if(missing(columns))
                  columns <- listTables(x)[[startTable]]
              columns <- cleanColumns(x, columns)
              ## Get the data.
              res <- getWhat(x, columns=columns, filter=filter, join="left join",
                             start.table=startTable, order.by=order.by, ...)
              return(res)
          })


####------------------------------------------------------------
####============================================================
##  Internal database stuff.
##
##  Mostly taken from mirhostgenes and ensembldb packages.
####============================================================

####============================================================
##  getWhat.
##
## That's the main entrance point for all stuff querying the database.
## just to add another layer; basically just calls buildQuery and executes the query
## return.all.columns: returns also columns added because of the filters.
setMethod("getWhat", "SimpleCompoundDb",
          function(x, columns, filter=list(), order.by="", order.type="asc",
                   skip.order.check=FALSE, join="join", start.table,
                   return.all.columns=TRUE){
              if(missing(columns)){
                  tabs <- listTables(x)
                  columns <- unlist(tabs[names(tabs) != "metadata"], use.names=FALSE)
              }
              Q <- buildQuery(x=x, columns=columns, filter=filter, order.by=order.by,
                              order.type=order.type,
                              skip.order.check=skip.order.check, join=join,
                              start.table=start.table, return.all.columns=return.all.columns)
              return(dbGetQuery(x@con, Q))
})

####============================================================
##  cleanFilter
##
##  loops through the filters and checks if the filters are supported
##  in the present database, i.e. if their column is present in the
##  database. The method removes all filters that are not supported
##  with a respective warning message.
##  x: the SimpleCompoundDb object.
##  filter: the list of filter objects.
####------------------------------------------------------------
setMethod("cleanFilter", "SimpleCompoundDb", function(x, filter){
    if(missing(filter))
        filter <- list()
    if(is(filter, "BasicFilter"))
        filter <- list(filter)
    if(is(filter, "list")){
        tabs <- listTables(x)
        if(length(filter) > 0){
            res <- lapply(filter, function(z){
                tmp <- try(column(z, db=x), TRUE)
                if(inherits(tmp, "try-error")){
                    warning("Removed filter ", class(z)[1], " as it is not supported",
                            " by the ", class(x)[1], " database.")
                    return(NULL)
                }else{
                    return(z)
                }
            })
            res <- res[lengths(res) > 0]
        }else{
            return(list())
        }
    }else{
        stop("Argument 'filter' has to be either a single filter object inheriting",
             " from 'BasicFilter' or a list of such objects!")
    }
})


####============================================================
##  buildQuery
##  x: SimpleCompoundDb
##  filter: list of BasicFilter objects. Should already be checked by the cleanFilter function.
##  order.by: if we should add an order by to the query. The provided column will also be checked
##            if it is available in the database.
##  order.type: asc or desc.
##  skip.order.check: whether the checking of the "order.by" should be omitted (i.e. if we're submitting
##                    a more complex order.by and not just the column name.
##  join: which join should be performed. join (default), left join, left outer join.
##  start.table: from which table should the join start?
##  return.all.columns: if all columns should be returned (also filter columns etc), or just
##                      the columns passed with argument columns.
##  Returns: a character string with the query.
####------------------------------------------------------------
setMethod("buildQuery", "SimpleCompoundDb",
          function(x, columns, filter=list(), order.by="", order.type="asc",
                   skip.order.check=FALSE, join="join",
                   start.table, return.all.columns=TRUE){
              resultcolumns <- columns    ## just to remember what we really want to give back
              join <- match.arg(join, c("join", "left join", "left outer join"))
              ## 1) get all column names from the filters also removing the prefix.
              filter <- cleanFilter(x, filter)
              if(length(filter) > 0){
                  ## Add the columns needed for the filter
                  filtercolumns <- unlist(lapply(filter, column, x))
                  ## Remove the prefix (column name for these)
                  filtercolumns <- .removePrefix(filtercolumns)
                  columns <- unique(c(columns, filtercolumns))
              }
              ## 2) Get all column names for the order.by, check if they exist in
              ##    the database:
              if(order.by!=""){
                  ## if we have skip.order.check set we use the order.by as is.
                  if(!skip.order.check){
                      order.by <- unlist(strsplit(order.by, split=",", fixed=TRUE))
                      order.by <- gsub(order.by, pattern=" ", replacement="", fixed=TRUE)
                      ## Check if order.by are valid columns.
                      order.by <- cleanColumns(x, order.by)
                      if(length(order.by)==0){
                          order.by <- ""
                      }else{
                          columns <- unique(c(columns, order.by))
                          order.by <- paste(unlist(prefixColumns(x, columns=order.by,
                                                                 with.tables=names(prefixColumns(x, columns)))
                                                   , use.names=FALSE), collapse=", ")
                      }
                  }
              }else{
                  order.by <- ""
              }
              ## Note: order by is now a vector!!!
              ## columns are now all columns that we want to fetch or that we need to
              ## filter or to sort.
              ## 3) check which tables we need for all of these columns:
              need.tables <- names(prefixColumns(x, columns))
              ##
              ## Now we can begin to build the query parts!
              ## a) the query part that joins all required tables.
              joinquery <- buildJoinQuery(x, columns=columns, join=join,
                                          start.table=start.table)
              ## b) the filter part of the query
              filterquery <- buildFilterQuery(x, filter=filter, with.tables=need.tables)
              ## c) the order part of the query
              if(order.by!=""){
                  orderquery <- paste(" order by", order.by, order.type)
              }else{
                  orderquery <- ""
              }
              ## And finally build the final query
              if(return.all.columns){
                  resultcolumns <- columns
              }
              finalquery <- paste0("select distinct ",
                                   paste(unlist(prefixColumns(x,
                                                              resultcolumns,
                                                              with.tables=need.tables),
                                                use.names=FALSE), collapse=","),
                                   " from ",
                                   joinquery,
                                   filterquery,
                                   orderquery
                                   )
              return(finalquery)
          })

####============================================================
##  buildJoinQuery
##
####------------------------------------------------------------
setMethod("buildJoinQuery", "SimpleCompoundDb", function(x, columns, join="join", start.table){
    ## Get the table names and prefixed column names.
    tabs <- names(prefixColumns(x, columns))
    ## Add tables that might be needed to join the specified tables.
    tabs <- addRequiredJoinTables(x, tabs)
    ## Sort the tables.
    tabs <- sortTablesByDegree(x, tabs)
    if(!missing(start.table)){
        if(any(tabs == start.table)){
            tabs <- c(start.table, tabs[tabs != start.table])
        }
    }
    query <- tabs[1]
    previousTable <- tabs[1]
    remainingTable <- tabs[-1]
    ## Now loop through the tables and add them sequentially.
    while(length(previousTable)!=length(tabs)){
        ## Repeat until I don't have joined all tables!
        ## check which joins are available for the last table(s)
        ## Basically we're extracting those rows from the .JOINS matrix that have
        ## the previousTable in the first or second column.
        previousIdx <- which(apply(.JOINS[ , 1:2, drop=FALSE ], MARGIN=1,
                                    function(z){ any(z %in% previousTable) }))
        if(length(previousIdx) == 0)
            stop("Argh, something weird happened...")
        ## Now check if in the remaining tables there is one that could be joined to this one.
        tmp <- .JOINS[ previousIdx, , drop=FALSE ]
        remainingIdx <- which(apply(tmp[ , 1:2, drop=FALSE ], MARGIN=1,
                                     function(z){ any(z %in% remainingTable) }))
        ## add this table to the previousTable vector
        previousTable <- unique(c(previousTable,
                                  tmp[remainingIdx[1], 1:2][tmp[remainingIdx[1], 1:2]!=previousTable[length(previousTable)]]))
        remainingTable <- remainingTable[remainingTable!=previousTable[length(previousTable)]]
        query <- paste(query, join, previousTable[length(previousTable)],
                       tmp[remainingIdx[1], 3])
    }
    return(query)
})

####============================================================
##  addRequiredTables
##
##  Add additional tables in case the submitted ones can not be joined
##  directly.
####------------------------------------------------------------
setMethod("addRequiredJoinTables", "SimpleCompoundDb", function(x, tables){
    ## Presently there is no need at all to do any joining... In case: check dbhelpers.R
    ## in mirhostgenes package, addRequiredTables function.
    return(tables)
})

.JOINS <- matrix(ncol=3, nrow=0)

####============================================================
##  buildFilterQuery
##
##  with.tables is optional to force usage of the specified tables.
####------------------------------------------------------------
setMethod("buildFilterQuery", "SimpleCompoundDb", function(x, filter, with.tables){
    filter <- cleanFilter(x, filter)
    if(missing(with.tables)){
        with.tables <- names(listTables(x))
    }
    if(length(filter) == 0)
        return("")
    query <- paste0(" where ",
                    paste(unlist(lapply(filter, where, x,
                                        with.tables=with.tables)),
                          collapse=" and "))
    return(query)
})

####============================================================
##  sortTablesByDegree
##
##  Sort the given table names by their degree, i.e. number of other
##  tables with which they are connected (via primary/foreign keys).
##  Returns a re-ordered character vector of table names.
####------------------------------------------------------------
setMethod("sortTablesByDegree", "SimpleCompoundDb",
          function(x, tables=names(listTables(x))){
              ## Well, at present we just have two tables anyway...
              tableOrder <- c(compound_basic=1, metadata=99)
              tables <- cleanTables(x, tables)
              return(tables[order(tableOrder[tables])])
          })

####============================================================
##  cleanColumns
##
##  Checks if the provided columns are available in the database
##  and removes those that are not in the database.
##  The method returns all columns that are in the database, removes all that
##  are not present (with a warning). If none are present, it returns NULL.
setMethod("cleanColumns", "SimpleCompoundDb", function(x, columns, excludes="metadata"){
    if(missing(columns))
        return(NULL)
    tabs <- listTables(x)
    availableCols <- unlist(tabs[!(names(tabs) %in% excludes)], use.names=FALSE)
    areOK <- columns %in% availableCols
    notOK <- columns[!areOK]
    if(length(notOK) > 0)
        warning("Columns ", paste(sQuote(notOK), collapse=", "), " are not available",
                " in the database and have thus been removed.")
    columns <- columns[areOK]
    if(length(columns) == 0){
        warning("None of the provided columns are known to the database!")
        return(NULL)
    }
    return(columns)
})

####============================================================
##  cleanTables
##
##  Checks if the provided table names are available in the database and removes
##  those that aren't
####------------------------------------------------------------
setMethod("cleanTables", "SimpleCompoundDb", function(x, tables){
    if(missing(tables))
        return(NULL)
    gotTabs <- names(listTables(x))
    areOK <- tables %in% gotTabs
    notOK <- tables[!areOK]
    if(length(notOK) > 0)
        warning("Tables ", paste(sQuote(notOK), collapse=", "), " are not available",
                " in the database and have thus been removed.")
    tables <- tables[areOK]
    if(length(tables) == 0){
        warning("None of the provided tables are known to the database!")
        return(NULL)
    }
    return(tables)
})


####
## Note: that's the central function that checks which tables are needed for the
## least expensive join!!! The names of the tables should then also be submitted
## to any other method that calls .prefixColumns (e.g. where of the Filter classes)
##
## this function checks:
## a) for multi-table attributes, selects the table with the highest degree (i.e.
##    the table connected to most other tables).
## b) pre-pend (inverse of append ;)) the table name to the attribute name.
## returns a list, names being the tables and the values being the attributes
## named: <table name>.<attribute name>
## clean: whether a cleanColumns should be called on the submitted attributes.
## with.tables: force the prefix to be specifically on the submitted tables.
##
## Return values:
## If none of the provided columns is present in the database, the
## function (given that clean=TRUE) returns NULL otherwise it returns
## a list, list names being the table names and elements the (prefixed) column names.
setMethod("prefixColumns", "SimpleCompoundDb", function(x, columns, clean=TRUE, with.tables){
    if(missing(columns))
        stop("columns is empty! No columns provided!")
    ## first get to the tables that contain these columns
    Tab <- listTables(x)   ## returns the tables by degree!
    if(!missing(with.tables)){
        with.tables <- with.tables[with.tables %in% names(Tab)]
        if(length(with.tables) > 0)
            Tab <- Tab[with.tables]
        if(length(Tab)==0)
            stop("None of the tables submitted with with.tables is present in the database!")
    }
    if(clean)
        columns <- cleanColumns(x, columns)
    if(length(columns)==0){
        return(NULL)
    }
    ## Group the columns by table.
    columns.bytable <- lapply(Tab, function(z){
        return(z[z %in% columns])
    })
    ## Kick out empty tables...
    columns.bytable <- columns.bytable[lengths(columns.bytable) > 0]
    if(length(columns.bytable)==0)
        stop("No columns available!")
    have.columns <- NULL
    ## Order the tables by the number of columns they have.
    columns.bytable <- columns.bytable[order(lengths(columns.bytable), decreasing=TRUE)]
    ## Remove redundant columns.
    ## Loop throught the columns by table and sequentially kick out columns
    ## for the current table if they where already
    ## in a previous (more relevant) table. That way we might reduce the number of tables
    ## from which we have to fetch data and thus speed up queries (given that different tables)
    ## have similar content.
    for(i in 1:length(columns.bytable)){
        currentCols <- columns.bytable[[i]]
        keepColumns <- currentCols[!(currentCols %in% have.columns)]   ## keep those
        if(length(keepColumns) > 0){
            have.columns <- c(have.columns, keepColumns)
            ## Add the <table name>.<column name>
            columns.bytable[[i]] <- paste(names(columns.bytable)[i], keepColumns, sep=".")
        }else{
            ## Just put in the empty stuff.
            columns.bytable[[i]] <- keepColumns
        }
    }
    return(columns.bytable[lengths(columns.bytable) > 0])
})

####============================================================
##  .removePrefix
##
##  remove the prefix from the column name (i.e. the table name).
####------------------------------------------------------------
.removePrefix <- function(x, split=".", fixed=TRUE){
    if(is(x, "list")){
        ## e.g. what we get back grom .prefixColumns
        return(lapply(x, function(z){
            tmp <- unlist(strsplit(z, split=split, fixed=fixed), use.names=FALSE)
            return(tmp[length(tmp)])
        }))
    }else{
        tmp <- strsplit(x, split=split, fixed=fixed)
        return(unlist(lapply(tmp, function(z){
            return(z[length(z)])
        })))
    }
}

