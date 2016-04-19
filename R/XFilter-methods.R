####============================================================
##  Implementations for CompoundidFilter
##
##  o column
##  o where
####------------------------------------------------------------
setMethod("column", signature(object="CompoundidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=CompoundidFilter, db=missing not implemented!")
              return(NULL)
          })
setMethod("column", signature(object="CompoundidFilter", db="SimpleCompoundDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="CompoundidFilter", db="SimpleCompoundDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, "accession", with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="CompoundidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              ## without a database we're just calling the where of BasicFilter
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="CompoundidFilter", db="SimpleCompoundDb",
                             with.tables="missing"),
          function(object, db, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="CompoundidFilter", db="SimpleCompoundDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- where(object) ## that returns the result from the BasicFilter...
              return(paste(column(object, db, with.tables=with.tables), suff))
          })

## Remove these for Bioconductor 3.3
## >>
setReplaceMethod("condition", "CompoundidFilter", function(x, value){
    if(x@.valueIsCharacter){
        allowed <- c("=", "!=", "in", "not in", "like")
        if(!any(allowed == value)){
            stop("Only ", paste(allowed, collapse=", "), " are allowed if the value from",
                 " the filter is of type character!")
        }
        if(value == "=" & length(x@value) > 1)
            value <- "in"
        if(value == "!=" & length(x@value) > 1)
            value <- "not in"
        if(value == "in" & length(x@value) == 1)
            value <- "="
        if(value == "not in" & length(x@value) == 1)
            value <- "!="
    }else{
        allowed <- c("=", ">", "<", ">=", "<=")
        if(!any(allowed == value)){
            stop("Only ", paste(allowed, collapse=", "), " are allowed if the value from",
                 " the filter is numeric!")
        }
    }
    x@condition <- value
    validObject(x)
    return(x)
})
setMethod("value", signature(x="CompoundidFilter", db="missing"),
          function(x, db, ...){
              return(x@value)
          })
setReplaceMethod("value", "CompoundidFilter", function(x, value){
    if(is.numeric(value)){
        x@.valueIsCharacter <- FALSE
    }else{
        x@.valueIsCharacter <- TRUE
    }
    x@value <- as.character(value)
    ## Checking if condition matches the value.
    if(length(value) > 1){
        if(x@condition == "=")
            x@condition <- "in"
        if(x@condition == "!=")
            x@condition <- "not in"
    }else{
        if(x@condition == "in")
            x@condition <- "="
        if(x@condition == "not in")
            x@condition <- "!="
    }
    ## Test validity
    validObject(x)
    return(x)
})
## <<
