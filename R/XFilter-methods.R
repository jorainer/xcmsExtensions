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

####============================================================
##  Implementations for MassrangeFilter
##
##  o column
##  o value: overwriting the BasicFilter method to enable validity check.
##  o condition: overwriting the BasicFilter method to enable validity check.
##  o where
####------------------------------------------------------------
setMethod("column", signature(object="MassrangeFilter", db="missing",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              return(object@column)
          })
## Here we can check if the column exists in the database.
setMethod("column", signature(object="MassrangeFilter", db="SimpleCompoundDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="MassrangeFilter", db="SimpleCompoundDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              ## Default to monoisotopic_molecular_weight.
              if(object@column == "mass")
                  object@column <- "monoisotopic_molecular_weight"
              cn <- columns(db)
              gotCol <- column(object)
              if(!any(gotCol == cn))
                  stop("Column '", gotCol, "' not present in any database table! Please use",
                       " one of the column names returned by the 'columns' method.")
              return(unlist(prefixColumns(db, gotCol, with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("value", signature(x="MassrangeFilter", db="missing"),
          function(x){
              ## Here we ensure again that x is a valid object!
              OK <- getValidity(getClassDef("MassrangeFilter"))(x)
              if(is(OK, "character"))
                  stop(OK)
              return(as.numeric(x@value))
          })
setReplaceMethod("value", "MassrangeFilter", function(x, value){
    ## Setting the value and then checking the validity of the object.
    x@value <- as.character(value)
    OK <- getValidity(getClassDef("MassrangeFilter"))(x)
    if(is(OK, "character"))
        stop(OK)
    return(x)
})
setMethod("condition", "MassrangeFilter",
          function(x){
              ## Here we ensure again that x is a valid object!
              OK <- getValidity(getClassDef("MassrangeFilter"))(x)
              if(is(OK, "character"))
                  stop(OK)
              return(x@condition)
          })
setReplaceMethod("condition", "MassrangeFilter", function(x, value){
    ## Setting the value and then checking the validity of the object.
    x@condition <- as.character(value)
    OK <- getValidity(getClassDef("MassrangeFilter"))(x)
    if(is(OK, "character"))
        stop(OK)
    return(x)
})
setMethod("where", signature(object="MassrangeFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              ## without a database we're just calling the where of BasicFilter
              return(.buildRangeQuery(value(object), column(object), condition(object)))
          }
          )
setMethod("where", signature(object="MassrangeFilter", db="SimpleCompoundDb",
                             with.tables="missing"),
          function(object, db, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="MassrangeFilter", db="SimpleCompoundDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              ## Check if the object is valid!
              OK <- getValidity(getClassDef("MassrangeFilter"))(object)
              if(is(OK, "character"))
                  stop(OK)
              quer <- .buildRangeQuery(value(object),
                                       column(object, db, with.tables=with.tables),
                                       condition(object))
              return(quer)
          })
## Little helper functions...
.buildRangeQuery <- function(value, column, condition="()"){
    condSplit <- unlist(strsplit(condition, split=""))
    quer <- paste(column, .conditionToSQL(condSplit[1]), value[1], "and",
                   column, .conditionToSQL(condSplit[2]), value[2])
    return(quer)
}
.conditionToSQL <- function(x){
    if(x == "(")
        return(">")
    if(x == "[")
        return(">=")
    if(x == ")")
        return("<")
    if(x == "]")
        return("<=")
    return(NULL)
}



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

