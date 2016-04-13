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
