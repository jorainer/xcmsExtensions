\name{CompoundidFilter-class}
\Rdversion{1.1}
\docType{class}
%% CompoundidFilter
\alias{CompoundidFilter}
\alias{CompoundidFilter-class}
\alias{column,CompoundidFilter,missing,missing-method}
\alias{column,CompoundidFilter,SimpleCompoundDb,missing-method}
\alias{column,CompoundidFilter,SimpleCompoundDb,character-method}
\alias{where,CompoundidFilter,missing,missing-method}
\alias{where,CompoundidFilter,SimpleCompoundDb,missing-method}
\alias{where,CompoundidFilter,SimpleCompoundDb,character-method}
\alias{condition<-,CompoundidFilter-method}
\alias{condition,CompoundidFilter-method}
\alias{value<-,CompoundidFilter-method}
\alias{value,CompoundidFilter-method}


\title{Filter objects to retrieve specific data from the database}

\description{
  The filter objects described on this page allow to filter the results
  returned from a \code{\link{SimpleCompoundDb}} database. They inherit
  from the \code{BasicFilter} class defined in the
  \code{ensembldb} database hence providing the same filtering
  framework. All filter objects provide the same methods.

  The individual filters are listed below:
  \describe{
    \item{CompoundidFilter}{
      Enables filtering by the compound id (or compound
      accession). These are provided in the \code{"accession"} column of
      the \code{"compound_basic"} database table.
    }
  }
}

\section{Objects of the class}{
  Filter objects should always be created by the respective constructor
  function, always named as the filter class
  (e.g. \code{CompoundidFilter}).
}
\section{Extends}{
  Class \code{BasicFilter} defined in the
  \code{ensembldb} package, directly.
}

\section{Methods}{

  (in alphabetic order)

  \describe{
    \item{column}{
      \code{signature(object="CompoundidFilter", db="SimpleCompoundDb",
      with.tables="character")}: returns the column (attribute name
      in the \code{\link{SimpleCompoundDb}} database) that is used by the
      filter. The method can be applied to any filter class defined
      in the \code{xcmsExtensions} packages. Argument \code{with.tables}
      allows to restrict the filter to certain database tables (by default
      all database tables are considered).

      The method returns a character string with the name of the column.
    }

    \item{condition, condition<-}{
      Get and set the \code{condition} of the filter (i.e. how values
      should be compared).

      The \code{condition} method returns a character string
      representing the condition.
    }

    \item{value, value<-}{
      Get and set the \code{value} of the filter.

      The \code{value} method returns a character vector of the
      value(s).
    }

    \item{where}{
      \code{signature(object="CompoundidFilter", db="SimpleCompoundDb",
      with.tables="character")}: returns the \emph{where} clause for the
      SQL query to retrieve the data. Argument \code{with.tables}
      allows to restrict the filter to certain database tables (by default
      all database tables are considered).

      The method returns a character string with the name of the column.
    }
  }
}

\value{
  Refer to the method and function description above for detailed
  information on the returned result object.
}
\author{
  Johannes Rainer.
}
\seealso{
  \code{\link{SimpleCompoundDb}}
}
\examples{

## With the package a SimpleCompoundDb object is provided which is bound
## to the variable name scDb.
scDb

## List all tables from the database along with the respective columns
listTables(scDb)

## List all column names
columns(scDb)

## Define a CompoundidFilter specifying some ids (accession ids).
cf <- CompoundidFilter(c("HMDB00010", "HMDB00002", "HMDB00011"))

## The default condition is an equality comparison, i.e. "=" for a single value
## and "in" for multiple values.
condition(cf)

## Use the filter to retrieve just data for these compounds from the database.
res <- compounds(scDb, filter=cf, columns=c("name", "inchikey"))
res

## Change the filter to "!="
condition(cf) <- "!="
condition(cf)

## The query will now return all entries except the data for the 3 compounds.
res <- compounds(scDb, filter=cf, columns=c("name", "inchikey"))
nrow(res)

}
\keyword{classes}

