\name{SimpleCompoundDb-class}
\Rdversion{1.1}
\docType{class}
\alias{scDb}
\alias{SimpleCompoundDb-class}
\alias{SimpleCompoundDb}
\alias{show,SimpleCompoundDb-method}
\alias{columns,SimpleCompoundDb-method}
\alias{as.data.frame,SimpleCompoundDb-method}
\alias{dbconn,SimpleCompoundDb-method}
\alias{listTables}
\alias{listTables,SimpleCompoundDb-method}
%% The /real/ stuff
\alias{compounds}
\alias{compounds,SimpleCompoundDb-method}
\alias{mzmatch,numeric,SimpleCompoundDb-method}
\alias{mzmatch,matrix,SimpleCompoundDb-method}

\title{A simple metabolite compound database}

\description{
  The \code{SimpleCompoundDb} represents simple database to store
  compound information and thus to assist and allow to perform a simple
  identification of compounds based on their mass.

  This package also provides a \code{SimpleCompoundDb object} bound to
  the variable name \code{scDb}.
}

\section{Objects of the class}{
  \code{SimpleCompoundDb} objects should only be created using the
  \code{SimpleCompoundDb} function, that requires the file name of the
  SQLite database providing the data.
}
\usage{

\S4method{as.data.frame}{SimpleCompoundDb}(x, row.names=NULL,
                                           optional=FALSE, ...)

\S4method{columns}{SimpleCompoundDb}(x)

\S4method{compounds}{SimpleCompoundDb}(x, columns, filter=list(),
                                       order.by="accession", ...)

\S4method{dbconn}{SimpleCompoundDb}(x)

\S4method{listTables}{SimpleCompoundDb}(x, ...)

\S4method{mzmatch}{numeric,SimpleCompoundDb}(x, mz, mzdev=0, ppm=10,
                                             column="monoisotopic_molecular_weight",
                                             ionAdduct=NULL)

\S4method{mzmatch}{matrix,SimpleCompoundDb}(x, mz, mzdev=0, ppm=10,
                                             column="monoisotopic_molecular_weight",
                                             ionAdduct=NULL)

SimpleCompoundDb(x)

}
\arguments{

  (in alphabetic order)

  \item{column}{
    For \code{mzmatch}: the name of the column in the database against
    which the provided M/Z values should be matched. Allowed are either
    \code{"avg_molecular_weight"} or
    \code{"monoisotopic_molecular_weight"}.
  }

  \item{columns}{
    For \code{compounds}: a character vector with names of columns that
    should be returned from the database. By default all columns from
    the \code{compound_basic} table are returned.
  }

  \item{filter}{
    For \code{compounds}: a single filter object of a list of filter
    objects. The filter objects should extend the
    \code{BasicFilter} class defined in the
    \code{ensembldb} package. See help for
    \code{\link{CompoundidFilter}} for supported filters.
  }

  \item{ionAdduct}{
    The name(s) of the suspected ion adducts of the measured M/Z
    values. By default (\code{ionAdduct=NULL}) it is assumed that the
    searched compount is already an ion. See
    \code{\link{supportedIonAdducts}} for a complete list of available
    ion adducts.
  }

  \item{mz}{
    For \code{mzmatch}: the \code{SimpleCompoundDb} instance.
  }

  \item{mzdev}{
    Numeric (length 1) specifying an allowed maximal difference of the
    M/Z values. See \code{\link{mzmatch}} for more details.
  }

  \item{optional}{
    For \code{as.data.frame}: not used.
  }

  \item{order.by}{
    For \code{compounds}: the column by which the result should be
    ordered.
  }

  \item{ppm}{
    Numeric (length 1) specifying the parts per million
    deviation/difference of the M/Z value. See \code{\link{mzmatch}} for
    more details.
  }

  \item{row.names}{
    For \code{as.data.frame}: not used.
  }

  \item{x}{
    For \code{SimpleCompoundDb}: a character string specifying the
    SQLite database file name providing the annotations.

    For \code{mzmatch}: a numeric vector with the M/Z values that should
    be matched against the database or a matrix with two columns (or
    columns named \code{"mzmin"} and \code{"mzmax"}) specifying the
    minimum and maximum M/Z of the peak.

    For all other methods: a \code{SimpleCompoundDb} object.
  }

  \item{...}{
    For \code{compounds}: additional arguments passed to the internal
    \code{getWhat} method.

    Ignored for all other functions.
  }

}
\section{Slots}{
  \describe{
    \item{con}{
      The database connection to the database.
    }

    \item{tables}{
      A \code{list} of the database tables. Names of the \code{list}
      represent the names of the database tables, the values their
      attributes (column names).
    }

    \item{.properties}{
      An internal \code{list} of optional properties.
    }

  }
}

\section{Constructors and alike}{
  \describe{

    \item{SimpleCompoundDb}{
      Constructor function to create a new \code{SimpleCompoundDb}
      instance providing a simplified access to the annotation data
      stored in the SQLite database, which file name has to be provided
      with the single argument \code{x}.
    }

  }
}

\section{Basic data retrieval and usage}{
  \describe{

    \item{as.data.frame}{
      Retrieve the full data stored in the database as a
      \code{data.frame}.
    }

    \item{compounds}{
      Retrieve data from the database starting at the
      \code{compound_basic} database table. By default this method
      returns data from that table, but, depending on the \code{columns}
      attribute, might also return data from other tables, that are
      joined to the \code{compound_basic} table. Tables are joined using
      a \emph{left join} starting from the \code{compound_basic} table,
      thus containing all data from that table.
      The results can be filtered by specifying one or more filters
      (defined with objects extending the
      \code{BasicFilter} from the \code{ensembldb}
      package, e.g. the \code{\link{CompoundidFilter}}).

      The method returns a \code{data.frame} with the results of the
      query.
    }

    \item{dbconn}{
      Returns the \code{SQLiteConnection} object
      providing the connection to the actual database.
    }

    \item{columns}{
      Get the column names from the database tables. Returns a character
      vector with the names of the columns.
    }

    \item{listTables}{
      Get the tables and their columns from the database. The method
      returns a \code{list}, the list names being the database table
      names and the values their attribute (database table column)
      names.
    }
  }
}

\section{Compound identification}{
  \describe{
    \item{mzmatch}{
      Matches the given M/Z values with the masses of
      compounds in the database and returns the matches if the
      \code{"deltaMz"} was smaller than the specified threshold. By
      default the method assumes that the M/Z corresponds to the mass of
      the compound, i.e. that the measured feature is already an
      ion. With the argument \code{ionAdduct} it is possible to specify
      the assumed ion adduct of the real compound that is measured in
      the MS. It is possible to specify single adduct names, or all
      possible (i.e. the most commonly found) ion adducts according to
      the Fiehn lab's ESI MS adduct calculator [Huang 1999]. In that case, the mass
      of the possible adduct is calculated for each input M/Z and this
      mass is looked up in the database. The \code{ppm} and \code{mzdev}
      arguments are used for the mass search.

      If the input argument \code{x} is a matrix specifying the minimal
      and maximal M/Z value for a peak matching is performed against
      this range, i.e. the minimum and maximum range is first converted
      to the corresponding mass and compounds in the database are
      searched with a mass within the such defined mass range which is
      extended by the specified ppm. In that case, the reported
      \code{"deltaMz"} in the result table is also the smallest distance
      of the compound mass to the mass of the peak (i.e. to the minimal
      and maximal mass for the \code{minmz} and \code{maxmz}).

      Like the
      generic \code{\link{mzmatch}} method, this method returns a
      \code{list}, each of the elements representing the result for one of
      the specified masses provided as a two-column matrix \code{"idx"}
      and \code{deltaMz}, \code{"idx"} containing the compound IDs,
      \code{"deltaMz"} the difference of the specified mass and the
      compound's M/Z. The ion adduct name is provided as third column
      (\code{"adduct"}). For more details see the \code{\link{mzmatch}}
      method help.
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
  \code{\link[xcms]{xcmsSet}}, \code{\link[xcms]{xcmsRaw}},
  \code{\link{MSdata}}, \code{\link{MSsliceList}},
  \code{\link{CompoundidFilter}}, \code{\link{mzmatch}},
  \code{\link{supportedIonAdducts}}
}
\references{
  Huang N, Siegel MM, Kruppa GH & Laukien FH (1999) Automation of a
  Fourier transform ion cyclotron resonance mass spectrometer for
  acquisition, analysis, and e-mailing of high-resolution exact-mass
  electrospray ionization mass spectral data. \emph{Journal of the
    American Society for Mass Spectrometry}. pp1166-1173.

}
\examples{

## With the package a SimpleCompoundDb object is provided which is bound
## to the variable name scDb.
scDb

## List all tables from the database along with the respective columns
listTables(scDb)

## List all column names
columns(scDb)

## Retrieve all data as a data.frame
tmp <- as.data.frame(scDb)
nrow(tmp)
head(tmp)

## Retrieve compound data for specific compounds only.
## Define the CompoundidFilter
cf <- CompoundidFilter(c("HMDB00010", "HMDB00002", "HMDB00011"))
res <- compounds(scDb, filter=cf, columns=c("name", "inchikey"))
res


## Perform a compound identification.
## Define masses which we would like to identify, i.e. match against the
## compounds in the database. In the default case it is expected that the
## measured features are ions and that the M/Z corresponds already to the
## mass.
comps <- c(300.1898, 298.1508, 491.2000, 169.13481, 169.1348)
Res <- mzmatch(comps, scDb)

## We get a list of results, each list element representing the result of one
## of the specified masses. If the provided mass does not match any compound
## in the database NA is returned.
Res

## Getting the compound name for the 2nd mass.
cf <- CompoundidFilter(Res[[2]][, 1])
res <- compounds(scDb, filter=cf, columns=c("name", "inchikey", "monoisotopic_molecular_weight"))
res

## Next we match the masses of all possible positive ion adducts that would result in the input M/Z
Res <- mzmatch(comps, scDb, ionAdduct=supportedIonAdducts(charge="pos"))
## Removing adducts that would not match.
lapply(Res, function(x){
    return(x[!is.na(x[, 1]), ])
})

## Perform some basic SQL queries.
require(RSQLite)
tmp <- dbGetQuery(dbconn(scDb), "select * from compound_basic where avg_molecular_weight < 100;")
nrow(tmp)
head(tmp)

}
\keyword{classes}

