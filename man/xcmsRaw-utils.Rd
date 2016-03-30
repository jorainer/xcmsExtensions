\name{xcmsRaw-utils}
\Rdversion{1.1}
\alias{scantimes}
\alias{scantimes,xcmsRaw-method}
\alias{getData}
\alias{getData,xcmsRaw,ANY-method}

\title{Various utility functions for xcmsRaw objects}

\description{
  These functions enable user-friendly access of the data stored in an
  \code{xcmsRaw} object.
}
\usage{

\S4method{scantimes}{xcmsRaw}(x)

\S4method{getData}{xcmsRaw}(x, mzrange=NULL, rtrange=NULL, intrange=NULL)

}
\arguments{

  (in alphabetic order)

  \item{mzrange}{
    For \code{getData}: an optional numeric vector (of length 2)
    to subset the data by m/z range. Alternatively a numeric matrix with
    two columns, each row specifying a m/z range to subset the data.
  }

  \item{rtrange}{
    For \code{getData}: an optional numeric vector (of length 2)
    to subset by retention time range. Alternatively a numeric matrix with
    two columns, each row specifying a retention time range to subset
    the data.
  }

  \item{intrange}{
    For \code{getData}: an optional numeric vector (of length 2)
    to subset by intensity range. Alternatively a numeric matrix with
    two columns, each row specifying an intensity range to subset the
    data.
  }

  \item{x}{
    An \code{xcmsRaw} object.
  }

}
\section{Accessing data in \code{xcmsRaw} objects}{

  \describe{

    \item{getData}{
      Extract retention time, mz and intensity values from the
      \code{xcmsRaw} object. If none of the 3 parameters \code{mzrange},
      \code{rtrange} and \code{intrange} are specified the method
      returns the full data. The method supports also definitions of
      multiple m/z, rt or intensity ranges, i.e. if \code{mzrange},
      \code{rtrange} and/or \code{intrange} is a numeric two-column
      matrix, each row specifying a range.

      The method returns a numeric matrix with columns \code{"time"},
      \code{"mz"} and \code{"intensity"} and rows depending on the input
      arguments (or a matrix with \code{nrow} equals 0 if the specified
      ranges do not match any data points). If \code{mzrange},
      \code{rtrange} or \code{intrange} is a matrix, the method returns
      a \code{list} of matrices, each representing the subset defined by
      one of the rows of the input matrix (matrices).
    }

    \item{scantimes}{
      Returns the sample (retention) time for each individual mz value
      stored the \code{xcmsRaw} object. The length of the returned
      numeric vector is thus equal to the number of available data
      points. This is in contrast to the \code{@scantime} slot, which
      stores only the discrete scan times.
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
  \code{\link{MSdata}}, \code{\link{MSslice}}, \code{\link{MSsliceList}}
}
\examples{
## Load the test data.
library(faahKO)
xset <- faahko
## Get the raw data for the first file
xraw <- getXcmsRaw(xset, 1)

## Get the scan time of the individual measured compounts.
head(scantimes(xraw))

## To get the unique scan times
head(xraw@scantime)

## Get the data matrix with retention time, mz and intensity
## for the defined retention time range.
datmat <- getData(xraw, rtrange=c(2509, 2530))
head(datmat)

## Subset also by mzrange.
datmat <- getData(xraw, rtrange=c(2509, 2530), mzrange=c(301, 302.003))
datmat

## Define 3 different ranges and extract the data subsets.
rtm <- rbind(c(2509, 2530),
             c(2534, 2540),
             c(2630, 2645))
mzm <- rbind(c(301, 302.003),
             c(304, 305),
             c(308, 308.3))
## Extracting the data we get a list of length 3, each list element
## representing the data subset defined by the corresponding rows in
## rtm and mzm.
getData(xraw, rtrange=rtm, mzrange=mzm)

}
\keyword{classes}

