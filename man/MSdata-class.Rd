\name{MSdata-class}
\Rdversion{1.1}
\docType{class}
\alias{MSdata-class}
\alias{MSdata}
\alias{msData}
\alias{as.matrix,}
\alias{msData,xcmsRaw-method}
\alias{show,MSdata-method}
\alias{rtime<-}
\alias{rtime<-,MSdata-method}
\alias{rtime,MSdata-method}
\alias{rtrange,MSdata-method}
\alias{intensity<-}
\alias{intensity<-,MSdata-method}
\alias{intensity,MSdata-method}
\alias{intrange}
\alias{intrange,MSdata-method}
\alias{mz<-,MSdata-method}
\alias{mz,MSdata-method}
\alias{mzrange,MSdata-method}


\title{Simple MS data object}

\description{
  The \code{MSdata} is a simple object to store intensity values for m/z
  - retention time tuples from a complete or only a subset of an MS run.
}

\section{Objects of the class}{
  \code{MSdata} objects are created by the constructor function
  \code{MSdata}, or the \code{msData} method e.g. to extract data from
  an \code{\link[xcms]{xcmsRaw}} object.
}
\usage{

\S4method{as.matrix}{MSdata}(x, ...)

\S4method{intensity}{MSdata}(object)

\S4method{intensity}{MSdata}(object) <- value

\S4method{intrange}{MSdata}(object)

MSdata(mz, rtime, intensity, mslevel=1)

\S4method{msData}{xcmsRaw}(object, mzrange=NULL, rtrange=NULL, ...)

\S4method{mz}{MSdata}(object)

\S4method{mz}{MSdata}(object) <- value

\S4method{mzrange}{MSdata}(object)

\S4method{rtime}{MSdata}(object, return.type="numeric")

\S4method{rtime}{MSdata}(object) <- value

\S4method{rtrange}{MSdata}(object)

}
\arguments{

  (in alphabetic order)

  \item{intensity}{
    A numeric vector specifying the intensity values.
  }

  \item{mz}{
    A numeric vector specifying the m/z values.
  }

  \item{mzrange}{
    For \code{msData}: a numeric of length 2 specifying the range of m/z
    values to be extracted from the object.
  }

  \item{object}{
    For all methods a \code{MSdata} object.
  }

  \item{return.type}{
    For \code{rtime}: either \code{"numeric"} or \code{"Rle"} to return
    the retention time as a numeric vector, or as a
    \code{\link[S4Vectors]{Rle}} object.
  }

  \item{rtime}{
    A numeric vector specifying the retention time values.
  }

  \item{rtrange}{
    For \code{msData}: a numeric of length 2 specifying the retention
    time range to be extracted from the object.
  }

  \item{value}{
    For \code{intensity<-}, \code{mz<-} and \code{rtime<-}, a numeric
    vector with the intensity, mz or retention time values,
    respectively.
  }

  \item{x}{
    For \code{as.matrix}: a \code{MSdata} object.
  }

  \item{...}{
    Presently not used.
  }

}
\section{Slots}{
  \describe{
    \item{intensity}{
      An integer vector with the intensities for each m/z - retention
      time tuple.
    }

    \item{intrange}{
      A numeric vector of length 2 representing the range of intensity
      values.
    }

    \item{mslevel}{
      The MS level.
    }

    \item{mz}{
      A numeric vector with the m/z values.
    }

    \item{mzrange}{
      A numeric vector of length 2 specifying the range of m/z
      values.
    }

    \item{rtime}{
      An \code{Rle} object with the retention time values.
    }

    \item{rtrange}{
      A numeric vector of length 2 specifying the retention time range.
    }

  }
}

\section{Constructors and alike}{
  \describe{
    \item{MSdata}{
      Constructor method to create a new \code{MSdata} instance.
    }

    \item{msData}{
      Extract MS data from an \code{\link[xcms]{xcmsRaw}} object and
      return the values as a \code{MSdata} object. Arguments
      \code{mzrange} and \code{rtrange} allow to define the data sub-set
      to be extracted. If not specified, the method will extract the
      full data.
    }
  }
}

\section{Accessing data of an \code{MSdata} object}{
  \describe{

    \item{as.matrix}{
      Extracts a matrix with column names \code{"rtime"}, \code{"mz"} and
      \code{"intensity"} with all data values of the \code{MSdata}
      object.
    }

    \item{intensity, intensity<-}{
      Get or set the intensity values of the \code{MSdata}
      object. \code{intensity} returns the intensity values.
    }

    \item{intrange}{
      Get the intensity value range. Returns a numeric vector of length
      2. For an empty \code{MSdata} object without data points it
      returns \code{c(NA, NA)}.
    }

    \item{mz, mz<-}{
      Get or set the mz values of the \code{MSdata} object. \code{mz}
      returns a numeric vector with the mz values.
    }

    \item{mzrange}{
      Get the m/z value range of the \code{MSdata} object. Returns a
      numeric vector of 2 or \code{c(NA, NA)} if the object is empty.
    }

    \item{rtime, rtime<-}{
      Get or set the retention times of the \code{MSdata}
      object. Depending on the \code{return.type} argument, the
      retention time is returned as a numeric vector or as an
      \code{\link[S4Vectors]{Rle}} object.
    }

    \item{rtrange}{
      Get the retention time range of the \code{MSdata} object. Returns a
      numeric vector of 2 or \code{c(NA, NA)} if the object is empty.
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
  \code{\link[xcms]{xcmsRaw}}, \code{\link{MSslice}}
}
\examples{

## Load the faahKO package with the test data for the xcms package.
library(faahKO)
xset <- faahko

## Get the raw data for the first sample.
xraw <- getXcmsRaw(xset, 1)

## Define a data sub-set and extract the corresponding values.
rtr <- c(2600, 2670)
mzr <- c(300, 400)
mssub <- msData(xraw, rtrange=rtr, mzrange=mzr)

## Get the mz value and retention time range
mzrange(mssub)
rtrange(mssub)

## The actual retention time range is different from the specifyed
## range as it depends on the acutally measured scan times.

## Get the intensity range
intrange(mssub)

## Get mz values, retention times and intensities.
head(mz(mssub))
head(rtime(mssub))
head(intensity(mssub))

## Extract a matrix with the full data.
head(as.matrix(mssub))

}
\keyword{classes}

