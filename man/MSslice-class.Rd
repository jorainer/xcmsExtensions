\name{MSslice-class}
\Rdversion{1.1}
\docType{class}
\alias{MSslice-class}
\alias{MSslice}
\alias{show,MSslice-method}
\alias{intrange,MSslice-method}
\alias{msSlice}
\alias{msSlice,list-method}
\alias{msSlice,xcmsSet-method}
\alias{msSlice,xcmsRaw-method}
\alias{msSlice,MSdata-method}
\alias{mzrange,MSslice-method}
\alias{rtrange,MSslice-method}
\alias{msData,MSslice-method}
\alias{length,MSslice-method}

\title{MS data subsets of several samples}

\description{
  Sub-sets of an MS run across several samples can be stored in a
  \code{MSslice} object. Such a slice of the MS data is defined by an
  m/z and a retention time range in each sample. The actual data of the
  slice in each sample is stored as \code{\link{MSdata}} object, i.e. as
  intensities for each mz/rt tuple measured in a sample.
  Depending on the measurements, the actual mz/rt ranges in each
  individual \code{\link{MSdata}} can be slightly different, but all are
  within the mz and rt range of the \code{MSslice} object (that can be
  returned by the \code{mzrange} and \code{rtrange} method).
}

\section{Objects of the class}{
  While \code{MSslice} objects can be created by the constructor
  function \code{MSslice} the default way to create such an object is by
  using the \code{msSlice} method on an \code{\link[xcms]{xcmsSet}} (or
  \code{\link[xcms]{xcmsRaw}}) object.
}
\usage{

\S4method{intrange}{MSslice}(object)

\S4method{length}{MSslice}(x)

MSslice(...)

\S4method{msData}{MSslice}(object, ...)

\S4method{msSlice}{list}(object, call=match.call(), ...)

\S4method{msSlice}{MSdata}(object, call=match.call(), ...)

\S4method{msSlice}{xcmsSet}(object, mzrange=NULL, rtrange=NULL,
                            rt="corrected", BPPARAM=bpparam())

\S4method{msSlice}{xcmsRaw}(object, mzrange=NULL, rtrange=NULL, ...)

\S4method{mzrange}{MSslice}(object)

\S4method{rtrange}{MSslice}(object)

}
\arguments{

  (in alphabetic order)

  \item{BPPARAM}{
    Parallel processing parameters; see
    \code{\link[BiocParallel]{bpparam}} for more information.
  }

  \item{call}{
    Optional argument allowing to specifying the call that used to
    generate the object.
  }

  \item{mzrange}{
    A numeric vector defining the m/z-value range of the slice.
  }

  \item{object}{
    For \code{msSlice}: either a \code{list} of \code{\link{MSdata}}
    objects, an \code{\link[xcms]{xcmsSet}} or
    \code{\link[xcms]{xcmsRaw}} object.

    For all other methods a \code{MSslice} object.
  }

  \item{rt}{
    For \code{msSlice}: either \code{"raw"} or \code{"corrected"} to
    specify whether the raw or retention-time corrected scan times
    should be used.
  }

  \item{rtrange}{
    A numeric vector defining the retention time range of the slice.
  }

  \item{x}{
    For \code{length}: an \code{MSslice} object.
  }

  \item{...}{
    For \code{MSslice}: \code{\link{MSdata}} objects for the
    \code{MSslice} object.
    Ignored for all other methods.
  }

}
\section{Slots}{
  \describe{
    \item{data}{
      A list of \code{\link{MSdata}} objects.
    }

    \item{call}{
      The call used to create the object.
    }

    \item{mzrange}{
      A numeric vector of length 2 specifying the range of m/z
      values of the slice.
    }

    \item{rtrange}{
      A numeric vector of length 2 specifying the retention time range
      of the slice.
    }

  }
}

\section{Constructors and alike}{
  \describe{

    \item{MSslice}{
      Constructor method to create a new \code{MSslice} instance
      accepting either a list of \code{\link{MSdata}} objects or
      individual \code{\link{MSdata}} objects.
    }

    \item{msSlice}{
      Extract an \code{msSlice} object from an
      \code{\link[xcms]{xcmsRaw}} or \code{\link[xcms]{xcmsSet}}
      object. The latter (i.e. extraction from an \code{xcmsSet} object)
      is the \emph{default} constructor that will be mostly used. It
      allows to extract intensity values for the m/z and rt range that
      defines the slice of MS data from each individual sample. The
      result is returned as an \code{MSslice} object. If multiple mz and
      rt ranges are specified (i.e. arguments \code{mzrange} and
      \code{rtrange} are 2-column matrices, each row specifying a
      different range) the result will be returned as an
      \code{\link{MSsliceList}} object.

      Note that the \code{mzrange} and \code{rtrange} of each of the
      extracted \code{\link{MSdata}} for each sample will be slightly
      different, but all are \emph{within} the \code{mzrange} and
      \code{rtrange} of the \code{MSslice} object.

      If \code{object} is an \code{MSdata} object it wraps that into an
      \code{MSslice} object.
    }

  }
}

\section{Accessing data of an \code{MSslice} object}{
  \describe{

    \item{getChromatogram}{
      Extracts a numeric matrix with the chromatogram (i.e. intensities
      along the retention time dimension) for each of the \code{MSdata}
      object and hence samples of the \code{MSslice} object. The
      rownames of the matrix represent the (unique) retention times and
      each column contains the spectrum of one sample/\code{MSdata}
      object (i.e. the intensities aggregated in m/z direction for the
      respective retention time using the aggregation function specified
      with argument \code{FUN}).

      In the default settings, without specifying any of the arguments
      \code{bins}, \code{nbin} or \code{binSize}, the method returns a
      matrix with one row for each unique retention time across all of
      the samples. Each row will only contain a numeric value for a
      sample in which a signal was measured for that specific retention
      time. In most instances it will be better to specify \code{nbin}
      or \code{binSize}, to bin
    }

    \item{intrange}{
      Get the range of intensity values of the slice. This represents the range of
      values across all samples (i.e. all \code{MSdata} objects within
      the \code{MSslice}). Returns a numeric of length
      2, or \code{c(NA, NA)} if the object is empty (or the internal
      \code{MSdata} object(s) don't contain any data).
    }

    \item{length}{
      Get the length of the object, i.e. the number of \code{MSdata}
      object stored within the object.
    }

    \item{msData}{
      Get the \code{list} of \code{\link{MSdata}} objects.
    }

    \item{mzrange}{
      Get the m/z value range of the slice. This represents the range of
      values across all samples (i.e. all \code{MSdata} objects within
      the \code{MSslice}). Returns a numeric of length
      2, or \code{c(NA, NA)} if the object is empty (or the internal
      \code{MSdata} object(s) don't contain any data).
    }

    \item{rtrange}{
      Get the retention time range of the slice. This represents the range of
      values across all samples (i.e. all \code{MSdata} objects within
      the \code{MSslice}). Returns a numeric of length
      2, or \code{c(NA, NA)} if the object is empty (or the internal
      \code{MSdata} object(s) don't contain any data).
    }
  }
}

\section{Data visualizations}{
  \describe{
    \item{plotChromatogram}{
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
  \code{\link{MSdata}}, \code{\link{MSsliceList}}
}
\examples{

## Load the faahKO package with the test data for the xcms package.
library(faahKO)
xset <- faahko

## Define a data sub-set and extract the corresponding values from
## the xcmsSet object. This will extract the intensity values of the
## defined slice from each individual sample.
rtr <- c(2600, 2670)
mzr <- c(300, 400)
mss <- msSlice(xset, rtrange=rtr, mzrange=mzr)

## Get the number of samples
length(mss)

## Get the mz, rt and intensity range. This corresponds to the range of
## values of the slice across all samples.
mzrange(mss)
rtrange(mss)
intrange(mss)

## Extract the list of MSdata objects
msData(mss)

## We can also extract a slice from an xcmsRaw object, which will then
## contain the data only of that individual sample. Without specifying an
## mzrange and rtrange we extract the full data.
msSlice(getXcmsRaw(xset, 1))

## Or from an MSdata object; extracting first a MSdata object containing
## the *full* MS data from a single file and creating then an MSslice object
## from that.
msd <- msData(getXcmsRaw(xset, 1))
msSlice(msd)

## Extracting the MSdata from the second sample and passing both to the
## MSslice constructor to generate a two-sample MSslice object.
msd2 <- msData(getXcmsRaw( xset, 2))
MSslice(msd, msd2)

}
\keyword{classes}

