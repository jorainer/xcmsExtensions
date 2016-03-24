\name{MSdata-class}
\Rdversion{1.1}
\docType{class}
\alias{MSdata-class}
\alias{MSdata}
\alias{show,MSdata-method}
\alias{rtime<-}
\alias{rtime<-,MSdata-method}
\alias{rtime,MSdata-method}
\alias{intensity<-}
\alias{intensity<-,MSdata-method}
\alias{intensity,MSdata-method}
\alias{mz<-,MSdata-method}
\alias{mz,MSdata-method}


\title{Simple MS data object}

\description{
  The \code{MSdata} is a simple object to store intensity values for m/z
  - retention time tuples from a complete or only a subset of an MS run.
}

\section{Objects of the class}{
  \code{MSdata} objects are created by the constructor function
  \code{MSdata}.
}
\usage{

\S4method{intensity}{MSdata}(object)

\S4method{intensity}{MSdata}(object) <- value

MSdata(mz, rtime, intensity, mslevel=1)

\S4method{mz}{MSdata}(object)

\S4method{mz}{MSdata}(object) <- value

\S4method{rtime}{MSdata}(object, return.type="numeric")

\S4method{rtime}{MSdata}(object) <- value

}
\arguments{

  (in alphabetic order)

  \item{intensity}{
    A numeric vector specifying the intensity values.
  }

  \item{mz}{
    A numeric vector specifying the m/z values.
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

  \item{value}{
    For \code{intensity<-}, \code{mz<-} and \code{rtime<-}, a numeric
    vector with the intensity, mz or retention time values,
    respectively.
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
  }
}

\section{Accessors}{
  \describe{

    \item{intensity, intensity<-}{
      Get or set the intensity values of the \code{MSdata}
      object.
    }

    \item{mz, mz<-}{
      Get or set the mz values of the \code{MSdata} object.
    }

    \item{rtime, rtime<-}{
      Get or set the retention times of the \code{MSdata}
      object. Depending on the \code{return.type} argument, the
      retention time is returned as a numeric vector or as an
      \code{\link[S4Vectors]{Rle}} object.
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
  \code{\link[xcms]{xcmsRaw}}
}
\examples{


}
\keyword{classes}

