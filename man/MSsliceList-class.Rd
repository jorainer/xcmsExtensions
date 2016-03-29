\name{MSsliceList-class}
\Rdversion{1.1}
\docType{class}
\alias{MSsliceList-class}
\alias{MSsliceList}
\alias{[,MSsliceList,ANY,ANY,ANY-method}
\alias{[[,MSsliceList,ANY,ANY-method}
\alias{show,MSsliceList-method}
\alias{intrange,MSsliceList-method}
\alias{intranges}
\alias{intranges,MSsliceList-method}
\alias{mzrange,MSsliceList-method}
\alias{mzranges}
\alias{mzranges,MSsliceList-method}
\alias{rtrange,MSsliceList-method}
\alias{rtranges}
\alias{rtranges,MSsliceList-method}
\alias{length,MSsliceList-method}
\alias{slices,MSsliceList-method}
\alias{slices<-,MSsliceList-method}
\alias{slices}
\alias{slices<-}

\title{List of MSslice objects}

\description{
  The \code{MSsliceList} allows to group several \code{\link{MSslice}}
  objects. Each of the \code{MSslice} objects is supposed to define a
  different slice of the MS data (defined by a m/z values and
  retention time range) and to contain the intensity values for the
  slice across several samples. The \code{MSsliceList} object can then
  contain the data of several different such slices. For example, an
  \code{\link{MSslice}} object could contain the intensity values for a
  peak (defined by its mz and rt range) across all samples of an
  experiment (i.e. an \code{\link[xcms]{xcmsSet}} object). The
  \code{MSsliceList} could then contain the values for several or all
  peaks identified in an experiment.
}

\section{Objects of the class}{
  \code{MSsliceList} objects can be created by the constructor
  function \code{MSsliceList} but by default they are supposed to be
  returned by the \code{\link{msSlice}} method for \code{\link[xcms]{xcmsSet}}
  objects with \code{mzrange} and \code{rtrange} parameters being
  matrices specifying different slices of data.

}
\usage{

\S4method{intrange}{MSsliceList}(object)

\S4method{intranges}{MSsliceList}(object)

\S4method{length}{MSsliceList}(x)

MSsliceList(...)

\S4method{mzrange}{MSsliceList}(object)

\S4method{mzranges}{MSsliceList}(object)

\S4method{rtrange}{MSsliceList}(object)

\S4method{rtranges}{MSsliceList}(object)

\S4method{slices}{MSsliceList}(object)

\S4method{slices}{MSsliceList}(object) <- value

}
\arguments{

  (in alphabetic order)

  \item{object}{
    An \code{MSsliceList} object.
  }

  \item{x}{
    For \code{length}: an \code{MSsliceList} object.
  }

  \item{value}{
    For \code{slices}: a list of \code{\link{MSslice}} objects.
  }

  \item{...}{
    For \code{MSsliceList}: \code{\link{MSslice}} objects for the
    \code{MSsliceList} object.
    Ignored for all other methods.
  }

}
\section{Slots}{
  \describe{
    \item{slices}{
      A list of \code{\link{MSslice}} objects.
    }
  }
}

\section{Constructors and alike}{
  \describe{

    \item{MSsliceList}{
      Constructor method to create a new \code{MSsliceList} instance
      accepting either a list of \code{\link{MSslice}} objects or
      individual \code{\link{MSslice}} objects.
    }

    \item{msSlice}{
      If the provided parameters \code{mzrange} and \code{rtrange}
      define multiple slices, the method returns a \code{MSsliceList}
      object from an \code{\link[xcms]{xcmsRaw}} or
      \code{\link[xcms]{xcmsSet}} object. See the help of
      \code{\link{msSlice}} for more information.
    }

  }
}

\section{Accessing and sub-setting}{
  \describe{

    \item{intrange}{
      Get the range of the intensity value ranges of the individual
      slices within the object.
    }

    \item{intranges}{
      Returns a list of the intensity range for each individual slice
      within the object.
    }

    \item{length}{
      Get the length of the object, i.e. the number of \code{MSslice}
      objects stored within the object.
    }

    \item{mzrange}{
      Get the range of m/z value ranges of the individual slices within
      the object.
    }

    \item{mzrange}{
      Returns a list of the m/z value ranges of the individual slices within
      the object.
    }

    \item{object[i]}{
      Subsets the \code{MSsliceList} object \emph{object} and returns a
      \code{MSsliceList} with only the elements selected by \code{i}.
    }

    \item{object[[i]]}{
      Extract a single \code{\link{MSslice}} object from the
      \code{MSsliceList}.
    }

    \item{rtrange}{
      Get the range of retention time ranges of the slices within the
      object.
    }

    \item{rtranges}{
      Returns a list of the retention time ranges of the individual
      slices within the object.
    }

    \item{slices, slices<-}{
      Getter and setter method for the \code{slices} slot.
      Returns the list of \code{\link{MSslice}} objects.
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
  \code{\link{MSdata}}, \code{\link{MSslice}}
}
\examples{

## Load the faahKO package with the test data for the xcms package.
library(faahKO)
xset <- faahko

## Define several data sub-sets and extract the corresponding values from
## the xcmsSet object. This will extract the intensity values of the
## defined slices from each individual sample.
rtr <- rbind(c(2600, 2670),
             c(2900, 2920)
             )
mzr <- rbind(c(300, 400),
             c(401, 402)
             )
msl <- msSlice(xset, rtrange=rtr, mzrange=mzr)
msl

## Get the number of slices
length(msl)

## Get the mz, rt and intensity range. This corresponds to the range of
## ranges of the individual slices across all samples.
mzrange(msl)
rtrange(msl)
intrange(msl)

## Get the ranges for each individual slice
mzranges(msl)
rtranges(msl)
intranges(msl)

## Sub-set the MSsliceList object.
msl[2]

## Extract an individual MSslice object
msl[[2]]

}
\keyword{classes}

