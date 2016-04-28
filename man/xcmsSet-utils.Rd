\name{xcmsSet-utils}
\Rdversion{1.1}
\alias{visualizeXcmsSet}
\alias{peakGroupSummary}
\alias{peakGroupSummary,xcmsSet-method}

\title{Various utility functions for xcmsSet objects}

\description{
  These functions provide additional functionality for \code{xcmsSet}
  objects from the \code{xcms} package.
}
\usage{

\S4method{peakGroupSummary}{xcmsSet}(object, class, value="into")

visualizeXcmsSet(...)

}
\arguments{

  (in alphabetic order)

  \item{class}{
    For \code{peakGroupSummary}: optional factor specifying groups of
    samples for which the summary should be performed. If not specified
    (the default) the factor returned by the \code{sampclass} method is
    used, i.e. peak group summaries are calculated for each sample
    class.
  }

  \item{object}{
    An \code{xcmsSet} object.
  }

  \item{value}{
    For \code{peakGroupSummary}: character string specifying the value
    on which the summary should be calculated. Defaults to
    \code{"into"}, i.e. summary is calculated on the original (raw) peak
    area (intensity). See \code{groupval} method for \code{xcmsSet}
    objects for more information.
  }

  \item{...}{
    Presently not used.
  }
}
\section{Data visualization}{

  \describe{

    \item{visualizeXcmsSet}{
      This function starts the \code{shiny} app to visualize
      \code{xcmsSet} objects present in the current workspace.
    }

  }
}

\section{Data summaries and quality controls}{
  \describe{

    \item{peakGroupSummary}{
      Calculates and returns peak group summaries, i.e. the minimal,
      maximal and mean intensity, the RSD (relative standard deviation
      of coefficient of variationy) and the proportion of samples in which the
      peak was identified. These summaries are calculated for groups of
      samples specified with the \code{class} argument, which defaults
      to the classes defined in the \code{xcmsSet} (returned by the
      \code{sampclass} method).

      The method returns a matrix, each row representing a peak group,
      with columns \code{"minInt"}, \code{"maxInt"}, \code{"meanInt"},
      \code{"RSD"} and \code{"propPresent"} for each of the sample
      classes. The ordering of the rows corresponds to the ordering of
      peak groups returned by the \code{groupval} or \code{groups}
      methods.
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
\keyword{classes}

