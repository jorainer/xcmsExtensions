\name{MSdata-class}
\Rdversion{1.1}
\docType{class}
\alias{MSdata-class}
\alias{MSdata}
\alias{msData}
\alias{as.matrix,MSdata-method}
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
\alias{getChromatogram}
\alias{getChromatogram,MSdata-method}
\alias{plotChromatogram}
\alias{plotChromatogram,MSdata-method}
\alias{getSpectrum}
\alias{getSpectrum,MSdata-method}
\alias{plotSpectrum}
\alias{plotSpectrum,MSdata-method}


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

\S4method{getChromatogram}{MSdata}(object, FUN=max, bins=NULL,
                                   nbin=NULL, binSize=NULL, ...)

\S4method{getSpectrum}{MSdata}(object, FUN=max, bins=NULL,
                               nbin=NULL, binSize=NULL, ...)

\S4method{intensity}{MSdata}(object)

\S4method{intensity}{MSdata}(object) <- value

\S4method{intrange}{MSdata}(object)

MSdata(mz, rtime, intensity, mslevel=1)

\S4method{msData}{xcmsRaw}(object, mzrange=NULL, rtrange=NULL, ...)

\S4method{mz}{MSdata}(object)

\S4method{mz}{MSdata}(object) <- value

\S4method{mzrange}{MSdata}(object)

\S4method{plotChromatogram}{MSdata}(object, FUN=max, bins=NULL,
                                    nbin=NULL, binSize=NULL, add=FALSE,
                                    main=paste(format(mzrange(object),
                                    2), collapse="-"),
              xlab="Retention time", ylab="Intensity", ...)

\S4method{plotSpectrum}{MSdata}(object, FUN=max, bins=NULL,
                                nbin=NULL, binSize=NULL, add=FALSE,
                                main=paste(format(rtrange(object),
                                    2), collapse="-"),
              xlab="M/Z", ylab="Intensity", ...)

\S4method{rtime}{MSdata}(object, return.type="numeric")

\S4method{rtime}{MSdata}(object) <- value

\S4method{rtrange}{MSdata}(object)

}
\arguments{

  (in alphabetic order)

  \item{add}{
    For \code{plotChromatogram} and \code{plotSpectrum}: if the plot should be added to an
    existing plot.
  }

  \item{bins}{
    For \code{getChromatogram}, \code{getSpectrum},
    \code{plotChromatogram} and \code{plotSpectrum}:
    Numeric vector specifying the bins to bin the retention time. See
    \code{getChromatogram} method description for details on the
    binning.
  }

  \item{binSize}{
     For \code{getChromatogram}, \code{getSpectrum},
    \code{plotChromatogram} and \code{plotSpectrum}:
    Numeric (length 1) specifying the bin width/size. See
    \code{getChromatogram} method description for details on the
    binning.
  }

  \item{FUN}{
     For \code{getChromatogram}, \code{getSpectrum},
    \code{plotChromatogram} and \code{plotSpectrum}: function
    that should be used to aggregate intensities within a bin. By
    default the largest intensity within a bin is selected. See
    \code{getChromatogram} method description for details on the
    binning.
  }

  \item{intensity}{
    A numeric vector specifying the intensity values.
  }

  \item{main}{
    For \code{plotChromatogram} and \code{plotSpectrum}: the title for the plot.
  }

  \item{mslevel}{
    For \code{MSdata}: the MS level.
  }

  \item{mz}{
    A numeric vector specifying the m/z values.
  }

  \item{mzrange}{
    For \code{msData}: a numeric of length 2 specifying the range of m/z
    values to be extracted from the object.
  }

  \item{nbin}{
    For \code{getChromatogram}, \code{getSpectrum},
    \code{plotChromatogram} and \code{plotSpectrum}: numeric
    (length 1) specifying the number of bins in which the retention time
    range should be binned. See \code{getChromatogram} method
    description for details on the binning.
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

  \item{xlab}{
    For \code{plotChromatogram} and \code{plotSpectrum}: the label for the x-axis.
  }

  \item{ylab}{
    For \code{plotChromatogram} and \code{plotSpectrum}: the label for the y-axis.
  }

  \item{...}{
    For \code{plotChromatogram} and \code{plotSpectrum}: additional graphical properties to be
    passed to the \code{plot} function.

    For all other methods: not used.
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

    \item{getChromatogram}{
      Extracts the chromatogram from the \code{MSdata} (i.e. intensities
      along the retention time dimension). The method returns a numeric
      matrix with two columns, \code{"rtime"} and
      \code{"intensity"}, the rows being ordered increasingly by
      retention time.

      Since the \code{MSdata} is in principle
      3-dimensional (intensities in the m/z - rt space), values have to
      be aggregated in the m/z dimension. By default the method only
      aggregates intensities along the m/z dimension with the same
      retention time value. Specifying either one of \code{bins},
      \code{nbin} and \code{binSize} it is in addition possible to
      bin the retention times into discrete, equal spaced, bins. Thus,
      all intensities falling into one bin (i.e. retention time within
      the bin and all intensities along the m/z dimension present in the
      \code{MSdata}) are aggregated using the function \code{FUN}, which
      by default just selects the maximal signal in each bin.

      The binning in rt-dimension is performed as follows: if
      \code{nbins} is specified, the range of retention times
      (\code{rtrange}) is split into \code{nbins} equal sized bins and
      all intensities with an retention time within a certain bin are
      aggregated (including all in m/z dimension; within
      \code{mzrange}). For \code{binSize}, the range is split into bins
      with the width/size \code{binSize}. Argument \code{bins} allows to
      define the bins externally and thus ensure an identical binning
      e.g. across \code{MSdata} objects from different samples.

      How the aggregation of intensity values along the m/z range (and
      eventually also for retention times within the same retention time
      bin) is performed can be specified with the \code{FUN}
      argument. By default, the maximal signal is selected
      (i.e. \code{FUN=max}). For a \code{MSdata} object representing a full MS
      run, this would correspond to a \emph{base peak chromatogram}
      (BPC), while using \code{FUN=sum} would result in an \emph{total
      ion chromatogram} (TIC) [Smith 2014].
    }

    \code{getSpectrum}{
      Extracts the spectrum from the \code{MSdata} (i.e. intensities
      along the M/Z dimension, usually for a single retention time or a
      small retention time range). The method returns a numeric matrix
      with two columns, \code{"mz"} and \code{"intensity"} with the rows
      being increasingly ordered by the mz value.

      By default, the method aggregates all intensity values along the
      retention time dimension (available in the \code{MSdata} object)
      for a single M/Zvalue. The argument \code{FUN} allows to specify
      the function that should be used to aggregate intensities.
      In addition, the arguments \code{bins}, \code{nbin} and
      \code{binSize} allow to bin M/Z values and thus to aggregate also
      intensity values of M/Z values falling within a certain bin. See
      the documentation of the \code{getChromatogram} method above for
      more details.
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

\section{Data visualizations}{
  \describe{
    \item{plotChromatogram}{
      The method plots the chromatogram of an \code{MSdata}
      object which is internally extracted by the \code{getChromatogram}
      method (refer to the description of that function for the data
      binning etc.).
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
\references{
  Smith R, Mathis AD, Ventura D and Prince JT. (2014) Proteomics,
  lipidomics, metabolomics: a mess spectrometry tutorial from a computer
  scientists's point of view. \emph{BMC Bioinformatics}, 15(Suppl 7):59.
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


## Extracting the chromatogram from an MSdata.
## Defining the mz and rt ranges. These correspond approx. to the
## peak in the xcms vignette (Figure 4).
mzr <- c(302, 302.5)
rtr <- c(2500, 2700)
msd <- msData(xraw, rtrange=rtr, mzrange=mzr)

## How many data points do we have?
length(rtime(msd))
## Hom many unique retention times?
length(unique(rtime(msd)))
## Thus, we have some intensities with the same retention time (but different
## m/z).

####-------------
## Extracting the chromatogram.
chr <- getChromatogram(msd)
nrow(chr)
## By default, the method aggregates intensities with identical retention time.

## Plotting the chromatogram
plotChromatogram(msd, type="l")

## Next we extract a chromatogram with binning on the retention time. Thus we
## set to split the retention time range into 50 equal spaced bins and to aggregate
## the intensities within each of the bins. As aggregation function we use "max", i.e.
## we select the maximal intensity for each bin.
chr2 <- getChromatogram(msd, FUN=max, nbin=50)
nrow(chr2)

## Adding that chromatogram to the existing plot.
plotChromatogram(msd, FUN=max, nbin=50, add=TRUE, col="blue", type="b")

## Generating the base peak chromatogram (BPC) and the total ion chromatogram (TIC)
fulld <- msData(xraw)
par(mfrow=c(1, 2))
plotChromatogram(fulld, FUN=max, type="l", main="BPC")
plotChromatogram(fulld, FUN=sum, type="l", main="TIC")


####-------------
## Extract a spectrum.
## Similar to the getChromatogram example above we can simply extract a
## spectrum using the getSpectrum method.
spc <- getSpectrum(msd)
nrow(spc)
## The spectrum for this peak is rather small, which is not surprising, since
## we restricted the m/z range. To get the spectrum for the compound corresponding
## to the peak, we restrict a narrow retention time range around the peak, and a
## larget mz range.
rtr <- c(2560, 2640)
mzr <- c(300, 330)
## Get the "slice" of MS data defined by this ranges
msd <- msData(xraw, mzrange=mzr, rtrange=rtr)

## Get the (unbinned) spectrum; intensities with the same m/z but different
## retention times are aggregated (i.e. the maximal signal is taken).
spc <- getSpectrum(msd)
nrow(spc)

## Next we plot the chromatogram and the spectrum for the MS data slice.
par(mfrow=c(1, 2))
plotChromatogram(msd, type="l")
plotSpectrum(msd, type="l")

}
\keyword{classes}

