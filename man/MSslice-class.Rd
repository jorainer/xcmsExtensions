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
\alias{getChromatogram,MSslice-method}
\alias{plotChromatogram,MSslice-method}
\alias{getSpectrum,MSslice-method}
\alias{plotSpectrum,MSslice-method}


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

\S4method{getChromatogram}{MSslice}(object, FUN=max, bins=NULL,
                           nbin=NULL, binSize=NULL, ...)

\S4method{getSpectrum}{MSslice}(object, FUN=max, bins=NULL,
                                nbin=NULL, binSize=NULL, ...)

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

\S4method{plotChromatogram}{MSslice}(object, FUN=max, bins=NULL,
                            nbin=NULL, binSize=NULL,
                            main=paste(format(mzrange(object), 2),
                            collapse="-"), xlab="Retention time",
                            ylab="Intensity", col=1, lty=1, ...)

\S4method{plotSpectrum}{MSslice}(object, FUN=max, bins=NULL,
                                 nbin=NULL, binSize=NULL,
                                 main=paste(format(rtrange(object), 2),
                                 collapse="-"), xlab="M/Z",
                                 ylab="Intensity", col=1, lty=1, ...)

\S4method{rtrange}{MSslice}(object)

}
\arguments{

  (in alphabetic order)


  \item{bins}{
    For \code{getChromatogram}, \code{getSpectrum},
    \code{plotChromatogram} and \code{plotSpectrum}: numeric
    vector specifying the bins to bin the retention time. See
    \code{getChromatogram} method description for \code{\link{MSdata}}
    objects for more details on the binning.
  }

  \item{binSize}{
    For \code{getChromatogram}, \code{getSpectrum},
    \code{plotChromatogram} and \code{plotSpectrum}: numeric
    (length 1) specifying the bin width/size. See
    \code{getChromatogram} method description for \code{\link{MSdata}}
    objects for more details on the binning.
  }

  \item{FUN}{
    For \code{getChromatogram}, \code{getSpectrum},
    \code{plotChromatogram} and \code{plotSpectrum}: function
    that should be used to aggregate intensities within a bin. By
    default the largest intensity within a bin is selected. See
    \code{getChromatogram} method description for \code{\link{MSdata}}
    objects for more details on the binning.
  }

  \item{BPPARAM}{
    Parallel processing parameters; see
    \code{\link[BiocParallel]{bpparam}} for more information.
  }

  \item{call}{
    Optional argument allowing to specifying the call that used to
    generate the object.
  }

  \item{col}{
    For \code{plotChromatogram} and \code{plotSpectrum}: a single color or vector specifying
    colors for each sample/\code{MSdata} object in the \code{MSslice}.
  }

  \item{lty}{
    For \code{plotChromatogram} and \code{plotSpectrum}: a numeric specifying the line type for
    each sample/\code{MSdata} object in the \code{MSslice}; if
    \code{lty} has length of 1, this line type will be used for all
    lines.
  }

  \item{main}{
    For \code{plotChromatogram} and \code{plotSpectrum}: the title of the plot.
  }

  \item{mzrange}{
    A numeric vector defining the m/z-value range of the slice.
  }

  \item{nbin}{
    For \code{getChromatogram}, \code{getSpectrum},
    \code{plotChromatogram} and \code{plotSpectrum}: numeric
    (length 1) specifying the number of bins in which the retention time
    range should be binned. See \code{getChromatogram} method
    description for \code{\link{MSdata}} objects for more details on the
    binning.
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

  \item{xlab}{
    For \code{plotChromatogram} and \code{plotSpectrum}: the label for the x-axis.
  }

  \item{ylab}{
    For \code{plotChromatogram} and \code{plotSpectrum}: the label for the y-axis.
  }

  \item{...}{
    For \code{MSslice}: \code{\link{MSdata}} objects for the
    \code{MSslice} object.

    For \code{plotChromatogram}: additional graphical parameters to be
    passed down to the base plot function.

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
      or \code{binSize}, to bin values also along the retention time
      dimension. Ideally, the binning should not be too coarse, but also
      not too detailed (in which case the matrix will again be
      relatively sparse); a binning with a \code{binSize=1} (i.e. one
      second), or slightly smaller, might be OK in most instances.

      For more details on the binning see the help for the
      \code{getChromatogram} method of \code{\link{MSdata}} objects.
    }

    \item{getSpectrum}{
      Extracts a numeric matrix with the spectrum (i.e. intensities
      along the M/Z dimension) for each of the \code{MSdata} objects and
      hence samples of the \code{MSslice} object. The rownames of the
      matrix represent the distinct M/Z values and columns the
      individual samples. For more details (e.g. also on binning etc)
      see the help of the \code{getChromatogram} method above.
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
      Plots the chromatogram of the rt-m/z slice across all of the
      \code{MSdata} objects (i.e. samples). The arguments \code{col} and
      \code{lty} can be used to use a different color for each line
      (i.e. sample) or a different line type. See description of the
      \code{getChromatogram} for details on the other parameters.
    }

    \item{plotSpectrum}{
      Plots the spectrum; see help of \code{plotChromatogram} above.
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

## Extract the chromatogram of a m/z-rt slice across all samples.
## Defining the mz and rt ranges for the peak in the xcms
## vignette (Figure 4).
mzr <- c(302, 302.5)
rtr <- c(2500, 2700)
mss <- msSlice(xset, rtrange=rtr, mzrange=mzr)

## We're getting a sparse matrix if we're not binning along the rt-dimension,
## since the retention times of the individual data points are slightly different
## between the samples:
head(getChromatogram(mss))

## We're next binning the values in retention time-bins of size 1 (1 second). By default
## (with FUN=max), the maximal signal is used for each bin.
chrM <- getChromatogram(mss, binSize=1)
dim(chrM)
head(chrM)

## We can also create an image from this matrix. Each column represents the chromatogram
## from one sample, each row one retention time bin with the intensity being color coded.
image(t(chrM))

## Next we plot the chromatogram of the slice across all samples.
## Without any binning:
plotChromatogram(mss, col="#00000080", type="l")

## And with binning.
plotChromatogram(mss, binSize=1, col="#00000080", type="l")

## Extracting the full MS data from the xcmsSet to create base peak chromatogram (BPC)
## and total ion chromatogram (TIC) plot of all samples.
msFull <- msSlice(xset)

## Base peak chromatogram with 1 second binning
plotChromatogram(msFull, binSize=1, col="#00000080", type="l")

## Total ion chromatogram with 1 second binning
plotChromatogram(msFull, FUN=sum, binSize=1, col="#00000080", type="l")

## Extending the M/Z range of the slice and extracting and plotting also the spectrum.
mzr <- c(300, 303)
mss <- msSlice(xset, rtrange=rtr, mzrange=mzr)

## Get the spectrum without binning:
head(getSpectrum(mss))

## We can also extract a slice from an xcmsRaw object, which will then
## contain the data only of that individual sample. Without specifying an
## mzrange and rtrange we extract the full data.
msSlice(getXcmsRaw(xset, 1))

par(mfrow=c(1, 2))
plotChromatogram(mss, binSize=1, col="#00000080", type="l")
plotSpectrum(mss, col="#00000040", type="h")

## Or from an MSdata object; extracting first a MSdata object containing
## the *full* MS data from a single file and creating then an MSslice object
## from that.
msd <- msData(getXcmsRaw(xset, 1))
msSlice(msd)

## Extracting the MSdata from the second sample and passing both to the
## MSslice constructor to generate a two-sample MSslice object.
msd2 <- msData(getXcmsRaw(xset, 2))
MSslice(msd, msd2)

}
\keyword{classes}

