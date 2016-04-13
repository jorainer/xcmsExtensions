\name{mzmatch}
\alias{mzmatch}
\alias{mzmatch,numeric,numeric-method}
\alias{mzmatch,matrix,numeric-method}
\title{
  Match M/Z values given a certain resolution
}
\description{
  The \code{mzmatch} function enables matching of M/Z values allowing
  also small differences between the values, which can be set with the
  \code{mzdev} and \code{ppm} arguments.
}
\usage{

\S4method{mzmatch}{numeric,numeric}(x, mz, mzdev=0, ppm=10)

\S4method{mzmatch}{matrix,numeric}(x, mz, mzdev=0, ppm=10)

}
\arguments{

  (in alphabetic order)

  \item{mz}{
    A numeric vector of M/Z values.
  }

  \item{mzdev}{
    A (single) numeric value specifying an allowed maximal
    shift/deviation of the M/Z values. See description of the return
    value for more details.
  }

  \item{ppm}{
    A (single) numeric value specifying the parts per million deviation
    of the M/Z value. See description of the return value for more
    details.
  }

  \item{x}{
    A numeric vector or matrix with the M/Z values that should be
    matched against \code{mz}. If a \code{matrix} it is expected to have
    two columns, specifying a lower and upper boundary for an M/Z
    range. See details below for more information.
  }

}
\value{
  The method always returns a \code{list} of length
  \code{length(x)}. Each list element is a two-column \code{matrix} with
  as many rows as there are matches of \code{x} in \code{mz}. If no
  match was found for an element in \code{x} the matrix will still have
  1 row, but contain only \code{NA} values. The columns of the
  \code{matrix} are \code{"idx"}, the index of the match in \code{mz},
  and \code{"deltaMz"}, the distance of the element in \code{x} to the
  matching \code{mz} value (i.e. the difference of the mass to the M/Z
  value of the compound).
}
\details{
  The method returns a match, if the distance between the masses is
  smaller than specified with the arguments \code{mzdev} and \code{ppm},
  i.e. if the distance \code{abs(x - mz)} is smaller or equal to
  \code{(mzdev + (x / 1000000) * ppm))}. Note that, if \code{x} is a
  \code{matrix} specifying an M/Z value range in each row, the midpoint
  of this range (i.e. the mean value) is matched against \code{mz} with
  the difference between that midpoint to the lower bound being added to
  the \code{mzdev} argument.
}
\author{
  Johannes Rainer
}
\examples{

## Simple example; matching numeric values.
gotMz <- c(2.23, 2.53)
## Define the mz argument
matchWith <- seq(1, 5, by=0.1)
matchWith
## Matching 'as.is' does not return any match
mzmatch(gotMz, matchWith)

## For this simple example we have increase the 'mzdev'; by setting
## mzdev=0.03 we get a match with the element 2.2 and 2.5
res <- mzmatch(gotMz, matchWith, mzdev=0.03)
res
res <- do.call(rbind, res)
matchWith[res[, "idx"]]

## If we increase the mzdev to 0.07 we get in addition a match to 2.3 and 2.6
res <- mzmatch(gotMz, matchWith, mzdev=0.07)
res

matchWith[res[[1]][, "idx"]]

## Alternatively, specify a M/Z range.
gotMz <- rbind(c(2.23, 2.25),
               c(2.53, 2.57))
## In this case, the midpoint of each range is matched against the mz values
## and the difference midpoint to min of range is added to the mzdev. In our
## example this means, we're matching 2.24 and 2.55 and adding 0.01 to the mzdev
## argument for the first, and 0.02 to the second. That way we get already two
## matches for the second mass.
res <- mzmatch(gotMz, matchWith, mzdev=0.03)
res

####
## Real life example: find a peak group matching a certain mass.

## We load the test data from the xcms package and perform the /default/ analysis
## as provided by the xcms vignette
library(faahKO)
xset <- group(faahko)
xset <- retcor(xset, family="symmetric")
xset <- group(xset, bw=10)
xset <- fillPeaks(xset)

## The groups method extracts the definition of the identified, and matched
## peaks
head(groups(xset))

## We can now use the matchmz method to find peak groups matching the peaks
## shown in Figure 3 of the xcms vignette.
gotMz <- rbind(c(300.1, 300.2),
               c(298.1, 298.2),
               c(491.1, 491.3))
idx <- mzmatch(gotMz, groups(xset)[, "mzmed"])

idx

groups(xset)[idx[[1]][, 1], ]

}
\keyword{data}
