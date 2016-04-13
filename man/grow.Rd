\name{grow}
\alias{grow}
\alias{grow,numeric-method}
\alias{grow,matrix-method}
\title{
  Expand/grow a numeric range
}
\description{
  The \code{grow} method allows to expand a numeric range by a given
  number, i.e. the number specified with the argument \code{by} is
  subtracted from the lower bound and added to the upper bound, such
  that the range is grown by \code{2 * by}.
}
\usage{

\S4method{grow}{numeric}(x, by=NULL, ...)

\S4method{grow}{matrix}(x, by=NULL, ...)

}
\arguments{
  \item{by}{
    Numeric vector of length 1 (or with a length matching the number of
    rows of \code{x} if \code{x} is a \code{matrix}).
  }

  \item{x}{
    Either a numeric vector of length 2 specifying the lower and upper
    bound of the range or a numeric matrix with two columns each row
    specifying a range.
  }

  \item{...}{
    Presently not used.
  }
}
\value{
  Returns the numeric vector or matrix with expanded ranges.
}
\details{
  The function subtracts \code{by} from the first element of \code{x}
  and adds it to the second.
}
\author{
  Johannes Rainer
}
\examples{

## Expand a numeric range.
grow(c(3, 4), by=4)

## Define a numeric matrix of ranges and expand them.
mat <- rbind(c(3, 4),
             c(11, 15),
             c(1, 4))
grow(mat, by=3)

## Specify a by for each row
grow(mat, by=c(1, 2, 3))

}
\keyword{data}
