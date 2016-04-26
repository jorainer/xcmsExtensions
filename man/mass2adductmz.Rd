\name{mass2adductmz}
\Rdversion{1.1}
\alias{mass2adductmz}
\alias{adductmz2mass}
\alias{mass2adductmz,numeric-method}
\alias{adductmz2mass,numeric-method}
\alias{supportedIonAdducts}
\alias{supportedIonAdducts,missing,missing-method}
\alias{supportedIonAdducts,missing,character-method}

\title{
  Conversions between real mass and M/Z of possible ion adducts
}
\description{
  The \code{mass2adductmz} and \code{adductmz2mass} enable the
  conversion between the mass of a compound and the M/Z value for a
  specific, or all possible, ion adducts for that compound. For
  conversion the data provided by the Fiehn lab (in their ESI-MS-adduct
  calculator) is used.
}
\usage{

\S4method{adductmz2mass}{numeric}(x, ionAdduct=supportedIonAdducts())

\S4method{mass2adductmz}{numeric}(x, ionAdduct=supportedIonAdducts())

\S4method{supportedIonAdducts}{missing,character}(x, charge="both")

}
\arguments{

  (in alphabetic order)

  \item{charge}{
    For \code{supportedIonAdducts}: either \code{"pos"}, \code{"neg"} or
    \code{"both"} (the default) to return the names of the supported
    positive, negative, or all ion adducts, respectively.
  }

  \item{ionAdduct}{
    A character string or vector representing the name of an ion
    adduct. By default, the most commonly found adducts in ESI are
    used. See \code{supportedIonAdducts} for all supported ion adduct
    names. If \code{NULL} is used, the compound is supposed to be
    already an ion and thus no conversion is performed.
  }

  \item{x}{
    For \code{adductmz2mass} and \code{mass2adductmz}: A numeric
    representing either the mass of a compound (for
    \code{mass2adductmz}) or the M/Z (for \code{adductmz2mass}) to be
    converted.
  }

}
\section{Methods}{

  \describe{
    \item{adductmz2mass}{
      Converts the M/Z value to the mass of the compound, assuming that
      the measured M/Z value corresponds to one of the most commonly seen
      ion adducts.  The name of the ion adduct (or of several
      adducts) can be specified with the \code{ionAdduct} argument.
      By default, the M/Z for the most common adducts are used (those
      returned by the \code{supportedIonAdducts} method). If
      \code{ionAdduct} is \code{NULL} the input values are returned,
      i.e. it is assumed that the mass corresponds already to an ion

      The method returns a list of masses, with the names of the list
      representing the ion adduct and the elements the mass for these ion
      adducts.
    }

    \item{mass2adductmz}{
      Converts the mass of a compound to M/Z values for the most commonly
      found ion adducts. The name of the ion adduct (or of several
      adducts) can be specified with the \code{ionAdduct} argument. By
      default, the M/Z for the most common adducts are used (those
      returned by the \code{supportedIonAdducts} method). If
      \code{ionAdduct} is \code{NULL} the input values are returned,
      i.e. it is assumed that the mass corresponds already to an ion.

      The method returns a list of M/Z values, the names corresponding to
      the name of the ion adduct and the elements being the M/Z values for
      these.
    }

    \item{supportedIonAdducts}{
      Returns a character vector with the names of the supported ion
      adducts.
    }

  }
}

\value{
  See method description above.
}
\details{
  The mass and adduct information are taken from the ESI MS adduct
  calculator from the Fiehn lab [Huang 1999].
}
\author{
  Johannes Rainer
}
\references{
  Huang N, Siegel MM, Kruppa GH & Laukien FH (1999) Automation of a
  Fourier transform ion cyclotron resonance mass spectrometer for
  acquisition, analysis, and e-mailing of high-resolution exact-mass
  electrospray ionization mass spectral data. \emph{Journal of the
    American Society for Mass Spectrometry}. pp1166-1173.

}
\examples{

## Get the names of all supported ion adducts
supportedIonAdducts()

## Get the names of the negative ion adducts
supportedIonAdducts(charge="neg")

## Specify some M/Z values.
mzs <- c(169.2870, 169.5650, 346.4605)

## Get the masses for these M/Z values assuming that they are from negative ion
## adducts of the real compound.
negMass <- adductmz2mass(mzs, ionAdduct=supportedIonAdducts(charge="neg"))

## The mass assuming that the M/Z corresponds to the M-H adduct.
negMass[["M-H"]]

## Calculating the mass form the 'M-H' adduct.
mass2adductmz(negMass[["M-H"]], ionAdduct="M-H")

}
\keyword{data}
