####============================================================
##  grow
##
##  Supposed, x is a numeric of length 2 representing a range,
##  grow expands that range by by.
####------------------------------------------------------------
setMethod("grow", "numeric", function(x, by=NULL, ...){
    if(length(x) != 2)
        stop("Argument 'x' is supposed to be a numeric of length 2",
             " representing a range.")
    if(is.null(by))
        return(x)
    if(length(by) > 1){
        warning("Length of 'by' is > 1; using only the first element.")
        by <- by[1]
    }
    ## Sorting the range... to ensure it's really lower, upper bound.
    x <- sort(x)
    x[1] <- x[1] - by
    x[2] <- x[2] + by
    return(x)
})
setMethod("grow", "matrix", function(x, by=NULL, ...){
    ## x has to be a numeric matrix!
    if(!is.numeric(x[, 1]))
        stop("Argument 'x' is supposed to be a numeric matrix!")
    if(ncol(x) != 2)
        stop("Argument 'x' is supposed to be a numeric matrix with",
             " 2 columns, each row representing a range!")
    if(is.null(by))
        return(x)
    if(length(by) > 1){
        if(length(by) != nrow(x)){
            warning("Length of 'by' is not 1 and does not match the",
                    " number of rows in 'x'! Using only the first element.")
            by <- rep(by[1], nrow(x))
        }
    }else{
        by <- rep(by, nrow(x))
    }
    ## Expanding
    x[, 1] <- x[, 1] - by
    x[, 2] <- x[, 2] + by
    return(x)
})


####============================================================
##  mzmatch
##
##  Performs matching of M/Z values against M/Z values given in a numeric vector 'mz'.
##
##  Returns a list of matrices, one matrix for each element in x with one row in the
##  matrix with the index of the corresponding match in argument 'mz' and the distance
##  (in M/Z dimension) between the two. A match is returned if the distance between the
##  masses is smaller than specified with the arguments 'mzdev' and 'ppm' (i.e. if the
##  distance abs(x - mz) is <= (mzdev + (x / 1000000) * ppm))
####------------------------------------------------------------
setMethod("mzmatch", signature(x="numeric", mz="numeric"),
          function(x, mz, mzdev=0, ppm=10){
              return(.multiMatchMzNumeric(x=x, mz=mz, mzdev=mzdev, ppm=ppm))
          })
## Same, but for x being a matrix.
setMethod("mzmatch", signature(x="matrix", mz="numeric"),
          function(x, mz, mzdev=0, ppm=10){
              ## if ncol(x) > 2 -> check if it has mzmin and mzmax (e.g. the groups or peak matrix)
              if(ncol(x) > 2){
                  if(all(c("mzmin", "mzmax") %in% colnames(x))){
                      ## Subset the matrix...
                      x <- x[, c("mzmin", "mzmax")]
                  }else{
                      stop("If 'x' is a matrix it should have 2 columns specifying the upper and",
                           " lower bounds of the mz ranges!")
                  }
              }
              return(.multiMatchMzMatrix(x=x, mz=mz, mzdev=mzdev, ppm=ppm))
          })
## Other way round: x numeric, mz: matrix

####============================================================
## .singleMatchMz
##
## Low level matching function. This is intended to work ONLY for a SINGLE peak
## or whatever specified with 'x'!
## We're accepting only:
## x of length 1: it's supposed to be the M/Z value.
## x of length 2: in this case x is supposed to be a range (start - end).
## mz: a numeric vector with the "reference" mz.
## mzdev: acceptable deviation of the mass, as absolute value.
## ppm: part per million.
## The function returns a matrix with columns idx and dist specifying the index of
## the match in argument mz and the distance to it. The matrix has as many rows as
## there are matches. If no match is found a matrix with a single row containing
## values NA for index and distance is returned.
####------------------------------------------------------------
.singleMatchMz <- function(x, mz, mzdev=0, ppm=10){
    ## FIXFIXFIX>> we're adding some tiny ppm to the ppm to avoid unexpected results
    ## by the comparison <= since abs(1-0.95) > 0.05 returns TRUE; although we expect
    ## 0.05 being identical to 0.05
    ## fixi <- .Machine$double.eps
    fixi <- 1e-9
    ## <<
    if(length(x) > 2)
        stop("'x' is supposed to be a single M/Z value, or a numeric of length",
             " 2 specifying an M/Z range.")
    ## Now, if x is a range:
    if(length(x) == 2){
        newX <- mean(x)
        ## Increase the mzdev by the difference mean to min.
        mzdev <- mzdev + abs(x[1] - newX)
        x <- newX
    }
    if(missing(mz))
        stop("Argument 'mz' is missing!")
    ## Done. Now calculate the distance of this x to all mz.
    dists <- abs(x - mz)
    ppm <- ppm / 1000000
    ## Now check which dists are "close enough"; replace all others with NA.
    closeIdx <- which(dists <= (mzdev + x*ppm + fixi))
    if(length(closeIdx) == 0){
        return(cbind(idx=NA, deltaMz=NA))
    }
    ## Else, reorder by distance.
    if(length(closeIdx) > 1){
        idx <- order(dists[closeIdx])
        closeIdx <- closeIdx[idx]
    }
    return(cbind(idx=closeIdx, deltaMz=dists[closeIdx]))
}

####============================================================
##  .multiMatchMz
##
##  Applies the .singleMatchMz to each elements of x.
####------------------------------------------------------------
.multiMatchMzNumeric <- function(x, mz, mzdev=0, ppm=10){
    ## Running the stuff with lapply
    return(lapply(x, function(z, themz){
        return(.singleMatchMz(z, mz=themz, mzdev=mzdev, ppm=ppm))
    }, themz=mz))
}
.multiMatchMzMatrix <- function(x, mz, mzdev=0, ppm=10){
    ## Using split lapply; it's slightly faster than an apply on the
    ## matrix.
    return(lapply(split(x, f=1:nrow(x)), function(z, themz){
          return(.singleMatchMz(z, mz=themz, mzdev=mzdev, ppm=ppm))
      }, themz=mz))
}



####============================================================
##  mass2adductmz
##
##  Calculates all possible adducts (based on data from the Fiehn lab) for the provided
##  masses.
##  ionAdduct: specify the name of the ion adduct; returns x if ionAdduct is NULL.
##  Returns a list of adduct mz values, names of the list representing the adduct, each
##  element being a numeric vector (same length than x) with the mz values for that adduct.
####------------------------------------------------------------
setMethod("mass2adductmz", "numeric",
          function(x, ionAdduct=supportedIonAdducts()){
              if(missing(ionAdduct))
                  ionAdduct <- supportedIonAdducts()
              return(.mass2adductmz(x, ionAdduct=ionAdduct))
          })
####============================================================
##  adductmz2mass
##
##  Calculates mass for all possible adducts (based on data
##  from the Fiehn lab) for the provided mz.
##  Returns a list of mass values, names of the list representing the adduct, each
##  element being a numeric vector (same length than x) with the mass if the provided
##  mz was this adduct ion.
####------------------------------------------------------------
setMethod("adductmz2mass", "numeric",
          function(x, ionAdduct=supportedIonAdducts()){
              if(missing(ionAdduct))
                  ionAdduct <- supportedIonAdducts()
              return(.adductmz2mass(x, ionAdduct=ionAdduct))
          })



####============================================================
##  supporteIonAdducts
##
##  Returns all supported ion names for adduct ion names.
setMethod("supportedIonAdducts", signature(x="missing", charge="missing"),
          function(x, charge){
              return(supportedIonAdducts(charge="both"))
          })
setMethod("supportedIonAdducts", signature(x="missing", charge="character"),
          function(x, charge="both"){
              charge <- match.arg(charge, c("both", "pos", "neg"))
              if(charge=="both")
                  return(ESIADDUCTS$ion_name)
              if(charge=="pos")
                  return(ESIADDUCTS$ion_name[ESIADDUCTS$charge > 0])
              if(charge=="neg")
                  return(ESIADDUCTS$ion_name[ESIADDUCTS$charge < 0])
          })

####============================================================
##  .mass2adductmz
##
##  internal function to get the adduct(s) for a given mass.
##  ionAdduct is just to optionally subset to specific ion names.
##
##  returns a list (always) the elements being the calculated adduct for one
##  specific ionAdduct.
####------------------------------------------------------------
.mass2adductmz <- function(x, ionAdduct=NULL){
    if(is.null(ionAdduct))
        return(list(M=x))
    ## Check if ionAdduct fit.
    notPresent <- !(ionAdduct %in% ESIADDUCTS$ion_name)
    if(all(notPresent))
        stop("None of the provided ion names are valid!")
    if(any(notPresent)){
        warning("Ion adduct names ", paste(ionAdduct[notPresent], collapse=", "),
                " are not valid",
                " and have been dropped.")
        ionAdduct <- ionAdduct[!notPresent]
    }
    tmp <- ESIADDUCTS[ESIADDUCTS$ion_name %in% ionAdduct, ]
    tmp <- split(tmp, tmp$ion_name)
    return(lapply(tmp, function(z, x){
        return(x * z[1, "mult"] + z[1, "mass"])
    }, x=x))
}
## And the reverse...
.adductmz2mass <- function(x, ionAdduct=NULL){
    if(is.null(ionAdduct))
        return(list(M=x))
    ## Check if ionAdduct fit.
    notPresent <- !(ionAdduct %in% ESIADDUCTS$ion_name)
    if(all(notPresent))
        stop("None of the provided ion names are valid!")
    if(any(notPresent)){
        warning("Ion adduct names ", paste(ionAdduct[notPresent], collapse=", "),
                " are not valid",
                " and have been dropped.")
        ionAdduct <- ionAdduct[!notPresent]
    }
    tmp <- ESIADDUCTS[ESIADDUCTS$ion_name %in% ionAdduct, ]
    tmp <- split(tmp, tmp$ion_name)
    return(lapply(tmp, function(z, x){
        return(x * z[1, "mult"] - z[1, "mass"])
    }, x=x))
}
