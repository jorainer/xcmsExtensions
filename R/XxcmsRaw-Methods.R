## Methods for xcmsRaw classes.
####============================================================
##  msSlice
##
##  Extract an MSslice object from the xcmsRaw object. This MSslice
##  contains (one or more) MSdata objects.
setMethod("msSlice", "xcmsRaw",
          function(object, mzrange=NULL, rtrange=NULL, ...){
              call <- match.call()
              ## Get the data from the raw object.
              ## We'll call the msData
              res <- msData(object, mzrange=mzrange, rtrange=rtrange, ...)
              if(is(res, "MSdata")){
                  return(msSlice(res))
              }else{
                  return(MSsliceList(lapply(res, MSslice)))
                  ## Got a list of objects, return a MSsliceList.
              }
          })

####============================================================
##  msMap
##
##  Gets either the full or a subset of the profile matrix as
##  MSmap object (defined in MSnbase), i.e. a two dimensional
##  matrix m/z against retention time, the values being the intensities
##  measured for each mz rt tuple.
####------------------------------------------------------------
setMethod("msMap", "xcmsRaw",
          function(object, mzrange=NULL, rtrange=NULL, intrange=NULL,
                   resMz=0.005, zeroIsNA=FALSE, ...){
              call <- match.call()
              datmat <- getData(x=object, mzrange=mzrange, rtrange=rtrange,
                             intrange=intrange)
              if(is(datmat, "list")){
                  res <- lapply(datmat, function(z){
                      return(.matrix2msmap(z, resMz=resMz, zeroIsNA=zeroIsNA,
                                           call=call))
                      })
              }else{
                  res <- .matrix2msmap(datmat, resMz=resMz, zeroIsNA=zeroIsNA,
                                       call=call)
              }
              return(res)
})
## The internal working fun... assumes that x is a matrix with
## columns time, mz and intensity such as the one returned by getData.
.matrix2msmap <- function(x, resMz=0.005, zeroIsNA=TRUE, call=match.call()){
    mzr <- range(x[, "mz"])
    minM <- round(mzr[1]/resMz) * resMz
    maxM <- round(mzr[2]/resMz) * resMz
    numC <- (maxM - minM) / resMz + 1
    ## Calculating scnIdx, i.e. the index (0-based) in which the retention
    ## time changes.
    rl <- runLength(Rle(x[, "time"]))
    scnIdx <- c(0, cumsum(rl)[-length(rl)])
    mat <- xcms:::profBinM(x[, "mz"], x[, "intensity"], scnIdx, numC,
                           minM, maxM, FALSE, list())
    msm <- new("MSmap", call=call, map=t(mat), rt=unique(x[, "time"]),
               mz=seq(minM, maxM, length.out=numC), t=FALSE)
    return(msm)
}

####============================================================
##  msData
##
##  Extract MS data from an xcmsRaw object for one or more specified
##  rt and/or mz ranges
####------------------------------------------------------------
setMethod("msData", "xcmsRaw",
          function(object, mzrange=NULL, rtrange=NULL, ...){
              ## Do input argument checking:
              ## mzrange and rtrange can be numeric of length 2.
              haverows <- unique(c(nrow(mzrange), nrow(rtrange)))
              call <- match.call()
              gotMatrix <- FALSE
              if(length(haverows) > 1)
                  stop("If 'mzrange' and 'rtrange' are matrices they have to",
                       " have the same number of rows!")
              if(is.null(mzrange) & is.null(rtrange)){
                  message("Both 'mzrange' and 'rtrange' are NULL, returning the full data.")
              }else{
                  if(!is.null(mzrange)){
                      if(is(mzrange, "matrix")){
                          ## Then rtrange has to NULL or a matrix
                          if(!is.null(rtrange)){
                              if(!is(rtrange, "matrix"))
                                  stop("If 'mzrange' is a matrix, also 'rtrange', if not NULL,",
                                       " has to be a matrix!")
                          }
                          gotMatrix <- TRUE
                          if(ncol(mzrange) != 2){
                              stop("'mzrange' is supposed to be either a numeric of",
                                   " length 2 or a matrix with two columns!")
                          }
                      }
                  }
                  if(!is.null(rtrange)){
                      if(is(rtrange, "matrix")){
                          ## Then rtrange has to NULL or a matrix
                          if(!is.null(mzrange)){
                              if(!is(mzrange, "matrix"))
                                  stop("If 'rtrange' is a matrix, also 'mzrange', if not NULL,",
                                       " has to be a matrix!")
                          }
                          gotMatrix <- TRUE
                          if(ncol(rtrange) != 2){
                              stop("'rtrange' is supposed to be either a numeric of",
                                   " length 2 or a matrix with two columns!")
                          }
                      }
                  }
              }
              ## Done, now we can proceed to the function.
              if(!gotMatrix){
                  ## Just get the data.
                  res <- .doGetMSdata(object, mzrange=mzrange, rtrange=rtrange, ...)
                  return(res)
              }else{
                  ## A little more to do...
                  ## We want to split the ranges
                  ## into a list (of lists) and run a lapply on that list element.
                  ## Note, funnily enough I can call NULL[1, ] and get NULL.
                  arglist <- vector("list", length=haverows)
                  sct <- scantimes(object)
                  for(i in 1:haverows){
                      arglist[[i]] <- list(mzrange=mzrange[i, ],
                                         rtrange=rtrange[i, ])
                  }
                  res <- lapply(arglist, function(z){
                      return(.doGetMSdata(object, mzrange=z$mzrange,
                                          rtrange=z$rtrange, sct=sct))
                  })
                  return(res)
              }
          })
## The workhorse.
.doGetMSdata <- function(x, sct=NULL, mzrange=NULL, rtrange=NULL, intrange=NULL){
    if(is.null(sct))
        sct <- scantimes(x)
    ## Return the full blown matrix if no subsetting is done.
    if(is.null(mzrange) & is.null(rtrange) & is.null(intrange)){
        return(MSdata(mz=x@env$mz, rtime=sct, intensity=x@env$intensity))
    }else{
        returnMe <- -9
        ## Subset based on mzrange
        if(!is.null(mzrange)){
            if(length(mzrange) != 2)
                stop("'mzrange' has to be a numeric vector of length 2.")
            mzrange <- sort(mzrange)
            returnMe <- which(x@env$mz >= mzrange[1] & x@env$mz <= mzrange[2])
            ## Return an empty matrix if we can not find any of these.
            if(length(returnMe) == 0)
                return(MSdata(rtime=numeric(), mz=numeric(), intensity=numeric()))
        }
        ## Subset based on rtrange
        if(!is.null(rtrange)){
            if(length(rtrange) != 2)
                stop("'rtrange' has to be a numeric vector of length 2.")
            rtrange <- sort(rtrange)
            returnMeRt <- which(sct >= rtrange[1] & sct <= rtrange[2])
            if(length(returnMeRt) == 0)
                return(MSdata(rtime=numeric(), mz=numeric(), intensity=numeric()))
            if(returnMe[1] > 0){
                returnMe <- intersect(returnMe, returnMeRt)
            }else{
                returnMe <- returnMeRt
            }
            ## Return an empty matrix if we can not find any of these.
            if(length(returnMe) == 0)
                return(MSdata(rtime=numeric(), mz=numeric(), intensity=numeric()))
        }
        ## Subset based on intrange
        if(!is.null(intrange)){
            if(length(intrange) != 2)
                stop("'intrange' has to be a numeric vector of length 2.")
            intrange <- sort(intrange)
            returnMeRt <- which(x@env$intensity >= intrange[1] & x@env$intensity <= intrange[2])
            if(length(returnMeRt) == 0)
                return(MSdata(rtime=numeric(), mz=numeric(), intensity=numeric()))
            if(returnMe[1] > 0){
                returnMe <- intersect(returnMe, returnMeRt)
            }else{
                returnMe <- returnMeRt
            }
        }
        return(MSdata(mz=x@env$mz[returnMe],
                      rtime=sct[returnMe],
                      intensity=x@env$intensity[returnMe]))
    }
}


####============================================================
##  getData
##
##  Get a matrix with the raw data, i.e. retention time, mz and
##  intensity. Allows to optionally define subsets of the data.
setMethod("getData", signature(x="xcmsRaw", mzrange="NumericOrNull",
                               rtrange="NumericOrNull",
                               intrange="NumericOrNull"),
          function(x, mzrange=NULL, rtrange=NULL, intrange=NULL){
              return(.doGetData(x, mzrange=mzrange, rtrange=rtrange,
                                intrange=intrange))
          })
## The same method but with multiple ranges.
setMethod("getData", signature(x="xcmsRaw", mzrange="MatrixOrNull",
                               rtrange="MatrixOrNull",
                               intrange="MatrixOrNull"),
          function(x, mzrange=NULL, rtrange=NULL, intrange=NULL){
              if(!is.null(mzrange) | !is.null(rtrange) | !is.null(intrange)){
                  ## Check input arguments:
                  haverows <- unique(c(nrow(mzrange), nrow(rtrange), nrow(intrange)))
                  if(length(haverows) > 1)
                      stop("If specified, 'mzrange', 'rtrange' and 'intrange' have to",
                           " have the same number of rows!")
                  if(!is.null(mzrange)){
                      if(ncol(mzrange) != 2)
                          stop("'mzrange' is supposed to be either a numeric of",
                               " length 2 or a matrix with two columns!")
                  }
                  if(!is.null(rtrange)){
                      if(ncol(rtrange) != 2)
                          stop("'rtrange' is supposed to be either a numeric of",
                               " length 2 or a matrix with two columns!")
                  }
                  if(!is.null(intrange)){
                      if(ncol(intrange) != 2)
                          stop("'intrange' is supposed to be either a numeric of",
                               " length 2 or a matrix with two columns!")
                  }
                  ## Finally done; now we can proceed. We want to split the ranges
                  ## into a list (of lists) and run a lapply on that list element.
                  ## Note, funnily enough I can call NULL[1, ] and get NULL.
                  arglist <- vector("list", length=haverows)
                  for(i in 1:haverows){
                      arglist[[i]] <- list(mzrange=mzrange[i, ],
                                         rtrange=rtrange[i, ],
                                         intrange=intrange[i, ])
                  }
                  res <- lapply(arglist, function(z){
                      return(.doGetData(x, mzrange=z$mzrange,
                                        rtrange=z$rtrange,
                                        intrange=z$intrange))
                  })
                  return(res)
              }else{
                  ## Well, all arguments are NULL, just call the function.
                  return(.doGetData(x, mzrange=mzrange, rtrange=rtrange,
                                    intrange=intrange))
              }
          })
## The workhorse.
.doGetData <- function(x, sct=NULL, mzrange=NULL, rtrange=NULL, intrange=NULL){
    if(is.null(sct))
        sct <- scantimes(x)
    ## Return the full blown matrix if no subsetting is done.
    if(is.null(mzrange) & is.null(rtrange) & is.null(intrange)){
        return(cbind(time=sct, mz=x@env$mz, intensity=x@env$intensity))
    }else{
        returnMe <- -9
        ## Subset based on mzrange
        if(!is.null(mzrange)){
            if(length(mzrange) != 2)
                stop("'mzrange' has to be a numeric vector of length 2.")
            mzrange <- sort(mzrange)
            returnMe <- which(x@env$mz >= mzrange[1] & x@env$mz <= mzrange[2])
            ## Return an empty matrix if we can not find any of these.
            if(length(returnMe) == 0)
                return(cbind(time=numeric(), mz=numeric(), intensity=numeric()))
        }
        ## Subset based on rtrange
        if(!is.null(rtrange)){
            if(length(rtrange) != 2)
                stop("'rtrange' has to be a numeric vector of length 2.")
            rtrange <- sort(rtrange)
            returnMeRt <- which(sct >= rtrange[1] & sct <= rtrange[2])
            if(length(returnMeRt) == 0)
                return(cbind(time=numeric(), mz=numeric(), intensity=numeric()))
            if(returnMe[1] > 0){
                returnMe <- intersect(returnMe, returnMeRt)
            }else{
                returnMe <- returnMeRt
            }
            ## Return an empty matrix if we can not find any of these.
            if(length(returnMe) == 0)
                return(cbind(time=numeric(), mz=numeric(), intensity=numeric()))
        }
        ## Subset based on intrange
        if(!is.null(intrange)){
            if(length(intrange) != 2)
                stop("'intrange' has to be a numeric vector of length 2.")
            intrange <- sort(intrange)
            returnMeRt <- which(x@env$intensity >= intrange[1] & x@env$intensity <= intrange[2])
            if(length(returnMeRt) == 0)
                return(cbind(time=numeric(), mz=numeric(), intensity=numeric()))
            if(returnMe[1] > 0){
                returnMe <- intersect(returnMe, returnMeRt)
            }else{
                returnMe <- returnMeRt
            }
        }
        return(cbind(time=sct[returnMe],
                     mz=x@env$mz[returnMe],
                     intensity=x@env$intensity[returnMe]))
    }
}

####============================================================
##  scantimes
##
##  Returns a vector with the retention time for the individual
##  measured mz values.
####------------------------------------------------------------
setMethod("scantimes", "xcmsRaw", function(x){
    return(rep.int(x@scantime, times=diff(c(x@scanindex, length(x@env$mz)))))
})

