## Methods for xcmsSet classes.

####============================================================
##  peakGroupSummary
##
##  Calculate summaries for peak groups per sample class.
####------------------------------------------------------------
setMethod("peakGroupSummary", "xcmsSet", function(object, class, value="into"){
    if(missing(class)){
        ## Using internal class
        class <- sampclass(object)
    }else{
        ## Check class.
        if(length(class) != length(sampnames(object)))
            stop("Length of 'class' has to match the number of data files.")
        if(!is.factor(class))
            class <- factor(class)
    }
    ## Internal function to perform the summary.
    sumfun <- function(mat){
        res <- apply(mat, MARGIN=1, function(z){
            zs <- z[!is.na(z)]
            mz <- mean(zs)
            return(c(minInt=min(zs),
                     maxInt=max(zs),
                     meanInt=mz,
                     RSD=sd(zs)/mz,
                     propPresent=length(zs)/length(z)))
        })
        return(t(res))
    }
    ## OK, now get the group data.
    gv <- groupval(object, value=value)
    ## Do the call separately for each level.
    Res <- do.call(cbind, lapply(levels(class), function(z){
        return(sumfun(gv[, class==z, drop=FALSE]))
    }))
    colnames(Res) <- paste(rep(levels(class), each=5),
                           c("minInt", "maxInt", "meanInt", "RSD",
                             "propPresent"), sep=".")
    return(Res)
})

####============================================================
##  msSlice
##
##  Extract an MSslice (or MSsliceList) object from the xcmsSet object.
##  Each MSslice object containing an MSdata object for each sample in the
##  xcmsSet object. NOTE: can NOT pass rt down; even BPPARAM does not work!
setMethod("msSlice", "xcmsSet",
          function(object, mzrange=NULL, rtrange=NULL, rt="corrected", BPPARAM=MulticoreParam()){
              rt <- match.arg(rt, c("corrected", "raw"))
              call <- match.call()
              ## Will run the msSlice separately on each xcmsRaw object. To reduce the required
              ## memory demand,
              ## use bplapply.
              ## We're setting profstep to 0, so we don't do the profile matrix calculation.
              sampidx <- 1:length(filepaths(object))
              resu <- bplapply(as.list(sampidx), function(z, object, rt, mzrange, rtrange){
                  ## Get the raw object. We don't want the profile matrix to be generated
                  ## thus we set profstep to 0.
                  xr <- getXcmsRaw(object, sampleidx=z, profstep=0, rt=rt)
                  ## Get the MSsliceList on that.
                  slices <- msSlice(xr, mzrange=mzrange, rtrange=rtrange)
                  return(slices)
              }, object=object, rt=rt, mzrange=mzrange, rtrange=rtrange, BPPARAM=BPPARAM)
              ## OK, resu will now be a list of MSslice or MSsliceList objects (depending
              ## on whether mzrange or rtrange were a matrix or not).
              ## Will now combine the MSdata for the individual samples into one MSslice.
              if(is(resu[[1]], "MSslice")){
                  ## Just ensure that all are MSslices...
                  if(!all(unlist(lapply(resu, function(z){is(z, "MSslice")}))))
                      stop("Not all of the returned objects are 'MSslice' objects!")
                  return(msSlice(object=unlist(lapply(resu, msData)),
                                 phenoData=AnnotatedDataFrame(phenoData(object))))
              }
              if(is(resu[[1]], "MSsliceList")){
                  if(!all(unlist(lapply(resu, function(z){is(z, "MSsliceList")}))))
                      stop("Not all of the returned objects are 'MSsliceList' objects!")
                  ## Check if the lengths are the same.
                  Ls <- unique(unlist(lapply(resu, length)))
                  if(length(Ls) != 1)
                      stop("Something wrong: got different lengths of MSsliceList objects!")
                  ## Do a for loop extracting each time the MSdata i and storing that into a
                  ## new MSslice
                  mslList <- vector("list", Ls)
                  for(i in 1:Ls){
                      mslList[[i]] <- msSlice(object=unlist(lapply(resu, function(z){
                          return(msData(z@slices[[i]]))
                      })), phenoData=AnnotatedDataFrame(phenoData(object)))
                  }
                  return(MSsliceList(mslList))
              }else{
                  stop("The method returned something wrong!")
              }

          })





## ## Extension methods and functions for the xcms package.
## ## Simple function to find the closest peak group.
## ## x: xcmsSet object.
## ## mz: the m/z to which the closest peaks should be returned.
## ## rt: the retention time to which the closest peaks should be returned.
## ## mzmax: the maximal allowed distance mz to mz of the peak.
## ## rtmax: the maximal allowed distance rt to rt of the peak.
## ## nres: the maximal number of peaks that should be returned. If Inf, it returns all.
## .findClosestGroup <- function(x, mz=NULL, rt=NULL, mzdev=NA, rtdev=NA, nres=1){
##     if(is.null(mz) & is.null(rt))
##         stop("Either 'mz' or 'rt' or both have to be specified!")
##     if(nrow(groups(x)) == 0)
##         stop("No groups found. You should run the 'group' command first!")
##     if(length(mz) > 1 | length(rt) > 1)
##         stop("mz or rt should be of length 1!")
##     ## Difference to mz
##     outsideRange <- numeric()
##     if(!is.null(mz)){
##         diffMz <- groups(x)[, "mzmed"] - mz
##         outsideRange <- which(abs(diffMz) > mzdev)
##     }
##     ## Difference to rt
##     if(!is.null(rt)){
##         diffRt <- groups(x)[, "rtmed"] - rt
##         outsideRange <- sort(c(outsideRange, which(abs(diffRt) > rtdev)))
##     }
##     ## Use Euclidian's distance if both are specified.
##     if(!is.null(mz) & !is.null(rt)){
##         theDist <- sqrt(diffMz^2 + diffRt^2)
##     }else{
##         ## Otherwise, just use the absolute distance.
##         if(!is.null(mz)){
##             theDist <- abs(diffMz)
##         }else{
##             theDist <- abs(diffRt)
##         }
##     }
##     ## "flag" those that are outside the range
##     if(length(outsideRange) > 0)
##         theDist[outsideRange] <- NA
##     idx <- order(theDist)
##     res <- data.frame(index=idx, distance=theDist[idx])
##     ## Only return nres elements.
##     res <- res[1:(min(c(nres, nrow(res)))), , drop=FALSE]
##     res <- res[!is.na(res$distance), , drop=FALSE]
##     return(res)
## }
## ####============================================================
## ##  closestGroup
## ##
## ##  Method to find and return indices of peak-groups that are closest
## ##  to the specified mz and or rt range.
## ##  x: xcmsSet object.
## ##  mz: the m/z to which the closest peaks should be returned.
## ##  rt: the retention time to which the closest peaks should be returned.
## ##  mzdev: the maximal allowed distance mz to mz of the peak.
## ##  rtdev: the maximal allowed distance rt to rt of the peak.
## ##  nres: the maximal number of peaks that should be returned. Set to Inf to return
## ##        all.
## ##  The method returns a list of data.frames, one list element for each specified
## ##  mz (or rt). The data.frame has columns "index" and "distance" indicating the
## ##  index of the group in the groups (or groupval) matrix and the distance of the
## ##  actual peak to the user provided coordinates. In case both mz and rt are specified
## ##  the distance represents the Euclidian distance between the points in the mz and rt
## ##  space.
## ##  sort: "closest": return results ordered by distance to specified region.
## ##        "intensity": return results (within range) ordered by signal intensity ("into").
## ##        "npeaks": return results (within range) ordered by number of peaks.
## ####------------------------------------------------------------
## setMethod("closestGroup", "xcmsSet", function(x, mz=NULL, rt=NULL,
##                                               mzdev=NA, rtdev=NA,
##                                               nres=1){
##     ## Allow length of mz, rt > 1.
##     if(is.null(mz) & is.null(rt))
##         stop("Either 'mz' or 'rt' (or both) have to be specified!")
##     if(!is.null(mz) & !is.null(rt)){
##         if(length(mz) != length(rt))
##             stop("If 'mz' and 'rt' are specified, their length has to match!")
##     }
##     ## Define the names for the result list.
##     mzn <- "mz:-"
##     rtn <- "rt:-"
##     if(!is.null(mz)){
##         mzn <- paste0("mz:", format(mz, digits=7))
##     }
##     if(!is.null(rt))
##         rtn <- paste0("rt:", format(rt, digits=7))
##     listNames <- paste0(mzn, ";", rtn, ";")
##     ## Prepare the actual call.
##     if(is.null(mz))
##         mz <- rep(NA, length(rt))
##     if(is.null(rt))
##         rt <- rep(NA, length(mz))
##     toProcess <- split(data.frame(mz=mz, rt=rt), f=1:length(mz))
##     res <- lapply(toProcess, function(z){
##         ## Calling the function.
##         localMz <- NULL
##         localRt <- NULL
##         if(!is.na(z[1, 1]))
##             localMz <- z[1, 1]
##         if(!is.na(z[1, 2]))
##             localRt <- z[1, 2]
##         return(.findClosestGroup(x=x, mz=localMz, rt=localRt,
##                                  mzdev=mzdev, rtdev=rtdev,
##                                  nres=nres))
##     })
##     ## Check if we could use the names.
##     if(length(unique(listNames)) == length(res))
##         names(res) <- listNames
##     return(res)
## })

## #### Test this
## allmz <- c(2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4)
## mzmat <- matrix(allmz, ncol=1)
## colnames(mzmat) <- "mz"
## library(RUnit)
## checkEquals(unname(.singlematchmz(mzmat, 2.001)["pos"]), 1)
## ## add intensity.
## mzmat <- cbind(mzmat, intensity=c(2, 3, 4, 5, 6, 2, 3, 4, 5, 8, 4, 3))
## checkEquals(unname(.singlematchmz(mzmat, 2.001)["pos"]), 10)
## ## ranges.
## mzmat <- cbind(mzmat, mzmin=mzmat[, "mz"]-1, mzmax=mzmat[, "mz"]+1)
## checkEquals(unname(.singlematchmz(mzmat, 2.001)["pos"]), 10)

## ## .matchmz
## checkEquals(unname(.matchmz(allmz, 2)[, "pos"]), 1)
## checkEquals(unname(.matchmz(allmz, c(2, 3))[, "pos"]), c(1, 2))
## checkEquals(unname(.matchmz(mzmat, c(2, 3))[, "pos"]), c(10, 11))
## checkEquals(unname(.matchmz(mzmat, c(3, 2))[, "pos"]), c(11, 10))



## ####============================================================
## ####============================================================
## ##  Test stuff
## ##
## ####------------------------------------------------------------
## spikes <- read.table("data/txt/_input-spiked-compounds.txt", sep="\t",
##                      as.is=TRUE, header=TRUE)


## ## testing stuff...
## which <- 2
## ## Get the peak group closest to the specified mass.
## spkGroup <- closestGroup(xsRt, mz=spikes[which, "m"], nres=10)
## spkGroup

## groupIdx <- spkGroup[[1]]$index[1]
## ## Now we're getting the raw file for the first control
## xraw <- getXcmsRaw(xsRt, 6)

## ## get the EIC:
## eics <- getEIC(xsRt, groupidx=groupIdx)

## par(mfrow=c(1, 2))
## plot(eics, col=paste0(groupCols[xsRt$group], "80"), groupidx=1)
## legend("topright", col=groupCols, lty=1, legend=names(groupCols))
## ## Now plot the rt range.
## abline(v=groups(xsRt)[groupIdx,
##                       c("rtmin", "rtmax")], col="darkgrey")

## groupval(xsRt, value="into")[groupIdx, ]
## groups(xsRt)[groupIdx, ]

## ## not so bad, but still...
## eicsMz <- getEIC(xsRt, mzrange=groups(xsRt)[groupIdx, c("mzmin", "mzmax"), drop=FALSE])
## plot(eicsMz, col=paste0(groupCols[xsRt$group], "80"), groupidx=1)
## legend("topright", col=groupCols, lty=1, legend=names(groupCols))
## ## Now plot the rt range.
## abline(v=groups(xsRt)[groupIdx,
##                       c("rtmin", "rtmax")], col="darkgrey")

## ##
## ## Get the raw matrix...
## rtr <- groups(xsRt)[groupIdx, c("rtmin", "rtmax"), drop=FALSE]
## mzr <- groups(xsRt)[groupIdx, c("mzmin", "mzmax"), drop=FALSE]

## rm <- rawMat(xraw,
##              mzrange=groups(xsRt)[groupIdx, c("mzmin", "mzmax"), drop=FALSE],
##              rtrange=groups(xsRt)[groupIdx, c("rtmin", "rtmax"), drop=FALSE])

## ## Do it completely manual...
## rtidxRange <- range(which(xraw@scantime >= rtr[1] & xraw@scantime <= rtr[2]))
## ## @env$mz : mz values.
## ## @env$intensity: intensity values.
## ## @scantime: scantime, shorter than both above.
## ## @scanindex: which index in mz and intensity matches the time. length scanindex=length scantime.
## scanIndexRange <- xraw@scanindex[rtidxRange] + 1 ## +1 since we're counting from 0!
## ## hm, not so sure if that's correct...
## masses <- xraw@env$mz[scanIndexRange[1]:scanIndexRange[2]]
## massRange <- which(masses >= mzr[1] & masses <= mzr[2])

## ints <- xraw@env$intensity[(scanIndexRange[1]:scanIndexRange[2])[massRange]]
## ## but that's the same than we get with rawMat
## library(RUnit)
## checkEquals(rm[, "intensity"], ints)

## ## rawEIC:
## rEic <- rawEIC(xraw, rtrange=rtr, mzrange=mzr)

## ## rawMat has also the mz (which is pretty nice).
## ## how could we get that for rawEIC too?
## ## get the mz index for the scans:
## sidx <- xraw@scanindex[rEic$scan[1] + 1]
## eidx <- xraw@scanindex[max(rEic$scan) + 1]
## ## Don't get it!!!

## ## getting stuff from profile matrix is something different...
## pEic <- getEIC(xraw, rtrange=rtr, mzrange=mzr)@eic[[1]][[1]]
## pEic[pEic[, "rt"] >= rtr[1] & pEic[, "rt"] <= rtr[2], ]
## ## ???

## ## scanindex, scantime.
## ## rEic$scan corresponds to the scantime >= rtr[1] and <= rtr[2]
## sct <- xraw@scantime[rEic$scan]
## plot(sct, rEic$intensity, type="b")
## ## scantimes match!
## plot(pEic[, "rt"], pEic[, "intensity"] ,type="b")
## ## OK; we're getting closer! profile matrix does not correspond to raw data!


## ##plotRaw(xraw, mzrange=mzr, rtrange=rtr, type="b")

## ## Summary:
## ## getEIC extracts stuff from the profile matrix, rawEIC and rawMat from the raw data.
## ## rawEIC uses a C call to extract stuff.
## ## rawMat uses R code to extract stuff from the raw matrix.

## ## timing.
