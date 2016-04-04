####============================================================
##  Testing MSslice
##
####------------------------------------------------------------
detach("package:xcmsExtensions", unload=TRUE)
library(xcmsExtensions)
library(faahKO)
library(RUnit)
xset <- faahko
suppressWarnings(
    xraw <- getXcmsRaw(xset, 1)
)

test_MSslice <- function(){
    ## Generate one from a simple subset.
    rtrange <- c(2600, 2670)
    mzrange <- c(300, 400)
    datmat <- getData(xraw, rtrange=rtrange, mzrange=mzrange)
    msd <- MSdata(mz=datmat[, "mz"], intensity=datmat[, "intensity"],
                  rtime=datmat[, "time"])

    ## Convert that into an MSslice.
    mss <- msSlice(msd)

    mss2 <- msSlice(list(msd))

    checkEquals(mss, mss2)
    checkEquals(rtrange(mss2), rtrange(msd))
    checkEquals(mzrange(mss2), mzrange(msd))
}

## test_plotChromatogram <- function(){
##     mzr <- c(302, 302.5)
##     rtr <- c(2500, 2700)
##     mss <- msSlice(xset, mzrange=mzr, rtrange=rtr)
## }

test_getChromatogram <- function(){
    do.plot <- FALSE
    mzr <- c(302, 302.5)
    rtr <- c(2500, 2700)
    suppressWarnings(
        mss <- msSlice(xset, mzrange=mzr, rtrange=rtr)
    )
    Test <- lapply(msData(mss), getChromatogram)

    Test <- xcmsExtensions:::.getChromList(mss)

    Test <- xcmsExtensions:::.getChromList(mss, nbin=20)
    checkEquals(unique(unlist(lapply(Test, nrow))), 20)
    Test2 <- xcmsExtensions:::.getChromList(mss, nbin=20, FUN=sum)
    checkTrue(all(Test[[1]][, 2] < Test2[[1]][, 2]))

    ## In that case also the rtime has to be identical!
    rtimes <- unique(unlist(lapply(Test, function(z)z[, 1])))
    checkEquals(length(rtimes), 20)

    ## getChromatogram for an MSslice is supposed to return a matrix with ncol
    ## being the number of samples, nrow the number of distinct scan times (across
    ## all samples) or number of bins (preferred way).
    ## Getting the full chromatogram.
    chrM <- getChromatogram(mss)
    checkEquals(ncol(chrM), length(mss))
    ## Do the binning in rt-space; bin it into 100 bins.
    chrM <- getChromatogram(mss, nbin=100)
    checkEquals(nrow(chrM), 100)
    ## Just visually checking...
    if(do.plot){
        plotChromatogram(msData(mss)[[1]], type="l")
        points(as.numeric(rownames(chrM)), chrM[, 1], col="blue", type="l")
        ## Wonderful!!!
        for(i in 1:ncol(chrM)){
            points(as.numeric(rownames(chrM)), chrM[, i], col="grey", type="l")
        }
        plotChromatogram(mss)
    }
    ## Check if we want to get rt spacing of 1 (sec?)
    chrM <- getChromatogram(mss, binSize=1)
    ## Minimal spacing between rts is 1. Can not assume that all spacings are 1, since we
    ## might have rts without measurements.
    checkTrue(min(diff(as.numeric(rownames(chrM)))) == 1)
}

test_getSpectrum <- function(){
    do.plot <- FALSE
    mzr <- c(300, 303)
    rtr <- c(2500, 2700)
    suppressWarnings(
        mss <- msSlice(xset, mzrange=mzr, rtrange=rtr)
    )
    spcM <- getSpectrum(mss)

    suppressWarnings(
        Test <- getSpectrum(msData(getXcmsRaw(xset, 3), mzrange=mzr, rtrange=rtr))
    )
    rownames(Test) <- Test[, "mz"]
    checkEquals(spcM[!is.na(spcM[, 3]), 3], Test[, "intensity"])
    if(do.plot){
        plotSpectrum(mss, type="h", col="#00000020")
    }
    ## Has to be increasingly.
    checkEquals(order(as.numeric(rownames(spcM))), 1:nrow(spcM))
}
