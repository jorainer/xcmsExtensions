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
globalMzr <- c(300, 320)
globalRtr <- c(2500, 2700)
suppressWarnings(
    globalMss <- msSlice(xset, mzrange=globalMzr, rtrange=globalRtr)
)
suppressWarnings(
    mssFull <- msSlice(xset)
)

test_names <- function(){
    mss <- globalMss
    checkEquals(names(mss), NULL)

    checkException(names(mss) <- c("a", "b"))
    names(mss) <- rev(1:12)
    checkEquals(names(mss), as.character(rev(1:12)))
}

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

test_bracketSubset <- function(){
    mss <- globalMss
    checkException(mss[3, 4])
    msss <- mss[c(1, 2, 5)]
    checkEquals(length(msss@data), 3)
    checkEquals(msss@data[[1]], mss@data[[1]])
    checkEquals(msss@data[[3]], mss@data[[5]])

    ## Repeating several idxs...
    msss <- mss[c(1, 1, 1)]
    checkEquals(msss@data[[1]], mss@data[[1]])
    checkEquals(msss@data[[2]], mss@data[[1]])

    ## Subsetting by name.
    checkException(mss[c("first", "second")])
    names(mss) <- sampnames(xset)
    msss <- mss[c("ko18", "ko22", "wt16")]
    checkEquals(length(msss), 3)
    msss2 <- mss[c(3, 6, 8)]
    checkEquals(msss, msss2)

    ## Use logical vector to subset.
    lv <- rep(FALSE, length(mss))
    lv[c(3, 8, 9)] <- TRUE
    msss <- mss[lv]
    checkEquals(msss, mss[c(3, 8, 9)])

    ## [[
    checkException(mss[[c(1, 2)]])
    checkException(mss[[1, 2]])

    ## OK, now something that works...
    msss <- mss[["ko16"]]
    checkEquals(msss, mss@data[[2]])
    msss <- mss[[12]]
    checkEquals(msss, mss@data[[12]])
}

## test_plotChromatogram <- function(){
##     mzr <- c(302, 302.5)
##     rtr <- c(2500, 2700)
##     mss <- msSlice(xset, mzrange=mzr, rtrange=rtr)
## }

test_getChromatogram <- function(){
    do.plot <- FALSE
    mzr <- globalMzr
    rtr <- globalRtr
    mss <- globalMss
    ## mzr <- c(302, 302.5)
    ## rtr <- c(2500, 2700)
    ## suppressWarnings(
    ##     mss <- msSlice(xset, mzrange=mzr, rtrange=rtr)
    ## )
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
    mzr <- globalMzr
    rtr <- globalRtr
    mss <- globalMss
    ## mzr <- c(300, 303)
    ## rtr <- c(2500, 2700)
    ## suppressWarnings(
    ##     mss <- msSlice(xset, mzrange=mzr, rtrange=rtr)
    ## )
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

test_binMz <- function(){
    do.plot <- FALSE
    mzr <- globalMzr
    rtr <- globalRtr
    mss <- globalMss

    ## Test the binning and also evaluate the @call
    bins <- xcmsExtensions:::.getBins(mzrange(mss), nbin=10)
    msd <- msData(xraw, rtrange=rtr, mzrange=mzr)
    msdBinned <- binMz(msd, bins=bins)

    mssBinned <- binMz(mss, nbin=10)
    checkEquals(mzrange(msdBinned), mzrange(mssBinned))

    ## Check if the binned data is the same...
    checkEquals(as.matrix(msdBinned), as.matrix(msData(mssBinned)[[1]]))
}

test_binRtime <- function(){
    do.plot <- FALSE
    mzr <- globalMzr
    rtr <- globalRtr
    mss <- globalMss

    ## Test the binning and also evaluate the @call
    bins <- xcmsExtensions:::.getBins(rtrange(mss), binSize=5)
    msd <- msData(xraw, rtrange=rtr, mzrange=mzr)
    msdBinned <- binRtime(msd, bins=bins)

    mssBinned <- binRtime(mss, binSize=5)
    checkEquals(mzrange(msdBinned), mzrange(mssBinned))
    checkEquals(rtrange(msdBinned), rtrange(mssBinned))

    ## Check if the binned data is the same...
    checkEquals(as.matrix(msdBinned), as.matrix(msData(mssBinned)[[1]]))
}


test_binMzRtime <- function(){
    do.plot <- FALSE
    mzr <- globalMzr
    rtr <- globalRtr
    mss <- globalMss

    mzB <- binMz(mss, nbin=10)
    mzRtB <- binRtime(mzB, binSize=5)

    mssB <- binMzRtime(mss, mzNbin=10, rtBinSize=5)
    checkEquals(as.matrix(msData(mzRtB)[[1]]), as.matrix(msData(mssB)[[1]]))
}

test_subset <- function(){

    ## Get the full Data.
    checkEquals(globalMss, subset(globalMss))

    ## Subset by rt range
    x1 <- subset(mssFull, rtrange=globalRtr)
    x2 <- msSlice(xset, rtrange=globalRtr)
    checkEquals(x1, x2)

    ## Subset by mz range
    x1 <- subset(mssFull, mzrange=globalMzr)
    x2 <- msSlice(xset, mzrange=globalMzr)
    checkEquals(x1, x2)

    ## Subset by both
    x1 <- subset(mssFull, mzrange=globalMzr, rtrange=globalRtr)
    x2 <- msSlice(xset, mzrange=globalMzr, rtrange=globalRtr)
    checkEquals(x1, x2)
}
