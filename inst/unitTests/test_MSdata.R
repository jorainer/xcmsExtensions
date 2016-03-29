####============================================================
##  Testing MSdata
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

test_MSdata <- function(){
    ## Get the full data.
    datmat <- getData(xraw)
    ## Generate the MSdata for the full data set.
    msd <- MSdata(mz=datmat[, "mz"], intensity=datmat[, "intensity"],
                  rtime=datmat[, "time"])
    object.size(xraw)
    object.size(datmat)
    object.size(msd)

    ## Generate one from a simple subset.
    rtrange <- c(2600, 2670)
    mzrange <- c(300, 400)
    datmat <- getData(xraw, rtrange=rtrange, mzrange=mzrange)
    msd <- MSdata(mz=datmat[, "mz"], intensity=datmat[, "intensity"],
                  rtime=datmat[, "time"])
    checkEquals(mz(msd), datmat[, "mz"])
    checkEquals(intensity(msd), datmat[, "intensity"])
    checkEquals(rtime(msd), datmat[, "time"])

    ## Exceptions:
    checkException(mz(msd) <- c(1, 2, 4))
    checkException(intensity(msd) <- c(1, 2, 4))
    checkException(rtime(msd) <- c(1, 2, 4))

    rtime(msd) <- rtime(msd) + 10

    checkEquals(rtrange(msd), range(rtime(msd)))
    checkEquals(mzrange(msd), range(mz(msd)))

    ## Convert that into an MSslice.
    mss <- msSlice(msd)
}

## Directly extract MSdata from a xcmsRaw object.
test_msData <- function(){
    ## Get the full stuff.
    full <- msData(xraw)

    ## Subset.
    rtrange <- c(2600, 2670)
    mzrange <- c(300, 400)
    datmat <- getData(xraw, rtrange=rtrange, mzrange=mzrange)
    mssub <- msData(xraw, rtrange=rtrange, mzrange=mzrange)
    mssubm <- as.matrix(mssub)
    colnames(mssubm) <- c("time", "mz", "intensity")
    checkEquals(datmat, mssubm)

    ## Error checking.
    checkException(msData(xraw, rtrange=c(1, 2, 3)))
    checkException(msData(xraw, rtrange=c(1, 2, 3), mzrange=2))
    checkException(msData(xraw, mzrange=2))
    checkException(msData(xraw, rtrange=rbind(c(1, 2), c(2,3)), mzrange=2))

    ## That's fine.
    ## Getting multiple stuff.
    rtm <- rbind(rtrange, rtrange+100)
    mzm <- rbind(mzrange, mzrange + 100)
    multis <- msData(xraw, rtrange=rtm)
    checkEquals(length(multis), 2)
    multis <- msData(xraw, mzrange=mzm)
    multis <- msData(xraw, mzrange=mzm, rtrange=rtm)

    ## What if we're outside???
    Test <- msData(xraw, mzrange=c(0, 10))
    Test <- msData(xraw, rtrange=c(0, 1))
}

test_binning <- function(){
    vals <- 1:2000
    bins <- .getBins(vals, nbin=20)
    ## expect to have 100 elements in each bin
    Test <- split(vals, findInterval(vals, bins, all.inside=TRUE))
    checkEquals(unique(unlist(lapply(Test, length), use.names=FALSE)), 100)
    bins <- .getBins(vals, binSize=100)
    ## expect to have 100 elements in each bin
    Test <- split(vals, findInterval(vals, bins, all.inside=TRUE))
    checkEquals(unique(unlist(lapply(Test, length), use.names=FALSE)), 100)

    ##
    checkTrue(all(findInterval(1:10, bins, all.inside=TRUE) == 1 ))

    checkTrue(findInterval(c(2000), bins, all.inside=TRUE) == 20)
    checkEquals(findInterval(c(300, 200, 1000), bins, all.inside=TRUE), c(3, 2, 10))

    system.time(Test <- .aggregateWithBins(vals, bins=bins, FUN=mean))
    checkEquals(Test, (50*seq(1, 40, by=2))+0.5)
    checkEquals(Test, .binMeans(bins))

    ## Check .binMeans.
    checkEquals(.binMeans(c(2, 4, 6, 8, 10)), c(3, 5, 7, 9))
}


test_plot_chrom <- function(){
    do.plot <- FALSE
    ## Plotting a chromatogram: intensity (y) vs rt (x)
    ## Use the peak from the example data.
    mzr <- c(302, 302.1)
    rtr <- c(2550, 2700)
    msd <- msData(xraw, mzrange=mzr, rtrange=rtr)

    if(do.plot)
        plot(rtime(msd), intensity(msd), type="b")

    ## Now doing a binning of the data.
    bins <- .getBins(rtrange(msd), binSize=1)
    ints <- findInterval(rtime(msd), bins, all.inside=TRUE)
    aggs <- .aggregateWithIntervals(intensity(msd), ints)
    midP <- .binMid(bins)

    points(midP[unique(ints)], aggs, col="blue", type="b")
    ## Plotting a spectrum: intensity (y) vs mz (x)
}
