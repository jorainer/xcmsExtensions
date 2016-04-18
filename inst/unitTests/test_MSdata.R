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
    bins <- xcmsExtensions:::.getBins(vals, nbin=20)
    ## expect to have 100 elements in each bin
    Test <- split(vals, findInterval(vals, bins, all.inside=TRUE))
    checkEquals(unique(unlist(lapply(Test, length), use.names=FALSE)), 100)
    bins <- xcmsExtensions:::.getBins(vals, binSize=100)
    ## expect to have 100 elements in each bin
    Test <- split(vals, findInterval(vals, bins, all.inside=TRUE))
    checkEquals(unique(unlist(lapply(Test, length), use.names=FALSE)), 100)

    ##
    checkTrue(all(findInterval(1:10, bins, all.inside=TRUE) == 1 ))

    checkTrue(findInterval(c(2000), bins, all.inside=TRUE) == 20)
    checkEquals(findInterval(c(300, 200, 1000), bins, all.inside=TRUE), c(3, 2, 10))

    system.time(Test <- xcmsExtensions:::.aggregateWithBins(vals, bins=bins, FUN=mean))
    checkEquals(Test, (50*seq(1, 40, by=2))+0.5)
    checkEquals(Test+0.5, xcmsExtensions:::.binMid(bins))

    ## Check .binMeans.
    checkEquals(xcmsExtensions:::.binMid(c(2, 4, 6, 8, 10)), c(3, 5, 7, 9))
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
    bins <- xcmsExtensions:::.getBins(rtrange(msd), binSize=1)
    ints <- findInterval(rtime(msd), bins, all.inside=TRUE)
    aggs <- xcmsExtensions:::.aggregateWithIntervals(intensity(msd), ints)
    midP <- xcmsExtensions:::.binMid(bins)
    if(do.plot)
        points(midP[unique(ints)], aggs, col="blue", type="b")
    ## Plotting a spectrum: intensity (y) vs mz (x)
}


test_getChrom <- function(){
    ## Testing the internal getChrom function.
    do.plot <- FALSE

    ## Get the data.
    mzr <- c(302, 302.1)
    rtr <- c(2550, 2700)
    msd <- msData(xraw, mzrange=mzr, rtrange=rtr)

    ## There should be no need to aggregate...
    chr <- xcmsExtensions:::.getChrom(msd)
    checkEquals(chr, as.matrix(msd)[, c("rtime", "intensity")])

    ## Do binning specifying nbin
    chrBin <- xcmsExtensions:::.getChrom(msd, nbin=20)
    ## Manually specify the bins.
    bns <- xcmsExtensions:::.getBins(rtrange(msd), nbin=20)
    chrBin2 <- xcmsExtensions:::.getChrom(msd, bins=bns)
    ## That should be the same.
    checkEquals(chrBin, chrBin2)

    ## the same with FUN=mean
    chrBin <- xcmsExtensions:::.getChrom(msd, nbin=20, FUN=mean)
    ## Manually specify the bins.
    bns <- xcmsExtensions:::.getBins(rtrange(msd), nbin=20)
    chrBin2 <- xcmsExtensions:::.getChrom(msd, bins=bns, FUN=mean)
    ## That should be the same.
    checkEquals(chrBin, chrBin2)

    ## The same with binSize
    chrBin <- xcmsExtensions:::.getChrom(msd, binSize=2)
    ## Manually specify the bins.
    bns <- xcmsExtensions:::.getBins(rtrange(msd), binSize=2)
    chrBin2 <- xcmsExtensions:::.getChrom(msd, bins=bns)
    ## That should be the same.
    checkEquals(chrBin, chrBin2)

    ## Do use a larger mz range.
    mzr <- c(200, 400)
    msd <- msData(xraw, mzrange=mzr, rtrange=rtr)

    ## Get the chrom.
    chr <- xcmsExtensions:::.getChrom(msd)
    ## Manually doing that.
    manVals <- unlist(lapply(split(intensity(msd), rtime(msd)), FUN=sum), use.names=FALSE)
    checkEquals(chr[, "intensity"], manVals)

    ## And with binning.
    chrBin <- xcmsExtensions:::.getChrom(msd, nbin=20)

    if(do.plot)
        plot(chr[, 1], chr[, 2])
    if(do.plot)
        points(chrBin[, 1], chrBin[, 2], col="blue", type="b")
    chrBinMean <- xcmsExtensions:::.getChrom(msd, nbin=20, FUN=mean)
    chr <- xcmsExtensions:::.getChrom(msd, FUN=mean)
    if(do.plot)
        plot(chr[, 1], chr[, 2])
    if(do.plot)
        points(chrBinMean[, 1], chrBinMean[, 2], col="green", type="b")
}

test_getChromatogram <- function(){
    do.plot <- FALSE
    mzr <- c(200, 300)
    rtr <- c(2400, 2700)

    datmat <- msData(xraw, mzrange=mzr, rtrange=rtr)
    suppressWarnings(
        datmat2 <- msData(getXcmsRaw(xset, 2), mzrange=mzr, rtrange=rtr)
    )
    ## getChromatogram
    chr1 <- getChromatogram(datmat)
    chr2 <- getChromatogram(datmat2)

    ## Test with 40 bins.
    chr1 <- getChromatogram(datmat, nbin=40)
    chr2 <- getChromatogram(datmat2, nbin=40)
    checkEquals(nrow(chr1), 40)

    ## Test getting the same rts.
    bns <- xcmsExtensions:::.getBins(rtr, nbin=50)
    chr1 <- getChromatogram(datmat, bins=bns)
    chr2 <- getChromatogram(datmat2, bins=bns)
    checkEquals(chr1[, 1], chr2[, 1])

    ## Test plotChromatogram
    if(do.plot){
        plotChromatogram(datmat, type="b")
        plotChromatogram(datmat2, type="b", add=TRUE, col="red")
        ## Using the max signal in m/z
        plotChromatogram(datmat, type="b", FUN=max, nbin=40)
        plotChromatogram(datmat2, type="b", add=TRUE, col="red", FUN=max, nbin=40)
    }
}


test_getSpectrum <- function(){
    do.plot <- FALSE
    ## Peak would be around rt 2620-2625
    mzr <- c(300, 600)
    rtr <- c(2560, 2625)
    ## Check if we could get the spectrum of that peak.
    msd <- msData(xraw, mzrange=mzr, rtrange=rtr)
    Test <- xcmsExtensions:::.getSpec(msd)
    spec <- getSpectrum(msd)
    ## mz has to be ordered increasing
    checkEquals(order(spec[, 1]), 1:nrow(spec))
    if(do.plot){
        par(mfrow=c(1, 2))
        plotChromatogram(msd, type="l")
        plotSpectrum(msd, type="l")
        ## At least looks fine
    }
    ## Now doing a binning along the M/Z too; with binSize=1
    specB <- getSpectrum(msd, binSize=1)
    checkEquals(order(specB[, 1]), 1:nrow(specB))
    ## Check that the binSize is 1.
    checkTrue(all(diff(as.numeric(rownames(specB))) == 1))
    ## Do the plot with and without binning.
    if(do.plot){
        plotSpectrum(msd, type="l")
        plotSpectrum(msd, binSize=1, type="l", col="blue")
    }
}

test_msData2mapMatrix <- function(){
    rtr <- c(2550, 2700)
    mzr <- c(300, 330)
    ## Extrac the MS data slice
    msd <- msData(xraw, mzrange=mzr, rtrange=rtr)

    mat <- xcmsExtensions:::.msData2mapMatrix(msd)

    ## Test also with a completely re-ordered MSdata; shouldn't be a big
    ## problem.
    msd2 <- msd
    idx <- sample(1:length(rtime(msd)))
    msd2@mz <- msd2@mz[idx]
    msd2@intensity <- msd2@intensity[idx]
    msd2@rtime <- S4Vectors::Rle(rtime(msd2)[idx])
    mat2 <- xcmsExtensions:::.msData2mapMatrix(msd2)
    checkEquals(mat, mat2)

    ## What with spectrum and chromatogram
    checkEquals(getSpectrum(msd), getSpectrum(msd2))
    checkEquals(getChromatogram(msd), getChromatogram(msd2))

    ## Check if the data is as we expect it! mat is a matrix, rows are mz, columns rt.
    gotMat <- as.matrix(msd)
    ## Values from the first column have to match corresponding entries in this matrix!
    for(i in 1:ncol(mat)){
        tmp <- unname(mat[, i])
        tmp <- tmp[!is.na(tmp)]
        checkEquals(tmp, gotMat[gotMat[, 1] == as.numeric(colnames(mat)[i]), "intensity"])
    }
    tmp <- as.numeric(mat)
    checkEquals(tmp[!is.na(tmp)], gotMat[, "intensity"])

    object.size(mat)
    object.size(msd)
}

test_mapMatrix <- function(){

    ## Test extracting the data as a 2d matrix.
    rtr <- c(2550, 2700)
    mzr <- c(300, 330)
    ## Extrac the MS data slice
    msd <- msData(xraw, mzrange=mzr, rtrange=rtr)
    ## Extract the matrix: rows are MZ, cols rt
    mat <- mapMatrix(msd)
    checkEquals(as.numeric(rownames(mat)), unique(xcmsExtensions:::mzOrdered(msd)))
    checkEquals(as.numeric(colnames(mat)), unique(xcmsExtensions:::rtimeOrdered(msd)))
    ## And the values.
    vals <- as.numeric(mat)
    vals <- vals[vals!=0]
    checkEquals(length(vals), length(intensity(msd)))
    checkEquals(vals, intensity(msd))
}


notrun_compare_matrix2sparseMatrix <- function(){
    rtr <- c(2550, 2700)
    mzr <- c(300, 330)
    ## Extrac the MS data slice
    msd <- msData(xraw, mzrange=mzr, rtrange=rtr)

    ## A) same content?
    mat <- xcmsExtensions:::.msData2mapMatrix(msd)
    smat <- xcmsExtensions:::.msData2mapSparseMatrix(msd)
    idxNotNa <- which(!is.na(as.numeric(mat)))
    checkEquals(as.numeric(mat[idxNotNa]), as.numeric(smat[idxNotNa]))
    checkEquals(rownames(mat), rownames(smat))
    checkEquals(colnames(mat), colnames(smat))

    ## B) Creation time
    msd <- msData(xraw)
    system.time(mat <- xcmsExtensions:::.msData2mapMatrix(msd)) ## Needs roughly 2 seconds.
    system.time(smat <- xcmsExtensions:::.msData2mapSparseMatrix(msd)) ## Needs roughly 0.6 secs.

    ## C) memory usage.
    object.size(mat)  ## about 20MB
    object.size(smat) ## about 6MB

    ## And the winner is: sparseMatrix!!!
}

test_binRtime <- function(){

    rtr <- c(2550, 2700)
    mzr <- c(300, 330)
    ## Extrac the MS data slice
    msd <- msData(xraw, mzrange=mzr, rtrange=rtr)
    system.time(rtbinned <- xcmsExtensions:::.binRtime(msd, nbin=20))
    system.time(rtbinned2 <- xcmsExtensions:::.binRtimeSlow(msd, nbin=20))
    ## Compare with the slow version, that does return true results
    mat <- as.matrix(rtbinned)
    mat2 <- as.matrix(rtbinned2)
    mat2 <- mat2[order(mat2[, 1], mat2[, 2]), ]
    checkEquals(mat, mat2)

    ## Randomly rearranged data.
    msd2 <- msd
    idx <- sample(1:length(rtime(msd)))
    msd2@mz <- msd2@mz[idx]
    msd2@intensity <- msd2@intensity[idx]
    msd2@rtime <- S4Vectors::Rle(rtime(msd2)[idx])
    rtbinned2 <- binRtime(msd2, nbin=20)
    checkEquals(rtbinned, rtbinned2)

    ## Repeat with binSize.
    system.time(rtbinned <- binRtime(msd, binSize=1))
    system.time(rtbinned2 <- xcmsExtensions:::.binRtimeSlow(msd, binSize=1))
    ## Compare with the slow version, that does return true results
    mat <- as.matrix(rtbinned)
    mat2 <- as.matrix(rtbinned2)
    mat2 <- mat2[order(mat2[, 1], mat2[, 2]), ]
    checkEquals(mat, mat2)

    ## Randomly rearranged data.
    msd2 <- msd
    idx <- sample(1:length(rtime(msd)))
    msd2@mz <- msd2@mz[idx]
    msd2@intensity <- msd2@intensity[idx]
    msd2@rtime <- S4Vectors::Rle(rtime(msd2)[idx])
    rtbinned2 <- binRtime(msd2, binSize=1)
    checkEquals(rtbinned, rtbinned2)
}

test_binMz <- function(){

    ## Testing the internal .binMz method and evaluate whether it is
    ## binning as expected.
    rtr <- c(2550, 2700)
    mzr <- c(300, 330)
    ## Extrac the MS data slice
    msd <- msData(xraw, mzrange=mzr, rtrange=rtr)
    system.time(mzbinned <- xcmsExtensions:::.binMz(msd, nbin=20))
    system.time(mzbinned2 <- xcmsExtensions:::.binMzSlow(msd, nbin=20))
    ## Expect 20 unique mz values
    checkEquals(length(unique(mz(mzbinned))), 20)
    ## Compare with the slow version, that does return true results
    mat <- as.matrix(mzbinned)
    mat2 <- as.matrix(mzbinned2)
    mat2 <- mat2[order(mat2[, 1], mat2[, 2]), ]
    checkEquals(mat, mat2)

    ## Test what happens if we have randomly ordered data.
    msd2 <- msd
    idx <- sample(1:length(rtime(msd)))
    msd2@mz <- msd2@mz[idx]
    msd2@intensity <- msd2@intensity[idx]
    msd2@rtime <- S4Vectors::Rle(rtime(msd2)[idx])
    mzbinned2 <- xcmsExtensions:::.binMz(msd2, nbin=20)
    checkEquals(mzbinned, mzbinned2)

    ## Repeat with specifying a binSize of 1
    system.time(mzbinned <- xcmsExtensions:::.binMz(msd, binSize=1))
    ## The same but not sorting.
    system.time(mzbinned2 <- xcmsExtensions:::.binMzUnsorted(msd, binSize=1))
    ## Check with the "slow" method
    system.time(mzbinned3 <- xcmsExtensions:::.binMzSlow(msd, binSize=1))
    ## Check if we've got the same:
    mat <- as.matrix(mzbinned)
    mat2 <- as.matrix(mzbinned2)
    mat3 <- as.matrix(mzbinned3)
    ## Order the matrices... no need for mat.
    mat2 <- mat2[order(mat2[, 1], mat2[, 2]), ]
    mat3 <- mat3[order(mat3[, 1], mat3[, 2]), ]
    checkEquals(mat, mat2)
    checkEquals(mat, mat3)
    ## OK, fine.

    ## Test what happens if we have randomly ordered data.
    msd2 <- msd
    idx <- sample(1:length(rtime(msd)))
    msd2@mz <- msd2@mz[idx]
    msd2@intensity <- msd2@intensity[idx]
    msd2@rtime <- S4Vectors::Rle(rtime(msd2)[idx])
    system.time(mzbinned2 <- xcmsExtensions:::.binMz(msd2, binSize=1))
    checkEquals(mzbinned, mzbinned2)
    ## OK; fine!

    ## At last check the method itself.
    mzb <- binMz(msd, binSize=1)
    checkEquals(mzb, mzbinned)

    ## Some error checking
    checkException(binMz(msd, binSize=1, nbin=3))

    ## Check performance on the full data set.
    ## msfull <- msData(xraw)
    ## system.time(mzbinned <- xcmsExtensions:::.binMz(msfull, binSize=1))   ## Amazingly slow... took me 10secs!
    ## ## If we don't sort?
    ## system.time(mzbinnedUnsorted <- xcmsExtensions:::.binMz(msfull, binSize=1, sort=FALSE))  ## took 8secs!
    ## ## Do the slow version
    ## system.time(mzbinned3 <- xcmsExtensions:::.binMzSlow(msfull, binSize=1))  ## 30secs!
    ## mat <- as.matrix(mzbinned)
    ## mat2 <- as.matrix(mzbinnedUnsorted)
    ## mat2 <- mat2[order(mat2[, 1], mat2[, 2]), ]
    ## mat3 <- as.matrix(mzbinned3)
    ## mat3 <- mat3[order(mat3[, 1], mat3[, 2]), ]
    ## checkEquals(mat, mat2)
    ## checkEquals(mat, mat3)
}

test_binMzRtime <- function(){
    do.plot <- FALSE
    ## Bin in both dimensions.
    rtr <- c(2550, 2700)
    mzr <- c(300, 330)
    ## Extrac the MS data slice
    msd <- msData(xraw, mzrange=mzr, rtrange=rtr)
    system.time(mzbinned <- xcmsExtensions:::.binMz(msd, nbin=20))
    mzNrt <- binRtime(mzbinned, binSize=5)

    ## All in one go:
    inOne <- binMzRtime(msd, mzNbin=20, rtBinSize=5)
    checkEquals(mzNrt, inOne)
    if(do.plot){
        par(mfrow=c(1, 2))
        plotSpectrum(msd, col="grey", type="l")
        plotSpectrum(inOne, col="blue", type="l", add=TRUE, lty=2)
        plotChromatogram(msd, col="grey", type="l")
        plotChromatogram(inOne, col="blue", type="l", add=TRUE, lty=2)
    }
}

notrun_test_2mat <- function(){
    ## want to transform a MSdata into a nice 2d map.
    rtr <- c(2550, 2700)
    mzr <- c(300, 330)
    ## Extrac the MS data slice
    msd <- msData(xraw, mzrange=mzr, rtrange=rtr)

    tmp <- as.matrix(msd)
    ## Unique rt and mz
    rtU <- unique(tmp[, "rtime"])
    mzU <- sort(unique(tmp[, "mz"]))
    ## The matrix will be a rtU * mzU matrix
    nrow(tmp)
    length(rtU) * length(mzU)
    ## Apparently we don't have for all rt and/or mz a value.

    mzIDf <- data.frame(mz=mz(msd), intensity=intensity(msd))
    mzIList <- split(mzIDf, rtime(msd))

    ## Now match all of mz values of each rt-data.frame to the
    ## mzU.
    matchList <- lapply(mzIList, function(z){
        tmp <- rep(NA, length(mzU))
        tmp[match(z[, "mz"], mzU)] <- z[, "intensity"]
        return(tmp)
    })
    fullInts <- unlist(matchList, use.names=FALSE)
    ## OK, now, the first length(mzU) entries correspond to the rtU[1]
    ## and so one. Thus, rows would be mz here, columns rt!
    twoDmat <- matrix(fullInts, ncol=length(rtU), nrow=length(mzU))
    rownames(twoDmat) <- mzU
    colnames(twoDmat) <- rtU
    object.size(twoDmat)

    ## The code actually went into the .msData2mapMatrix.
}

test_subset <- function(){
    rtr <- c(2550, 2700)
    mzr <- c(300, 330)

    ## Get the full data
    msdFull <- msData(xraw)

    checkEquals(msdFull, subset(msdFull))

    ## subset by rt range
    x1 <- subset(msdFull, rtrange=rtr)
    x2 <- msData(xraw, rtrange=rtr)
    checkEquals(x1, x2)
    ## subset by mz range
    x1 <- subset(msdFull, mzrange=mzr)
    x2 <- msData(xraw, mzrange=mzr)
    checkEquals(x1, x2)
    ## subset by both
    x1 <- subset(msdFull, mzrange=mzr, rtrange=rtr)
    x2 <- msData(xraw, mzrange=mzr, rtrange=rtr)
    checkEquals(x1, x2)
}

