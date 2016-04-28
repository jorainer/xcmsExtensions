detach("package:xcmsExtensions", unload=TRUE)
library(xcmsExtensions)
library(faahKO)
library(RUnit)
xset <- faahko
suppressWarnings(
    xraw <- getXcmsRaw(xset, 1)

)

test_peakGroupSummary <- function(){
    checkException(grpsum <- peakGroupSummary(xset))
    ## group that
    xs <- group(xset)
    grpsum <- peakGroupSummary(xs)
    checkEquals(nrow(grpsum), nrow(groups(xs)))
    ## Check if the numbers match.
    checkEquals(unname(grpsum[, "KO.propPresent"]*6), unname(groups(xs)[, "KO"]))
    checkEquals(unname(grpsum[, "WT.propPresent"]*6), unname(groups(xs)[, "WT"]))
}

test_msSlice <- function(){
    all <- msSlice(xset, rt="raw")

    ## Define some ranges...
    rtm <- rbind(c(2509, 2530),
                 c(2534, 2540),
                 c(2630, 2645),
                 c(2740, 2760),
                 c(3200, 3300),
                 c(14000, 16000))
    mzm <- rbind(c(301, 302.003),
                 c(304, 305),
                 c(308, 308.3),
                 c(308, 309),
                 c(350, 354),
                 c(309, 309.9))

    ## Check some exceptions
    checkException(msSlice(xset, rtrange=rtm[1, ], mzrange=mzm))
    ## What if we specify a wrong range?
    mss <- msSlice(xset, rtrange=c(1, 2))
    checkEquals(length(msData(mss)), length(filepaths(xset)))

    ## Just get a single MSslice object from a subset.
    system.time(mss <- msSlice(xset, rtrange=rtm[1, ], rt="raw"))
    bpp <- bpparam()
    bpworkers(bpp) <- 1
    system.time(mss <- msSlice(xset, rtrange=rtm[1, ], rt="raw", BPPARAM=bpp))
    system.time(mss <- msSlice(xset, rtrange=rtm[1, ], mzrange=mzm[1, ],
                               rt="raw", BPPARAM=bpp))

    ## Getting a MSsliceList.
    msl <- msSlice(xset, rtrange=rtm)
    msl <- msSlice(xset, rtrange=rtm, mzrange=mzm)

    ## Test sub-setting
    checkException(msl[1:10])
    mslsub <- msl[c(1, 4, 5)]

    checkEquals(slices(msl)[1], slices(mslsub)[1])
    checkEquals(slices(msl)[4], slices(mslsub)[2])
    checkEquals(slices(msl)[5], slices(mslsub)[3])

    mss <- msl[[1]]
    checkTrue(is(mss, "MSslice"))

    mss <- msl[[3]]
    checkEquals(mss, slices(msl)[[3]])

    ## Get a single one; that should be a MSslice
    mss <- msSlice(xset, rtrange=rtm[1, ], mzrange=mzm[1, ])

    ## Plot that.
    ## GO ON CHECKING HERE.
}

