####============================================================
##  Test MSsliceList related things.
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

test_MSsliceList <- function(){
    ## Basic constructor.
    dummy <- MSsliceList()
    checkEquals(length(dummy), 0)

    checkException(MSsliceList(2))
    ## Get an MSslice:
    mss <- msSlice(xraw, mzrange=c(200, 300))
    dummy <- MSsliceList(mss)
    checkEquals(length(dummy), 1)

    dummy <- MSsliceList(mss, mss)
    checkEquals(length(dummy), 2)

    slices(dummy) <- mss
    checkEquals(length(dummy), 1)

    ## I should get a slice list if specifying a matrix of ranges...
    rtm <- rbind(c(2509, 2530),
                 c(2534, 2540),
                 c(2630, 2645))
    mzm <- rbind(c(301, 302.003),
                 c(304, 305),
                 c(308, 308.3))

    msl <- msSlice(xraw, rtrange=rtm)
    ## The rtranges have to be within the specified range.

    msl <- msSlice(xraw, mzrange=mzm)
}

