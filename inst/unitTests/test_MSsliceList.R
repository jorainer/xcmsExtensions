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
    checkEquals(unlist(rtranges(dummy)), rtrange(mss))
    checkEquals(unlist(mzranges(dummy)), mzrange(mss))
    checkEquals(unlist(intranges(dummy)), intrange(mss))

    dummy <- MSsliceList(mss, mss)
    checkEquals(length(dummy), 2)
    checkEquals(unlist(mzranges(dummy)), rep(mzrange(mss), 2))

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
    rtrs <- rtranges(msl)
    for(i in 1:length(rtrs)){
        checkTrue(rtrs[[i]][1] >= rtm[i, 1])
        checkTrue(rtrs[[i]][2] <= rtm[i, 2])
    }

    msl <- msSlice(xraw, mzrange=mzm)
    mzrs <- mzranges(msl)
    for(i in 1:length(mzrs)){
        checkTrue(mzrs[[i]][1] >= mzm[i, 1])
        checkTrue(mzrs[[i]][2] <= mzm[i, 2])
    }

    msl <- msSlice(xraw, mzrange=mzm, rtrange=rtm)
    mzrs <- mzranges(msl)
    for(i in 1:length(mzrs)){
        checkTrue(mzrs[[i]][1] >= mzm[i, 1])
        checkTrue(mzrs[[i]][2] <= mzm[i, 2])
    }
    rtrs <- rtranges(msl)
    for(i in 1:length(rtrs)){
        checkTrue(rtrs[[i]][1] >= rtm[i, 1])
        checkTrue(rtrs[[i]][2] <= rtm[i, 2])
    }
}

