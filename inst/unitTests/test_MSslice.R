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

