detach("package:xcmsExtensions", unload=TRUE)
library(xcmsExtensions)
library(faahKO)
library(RUnit)
xset <- faahko
suppressWarnings(
    xraw <- getXcmsRaw(xset, 1)
)

test_scantimes <- function(){
    sct <- scantimes(xraw)
    checkEquals(unique(sct), xraw@scantime)
    ## Check if I've got the correct length:
    checkEquals(length(sct), length(xraw@env$mz))
    ## Correct number of times...
    lengths <- unlist(lapply(split(sct, sct), length), use.names=FALSE)
    checkEquals(diff(xraw@scanindex), lengths[-length(lengths)])
}

test_getData <- function(){
    ## Get the full matrix
    system.time(fullmat <- getData(xraw))
    ## Get the full matrix with rawMat
    system.time(rawmat <- rawMat(xraw))
    ## So, we're 5times faster.
    checkEquals(fullmat, rawmat)

    rtr <- c(2509, 2530)
    system.time(datmat <- getData(xraw, rtrange=rtr))
    system.time(rawmat <- rawMat(xraw, rtrange=rtr))
    ## Dash it; so they are 6times faster.
    ## What the heck! Results are NOT the same!!!
    ## That's essentially a bug in rawMat!!!
    checkEquals(datmat,
                fullmat[fullmat[, "time"] >= rtr[1] & fullmat[, "time"] <= rtr[2], ])


    mzr <- c(301, 302.003)
    system.time(datmat <- getData(xraw, mzrange=mzr))
    system.time(rawmat <- rawMat(xraw, mzrange=mzr))
    ## OK, we're faster again.
    checkEquals(datmat, rawmat)
    checkEquals(datmat,
                fullmat[fullmat[, "mz"] >= mzr[1] & fullmat[, "mz"] <= mzr[2], ])

    ## Combine both
    system.time(datmat <- getData(xraw, mzrange=mzr, rtrange=rtr))
    system.time(rawmat <- rawMat(xraw, mzrange=mzr, rtrange=rtr))
    checkEquals(datmat,
                fullmat[fullmat[, "mz"] >= mzr[1] &
                        fullmat[, "mz"] <= mzr[2] &
                        fullmat[, "time"] >= rtr[1] &
                        fullmat[, "time"] <= rtr[2], ])

    ## Intensity range:
    system.time(datmat <- getData(xraw, intrange=c(300, Inf)))
    checkTrue(all(datmat[, "intensity"] >= 300))
}

## Test what happens if we extract data and return it as MSmap objects.
test_msMap <- function(){
    do.plot <- FALSE
    ## use the peak from the vignette (Figure 4)
    rtr <- c(2550, 2700)
    mzr <- c(301, 304.1)
    mzRes <- 0.01
    datm <- getData(xraw, rtrange=rtr, mzrange=mzr)
    if(do.plot)
        plot(datm[, "time"], datm[, "intensity"], type="b")
    ## obviously not equally spaced.
    ## Get the EIC for that:
    profStep(xraw) <- mzRes
    eicm <- getEIC(xraw, rtrange=matrix(rtr, nrow=1), mzrange=matrix(mzr, nrow=1))
    if(do.plot)
        plot(eicm)
    ## Looks similar, but have far more data points!
    ## Get the MSmap
    msm <- msMap(xraw, rtrange=rtr, mzrange=mzr, resMz=mzRes)
    checkEquals(rtime(msm), unique(datm[, "time"]))

    reic <- rawEIC(xraw, rtrange=matrix(rtr, nrow=1), mzrange=matrix(mzr, nrow=1))
    ## Compare that with the EIC:
    eic <- eicm@eic[[1]][[1]]
    if(do.plot)
        plot(msm, allTicks=FALSE)
    ## if(do.plot)
    ##     plot3D(msm, rgl=TRUE)
}


notrun_test_data2MSmap <- function(){
    ## Get the full matrix.
    datmat <- getData(xraw)
    ## transform to a 2dimensional data.
    ## Means however that I have
    profStep <- 0.01
    profStep(xraw) <- profStep
    mzr <- range(xraw@env$mz)
    minM <- round(mzr[1]/profStep) * profStep
    maxM <- round(mzr[2]/profStep) * profStep
    numC <- (maxM-minM)/profStep + 1
    TestProf <- xcms:::profBinM(xraw@env$mz, xraw@env$intensity, xraw@scanindex,
                                numC, minM, maxM, FALSE, list())
    checkEquals(TestProf, xraw@env$profile)
    ## Roll my own:
    system.time(
        rl <- runLength(Rle(datmat[, "time"]))
    )
    system.time(
        myidx <- c(0, cumsum(rl)[-length(rl)])
    )
    checkEquals(myidx, xraw@scanindex)
    myprof <- xcms:::profBinM(datmat[, "mz"], datmat[, "intensity"], myidx, numC, minM, maxM, FALSE, list())
    checkEquals(myprof, TestProf)
}


## Test with arguments being matrices
test_getData_matrix <- function(){
    rtm <- rbind(c(2509, 2530),
                 c(2534, 2540),
                 c(2630, 2645))
    mzm <- rbind(c(301, 302.003),
                 c(304, 305),
                 c(308, 308.3))
    intm <- rbind(c(100, Inf),
                  c(300, 900),
                  c(-Inf, Inf))
    ## Check exceptions
    checkException(getData(xraw, mzrange=mzm, rtrange=c(2, 5)))
    checkException(getData(xraw, mzrange=mzm, rtrange=rtm[1:2, ]))
    ## Just rtm
    datmat <- getData(xraw, rtrange=rtm)
    checkEquals(datmat[[1]], getData(xraw, rtrange=rtm[1, ]))
    checkEquals(datmat[[2]], getData(xraw, rtrange=rtm[2, ]))
    checkEquals(datmat[[3]], getData(xraw, rtrange=rtm[3, ]))
    ## Just mzm
    datmat <- getData(xraw, mzrange=mzm)
    checkEquals(datmat[[1]], getData(xraw, mzrange=mzm[1, ]))
    checkEquals(datmat[[2]], getData(xraw, mzrange=mzm[2, ]))
    checkEquals(datmat[[3]], getData(xraw, mzrange=mzm[3, ]))
    ## Just intm
    datmat <- getData(xraw, intrange=intm)
    checkEquals(datmat[[1]], getData(xraw, intrange=intm[1, ]))
    checkEquals(datmat[[2]], getData(xraw, intrange=intm[2, ]))
    checkEquals(datmat[[3]], getData(xraw, intrange=intm[3, ]))
    ## Combinations.
    datmat <- getData(xraw, rtrange=rtm, mzrange=mzm)
    checkEquals(datmat[[1]], getData(xraw, rtrange=rtm[1, ], mzrange=mzm[1, ]))
    checkEquals(datmat[[2]], getData(xraw, rtrange=rtm[2, ], mzrange=mzm[2, ]))
    checkEquals(datmat[[3]], getData(xraw, rtrange=rtm[3, ], mzrange=mzm[3, ]))
}


notrun_xcms_bug_rawMat <- function(){
    ## Get the full matrix
    rawmat <- rawMat(xraw)
    ## Define a rtrange
    rtr <- c(2509, 2530)
    submat <- rawMat(xraw, rtrange=rtr)
    ## Subsetting the full matrix:
    submat2 <- rawmat[rawmat[, "time"] >= rtr[1] & rawmat[, "time"] <= rtr[2], ]
    ## These have to be identical!
    checkEquals(submat, submat2)
}


