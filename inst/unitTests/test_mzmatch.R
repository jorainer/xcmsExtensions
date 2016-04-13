####============================================================
##  Test matching of mz values...
##
####------------------------------------------------------------
detach("package:xcmsExtensions", unload=TRUE)
library(xcmsExtensions)
## library(faahKO)
## library(RUnit)
## xset <- faahko
## suppressWarnings(
##     xraw <- getXcmsRaw(xset, 1)
## )


test_matching <- function(){

    ## Testing the lowest level function.
    gotMz <- 2.23
    matchWith <- seq(1, 10, by=0.1)
    ## Will not find anything with the "standard" settings.
    res <- xcmsExtensions:::.singleMatchMz(gotMz, matchWith)
    checkTrue(is.na(res[, 1]))
    ## Tune the ppm so that we find 2.2.
    accPpm <- 0.03 * 1000000 / gotMz
    res <- xcmsExtensions:::.singleMatchMz(gotMz, matchWith, mzdev=0, ppm=accPpm)
    checkEquals(nrow(res), 1)
    checkEquals(unname(res[1, 1]), which(matchWith == 2.2))
    ## Increasing the search range to include also 2.3; setting mzdev to 0.4
    res <- xcmsExtensions:::.singleMatchMz(gotMz, matchWith, mzdev=0.04, ppm=accPpm)
    checkEquals(nrow(res), 2)
    checkEquals(unname(res[, 1]), which(matchWith %in% c(2.2, 2.3)))

    ## Just using the mzdev and skip the ppm
    res <- xcmsExtensions:::.singleMatchMz(gotMz, matchWith, mzdev=0.04)
    checkEquals(nrow(res), 1)
    checkEquals(unname(res[1, 1]), which(matchWith == 2.2))

    res <- xcmsExtensions:::.singleMatchMz(gotMz, matchWith, mzdev=0.07)
    checkEquals(nrow(res), 2)
    checkEquals(unname(res[, 1]), which(matchWith %in% c(2.2, 2.3)))

    ## Just trying a second value:
    gotMz <- 2.53
    res <- xcmsExtensions:::.singleMatchMz(gotMz, matchWith)
    checkTrue(is.na(res[, 1]))
    res <- xcmsExtensions:::.singleMatchMz(gotMz, matchWith, ppm=0, mzdev=0.03)
    checkEquals(nrow(res), 1)
    checkEquals(unname(res[1, 1]), 16)
    ## with mzdev=0.07 I should get 2.5 and 2.6
    res <- xcmsExtensions:::.singleMatchMz(gotMz, matchWith, ppm=0, mzdev=0.07)
    checkEquals(nrow(res), 2)

    ## OK, now testing for x being a RANGE: in this case we're looking for a mass
    ## that's in the middle of these, i.e. at 2.24, but using an mzdev of 0.01.
    gotMz <- c(2.23, 2.25)
    res <- xcmsExtensions:::.singleMatchMz(gotMz, matchWith, mzdev=0)
    checkTrue(is.na(res[, 1]))
    ## That way we find only 2.2
    res <- xcmsExtensions:::.singleMatchMz(gotMz, matchWith, mzdev=0, ppm=accPpm)
    checkEquals(nrow(res), 1)
    checkEquals(unname(res[1, 1]), which(matchWith == 2.2))
    ## To get also 2.3 we need to increase the mzdev a little, i.e. from 2.25 to 2.3
    res <- xcmsExtensions:::.singleMatchMz(gotMz, matchWith, mzdev=0.05, ppm=0)
    checkEquals(nrow(res), 2)
    checkEquals(unname(res[, 1]), which(matchWith %in% c(2.2, 2.3)))

    ## While I wouldn't find the 2.3 if I used just the midpoint.
    res <- xcmsExtensions:::.singleMatchMz(2.24, matchWith, mzdev=0.05, ppm=0)
    checkEquals(nrow(res), 1)


    ## Test the multi matching stuff...
    gotMz <- c(2.23, 2.53)
    res <- xcmsExtensions:::.multiMatchMzNumeric(gotMz, matchWith)
    checkEquals(length(res), 2)
    res <- do.call(rbind, res)
    checkTrue(all(is.na(res[, 1])))

    ## Now doing serious stuff.
    res <- xcmsExtensions:::.multiMatchMzNumeric(gotMz, matchWith, ppm=accPpm)
    res <- do.call(rbind, res)
    checkEquals(unname(res[, 1]), which(matchWith %in% c(2.2, 2.5)))

    ## Just check what happens if we revert...
    res <- xcmsExtensions:::.multiMatchMzNumeric(gotMz[c(2, 1)], matchWith, ppm=accPpm)
    res <- do.call(rbind, res)
    checkEquals(unname(res[, 1]), c(16, 13))

    ## Get multiple stuff...
    res <- xcmsExtensions:::.multiMatchMzNumeric(gotMz, matchWith, mzdev=0.07, ppm=0)
    res <- do.call(rbind, res)
    checkEquals(unname(res[, 1]), c(13, 14, 16, 17))
}

test_matchmz_numeric <- function(){
    ## Repeating some of the stuff from above.

    gotMz <- 2.23
    matchWith <- seq(1, 10, by=0.1)
    ## Will not find anything with the "standard" settings.
    res <- mzmatch(gotMz, matchWith)
    res <- do.call(rbind, res)
    checkTrue(is.na(res[, 1]))
    ## Tune the ppm so that we find 2.2.
    accPpm <- 0.03 * 1000000 / gotMz
    res <- mzmatch(gotMz, matchWith, mzdev=0, ppm=accPpm)
    res <- do.call(rbind, res)
    checkEquals(nrow(res), 1)
    checkEquals(unname(res[1, 1]), which(matchWith == 2.2))
    ## Increasing the search range to include also 2.3; setting mzdev to 0.4
    res <- mzmatch(gotMz, matchWith, mzdev=0.04, ppm=accPpm)
    res <- res[[1]]
    checkEquals(nrow(res), 2)
    checkEquals(unname(res[, 1]), which(matchWith %in% c(2.2, 2.3)))

    gotMz <- c(2.23, 2.53)
    res <- mzmatch(gotMz, matchWith)
    checkEquals(length(res), 2)
    res <- do.call(rbind, res)
    checkTrue(all(is.na(res[, 1])))

    ## Now doing serious stuff.
    res <- mzmatch(gotMz, matchWith, ppm=accPpm)
    res <- do.call(rbind, res)
    checkEquals(unname(res[, 1]), which(matchWith %in% c(2.2, 2.5)))

    ## Just check what happens if we revert...
    res <- mzmatch(gotMz[c(2, 1)], matchWith, ppm=accPpm)
    res <- do.call(rbind, res)
    checkEquals(unname(res[, 1]), c(16, 13))

    ## Get multiple stuff...
    res1 <- mzmatch(gotMz, matchWith, mzdev=0.07, ppm=0)
    res2 <- mzmatch(gotMz, matchWith, mzdev=0.07)
    checkEquals(res1, res2)

    res <- do.call(rbind, res1)
    checkEquals(unname(res[, 1]), c(13, 14, 16, 17))
}

test_matchmz_matrix <- function(){

    ## Jump directly into the multi stuff....
    matchWith <- seq(1, 10, by=0.1)
    xM <- rbind(c(2.23, 2.25),
                c(3.23, 3.25),
                c(2.53, 2.55))
    res <- mzmatch(xM, matchWith)
    res <- do.call(rbind, res)
    checkTrue(all(is.na(res[, 1])))
    ## Expanding by 0.03 -> should get 1 each
    res <- mzmatch(xM, matchWith, mzdev=0.03)
    checkTrue(unique(unlist(lapply( res, nrow))) == 1)
    res <- do.call(rbind, res)
    checkEquals(unname(res[, 1]), c(13, 23, 16))
    ## Expanding by 0.05 -> should get 2 each
    res <- mzmatch(xM, matchWith, mzdev=0.05, ppm=0)
    checkTrue(unique(unlist(lapply( res, nrow))) == 2)
    res <- do.call(rbind, res)
    checkEquals(unname(res[, 1]), c(13, 14, 23, 24, 16, 17))

}



## testApply <- function(x){
##     tmp <- apply(x, MARGIN=1, FUN=function(z){
##         return(xcmsExtensions:::.singleMatchMz(z, matchWith, mzdev=0, ppm=accPpm))
##     })
##     return(tmp)
## }

## testLapply <- function(x){
##     tmp <- lapply(split(x, f=1:nrow(x)), function(z){
##         return(xcmsExtensions:::.singleMatchMz(z, matchWith, mzdev=0, ppm=accPpm))
##     })
##     return(tmp)
## }

## Test <- matrix(rep(c(2.23, 2.25), each=10000), ncol=2)
## system.time(testApply(Test))
## system.time(testLapply(Test))
## ## lapply slightly faster!
