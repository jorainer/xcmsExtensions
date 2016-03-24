####============================================================
##  Testing utility methods and functions
##
##
####------------------------------------------------------------
#### Test cases...
library(faahKO)
xset <- faahko
xset <- group(xset)
xset <- retcor(xset)
xset <- group(xset)

## Testing if we get similar/same results as matchpeaks from xcms.
test_matchmz <- function(){
    ## Get
}

