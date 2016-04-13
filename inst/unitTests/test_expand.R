####============================================================
##  Testing the expand method.
##
####------------------------------------------------------------
test_grow_num <- function(){
    vals <- c(3, 4)
    res <- grow(vals, by=2)
    checkEquals(res, c(1, 6))
    checkException(grow(c(1, 2, 3)))

    ## Matrix
    mat <- rbind(c(3, 4),
                 c(1, 5),
                 c(11, 12))
    res <- grow(mat, by=2)
    checkEquals(res, rbind(c(1, 6),
                           c(-1, 7),
                           c(9, 14)))
    res <- grow(mat, by=c(1, 2, 3))
    checkEquals(res, rbind(c(2, 5),
                           c(-1, 7),
                           c(8, 15)))
}
