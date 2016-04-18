notrun_test_shiny <- function(){

    library(xcmsExtensions)
    library(xcms)
    library(faahKO)
    cdfpath <- system.file("cdf", package="faahKO")
    cdffiles <- dir(cdfpath, recursive=TRUE, full.names=TRUE)
    ## Extract the genotype from the file name
    genot <- rep("KO", length(cdffiles))
    genot[grep(cdffiles, pattern="WT")] <- "WT"
    ## And the sample name.
    tmp <- strsplit(cdffiles, split=.Platform$file.sep)
    sampn <- unlist(lapply(tmp, function(z){
        return(gsub(z[length(z)], pattern=".CDF", replacement="", fixed=TRUE))
    }))
    ## And the phenodata table
    pheno <- data.frame(file=cdffiles, genotype=genot, name=sampn)
    xset <- xcmsSet(files=cdffiles, phenoData=pheno)
    ## Setting the sample class; that's important for the peak grouping
    ## algorithm
    sampclass(xset) <- xset$genotype
    ## At last we define also a color for each of the two genotypes.
    genoColor <- c("#E41A1C80", "#377EB880")
    names(genoColor) <- c("KO", "WT")

    visualizeXcmsSet()
    library(shiny)
    library(shinyjs)
    runApp("../shinyHappyPeople/")
}
