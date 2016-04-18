####============================================================
##  .getXsetFromWs
##
##  Get all the xcmsSet objects in the present working space.
####------------------------------------------------------------
.getXsetFromWs <- function(){
    message(".getXsetFromWs...", appendLF=FALSE)
    allVars <- ls(envir=globalenv())
    xses <- lapply(allVars, function(z){
        tmp <- get(z, envir=globalenv())
        if(is(tmp, "xcmsSet"))
            return(z)
    })
    message("OK")
    return(unlist(xses, use.names=FALSE))
}

####============================================================
##  .getMSsliceForXcmsSet
##
##  Takes the name of an xcmsSet object, extracts the full MSslice
##  object and stores that in a temporary variable in the current
##  workspace. If such a temporary variable does already exist it
##  returns the existing one.
####------------------------------------------------------------
.getMSsliceForXcmsSet <- function(x){
    message(".getMSsliceForXcmsSet...", appendLF=FALSE)
    varName <- paste0(".XCMSEXT_TMP_", x)
    if(any(ls(envir=globalenv(), all.names=TRUE) == varName)){
        message("OK")
        return(get(varName, envir=globalenv()))
    }else{
        message("Extracting MSslice object from ", x, "...", appendLF=FALSE)
        suppressMessages(
            mss <- msSlice(get(x, envir=globalenv()))
        )
        ## Remove all potentially existing other temporary vars.
        .cleanupWs()
        message("OK")
        assign(varName, mss, envir=globalenv())
        message("OK")
        return(mss)
    }
}

####============================================================
##  .cleanupWs
##
##  Clean the workspace up, i.e. remove temporary variables.
####------------------------------------------------------------
.cleanupWs <- function(){
    message(".cleanupWs...", appendLF=FALSE)
    vars <- ls(envir=globalenv(), all.names=TRUE)
    idx <- grep(vars, pattern="^.XCMSEXT_TMP_")
    if(length(idx) > 0){
        rm(list=vars[idx], envir=globalenv())
    }
    message("OK")
}

####============================================================
##  textInputSmall
##
####------------------------------------------------------------
textInputSmall <- function(inputId, label, value="", ...){
    div(style="display:inline-block",
        tags$label(label, `for`=inputId),
        tags$input(id=inputId, type="text", value=value, ...))
}


####============================================================
##  shiny server.
##
##
####------------------------------------------------------------
shinyServer(
    function(input, output){
        ## Idea: store the MSslice object somewhere (global environment)
        ## and always "only" subset that one if there is a refresh.
        ## That should help avoiding to always completely reloading the data.
        ## Need then also to remove other objects when there is an object change.
        ##
        ## Set the list of available xcmsSet objects.
        xses <- c("- select -", .getXsetFromWs())
        if(length(xses) == 1)
            stop("\nNo xcmsSet objects found in the current workspace!\n")
        ## Extract the MS slice object from the selected xset given
        ## the settings; if now range is available, just return all.
        msSliceData <- reactive({
            if(length(input$xcmsObject) > 0){
                if(input$xcmsObject != "- select -"){
                    mss <- .getMSsliceForXcmsSet(input$xcmsObject)
                    ## Eventually subset the slice object.
                    ## Note: if we're doing this here, we're also changing the slider
                    ##       range; if we're doing it below, the ranges will remain the
                    ##       global range.
                    ## if(!is.null(mss)){
                    ##     ## Get data from other inputs:
                    ##     rtr <- input$rtrange
                    ##     message("Got rtrange: (", rtr[1], ", ", rtr[2], ")")
                    ##     if(rtr[1] > 0){
                    ##         mss <- subset(mss, rtrange=rtr)
                    ##     }
                    ##     mzr <- input$mzrange
                    ##     message("Got mzrange: (", mzr[1], ", ", mzr[2], ")")
                    ##     if(mzr[1] > 0){
                    ##         mss <- subset(mss, mzrange=mzr)
                    ##     }
                    ## }
                    return(mss)
                }
            }
            return(NULL)
        })
        ## Drop down menu for the selection of the xcmsSet object..
        output$xcmsSetObjects <-  renderUI(
            selectInput("xcmsObject", "Select the xcmsSet object", xses)
        )
        ## Slider to select the retention time range.
        output$rtrange <- renderUI({
            mss <- msSliceData()
            if(!is.null(mss)){
                rngs <- round(rtrange(mss))
                ## Do a rtrange stepping of 1 (sec); eventually do it based on
                ## the rtrange? e.g. diff(rtrange) / 100?
                sliderInput("rtrange", "rtrange", rngs[1], rngs[2], value=c(rngs[1], rngs[2]),
                            step=1)
            }else{
                sliderInput("rtrange", "rtrange", -10, -5, value=c(-10, -5), step=0.1)
            }
        })
        ## Slider to select the M/Z range.
        output$mzrange <- renderUI({
            mss <- msSliceData()
            if(!is.null(mss)){
                rngs <- round(mzrange(mss))
                sliderInput("mzrange", "rtrange", rngs[1], rngs[2], value=c(rngs[1], rngs[2]),
                            step=0.1)
            }else{
                sliderInput("mzrange", "rtrange", -10, -5, value=c(-10, -5), step=0.1)
            }
        })
        ## Dynamic color selection / sample.
        output$colorPicker <- renderUI({
            mss <- msSliceData()
            if(!is.null(mss)){
                sampIdx <- 1:length(mss@data)
                sampNames <- as.character(sampIdx)
                if(!is.null(names(mss))){
                    sampNames <- names(mss)
                }
                lapply(sampIdx, function(z){
                    ## textInput(paste0("sampleColor", z), label=sampNames[z], value="black",
                    ##           width="30px")
                    ## textInputSmall(paste0("sampleColor", z), label=sampNames[z], value="black",
                    ##                class="input-mini")
                    ## That requires shinyjs!
                    colourInput(inputId=paste0("sampleColor", z), label=sampNames[z],
                                value="black", allowTransparent=FALSE)
                })
            }
        })

        ## Results panel stuff.
        output$selectedObjectName <- renderText({
            ifelse(input$xcmsObject == "- select -", yes="", no=input$xcmsObject)
        })
        ## Chromatogram
        output$chromPlot <- renderPlot({
            ## Making a progress...
            progress <- shiny::Progress$new()
            on.exit(progress$close())
            progress$set(message="Getting data", value=0)
            mss <- msSliceData()
            progress$set(message="Subsetting data", value=1/3)
            if(!is.null(mss)){
                ## Subset
                rtr <- input$rtrange
                message("Got rtrange: (", rtr[1], ", ", rtr[2], ")")
                if(rtr[1] > 0){
                    mss <- subset(mss, rtrange=rtr)
                }
                mzr <- input$mzrange
                message("Got mzrange: (", mzr[1], ", ", mzr[2], ")")
                if(mzr[1] > 0){
                            mss <- subset(mss, mzrange=mzr)
                }
                ## Get the colors...
                sampColors <- lapply(1:length(mss@data), function(z){
                    message("Sample ", z, ":", appendLF=FALSE)
                    theCol <- input[[paste0("sampleColor", z)]]
                    if(is.null(theCol))
                        theCol <- "#000000"
                    message(" ", theCol)
                    return(paste0(theCol, "80"))
                })
                ## Plot chromatogram
                progress$set(message="Subsetting data", value=2/3)
                message("Plotting the chromatogram...", appendLF=FALSE)
                plotChromatogram(mss, type="l", col=unlist(sampColors, use.names=FALSE))
                message("OK")
            }else{
                ##par()
            }
        })
        ## Spectrum
        output$specPlot <- renderPlot({
            ## Making a progress...
            progress <- shiny::Progress$new()
            on.exit(progress$close())
            progress$set(message="Getting data", value=0)
            mss <- msSliceData()
            progress$set(message="Subsetting data", value=1/3)
            if(!is.null(mss)){
                ## Subset
                rtr <- input$rtrange
                message("Got rtrange: (", rtr[1], ", ", rtr[2], ")")
                if(rtr[1] > 0){
                    mss <- subset(mss, rtrange=rtr)
                }
                mzr <- input$mzrange
                message("Got mzrange: (", mzr[1], ", ", mzr[2], ")")
                if(mzr[1] > 0){
                            mss <- subset(mss, mzrange=mzr)
                }
                ## Get the colors...
                sampColors <- lapply(1:length(mss@data), function(z){
                    message("Sample ", z, ":", appendLF=FALSE)
                    theCol <- input[[paste0("sampleColor", z)]]
                    if(is.null(theCol))
                        theCol <- "#000000"
                    message(" ", theCol)
                    return(paste0(theCol, "80"))
                })
                ## Plot spectrum
                message("Plotting the spectrum...", appendLF=FALSE)
                message("Plotting the spectrum...", appendLF=FALSE)
                plotSpectrum(mss, type="l", col=unlist(sampColors, use.names=FALSE))
                message("OK")
            }else{
                ##par()
            }
        })
        ## Check if we hit the close button; if so return.
        observe({
            if(input$closeButton > 0){
                message("See you next time.")
                ## OK, now, get the MSsliceObject and return that.
                mss <- msSliceData()
                if(!is.null(mss)){
                    ## Subset
                    rtr <- input$rtrange
                    message("Got rtrange: (", rtr[1], ", ", rtr[2], ")")
                    if(rtr[1] > 0){
                        mss <- subset(mss, rtrange=rtr)
                    }
                    mzr <- input$mzrange
                    message("Got mzrange: (", mzr[1], ", ", mzr[2], ")")
                    if(mzr[1] > 0){
                        mss <- subset(mss, mzrange=mzr)
                    }
                }
                .cleanupWs()
                stopApp(mss)
            }
        })
    })


#### OLD STUFF FROM ensembldb BELOW
##
## ## list all packages...
## packs <- installed.packages()
## epacks <- packs[grep(packs, pattern="^Ens")]

## ## library(EnsDb.Hsapiens.v75)
## ## edb <- EnsDb.Hsapiens.v75

## TheFilter <- function(input){
##     Cond <- input$condition
##     ## check if we've got something to split...
##     Vals <- input$geneName
##     ## check if we've got ,
##     if(length(grep(Vals, pattern=",")) > 0){
##         ## don't want whitespaces here...
##         Vals <- gsub(Vals, pattern=" ", replacement="", fixed=TRUE)
##         Vals <- unlist(strsplit(Vals, split=","))
##     }
##     if(length(grep(Vals, pattern=" ", fixed=TRUE)) > 0){
##         Vals <- unlist(strsplit(Vals, split=" ", fixed=TRUE))
##     }
##     if(input$type=="Gene name"){
##         return(GenenameFilter(Vals, condition=Cond))
##     }
##     if(input$type=="Chrom name"){
##         return(SeqnameFilter(Vals, condition=Cond))
##     }
##     if(input$type=="Gene biotype"){
##         return(GenebiotypeFilter(Vals, condition=Cond))
##     }
##     if(input$type=="Tx biotype"){
##         return(TxbiotypeFilter(Vals, condition=Cond))
##     }
## }

## ## checkSelectedPackage <- function(input){
## ##     if(is.null(input$package)){
## ##         return(FALSE)
## ##     }else{
## ##         require(input$package, character.only=TRUE)
## ##         message("Assigning ", input$package, " to variable edb.")
## ##         assign("edb", get(input$package), envir=globalenv())
## ##         return(TRUE)
## ##     }
## ## }

## ## Based on the given EnsDb package name it loads the library and returns
## ## the object.
## getEdb <- function(x){
##     require(x, character.only=TRUE)
##     return(get(x))
## }

## ## Define server logic required to draw a histogram
## shinyServer(function(input, output) {

##     ## Generate the select field for the package...
##     output$packages <- renderUI(
##         selectInput("package", "Select installed EnsDb package", as.list(epacks))
##     )

##     selectedPackage <- reactive({
##         names(epacks) <- epacks
##         ## epacks <- sapply(epacks, as.symbol)
##         ## load the package.
##         if(length(input$package) > 0){
##             require(input$package, character.only=TRUE)
##             ## Actually, should be enough to just return input$package...
##             ##return(switch(input$package, epacks))
##             return(getEdb(input$package))
##         }else{
##             return(NULL)
##         }
##     })

##     ## Metadata infos
##     output$metadata_organism <- renderText({
##         edb <- selectedPackage()
##         if(!is.null(edb)){
##             ## db <- getEdb(edb)
##             paste0("Organism: ", organism(edb))
##         }

##     })
##     output$metadata_ensembl <- renderText({
##         edb <- selectedPackage()
##         if(!is.null(edb)){
##             ## db <- getEdb(edb)
##             md <- metadata(edb)
##             rownames(md) <- md$name
##             paste0("Ensembl version: ", md["ensembl_version", "value"])
##         }
##     })
##     output$metadata_genome <- renderText({
##         edb <- selectedPackage()
##         if(!is.null(edb)){
##             ## db <- getEdb(edb)
##             md <- metadata(edb)
##             rownames(md) <- md$name
##             paste0("Genome build: ", md["genome_build", "value"])
##         }
##     })

##     output$genename <- renderText({
##         if(length(input$geneName) > 0){
##             input$geneName
##         }else{
##             return()
##         }
##     })
##     ## That's the actual queries for genes, transcripts and exons...
##     output$Genes <- renderDataTable({
##         ## if(!checkSelectedPackage(input))
##         ##     return()
##         if(length(input$package) == 0)
##             return
##         if(!is.na(input$geneName) & length(input$geneName) > 0 & input$geneName!=""){
##             edb <- selectedPackage()
##             res <- genes(edb, filter=TheFilter(input),
##                          return.type="data.frame")
##             assign(".ENS_TMP_RES", res, envir=globalenv())
##             return(res)
##         }
##     })
##     output$Transcripts <- renderDataTable({
##         if(length(input$package) == 0)
##             return
##         if(!is.na(input$geneName) & length(input$geneName) > 0 & input$geneName!=""){
##             edb <- selectedPackage()
##             res <- transcripts(edb, filter=TheFilter(input),
##                          return.type="data.frame")
##             assign(".ENS_TMP_RES", res, envir=globalenv())
##             return(res)
##         }
##     })
##     output$Exons <- renderDataTable({
##         if(length(input$package) == 0)
##             return
##         if(!is.na(input$geneName) & length(input$geneName) > 0 & input$geneName!=""){
##             edb <- selectedPackage()
##             res <- exons(edb, filter=TheFilter(input),
##                          return.type="data.frame")
##             assign(".ENS_TMP_RES", res, envir=globalenv())
##             return(res)
##         }
##     })
##     observe({
##         if(input$closeButton > 0){
##             ## OK, now, gather all the data and return it in the selected format.
##             edb <- selectedPackage()
##             resType <- input$returnType
##             resTab <- input$resultTab
##             res <- NULL
##             ## If result type is data.frame we just return what we've got.
##             if(resType == "data.frame"){
##                 res <- get(".ENS_TMP_RES")
##             }else{
##                 ## Otherwise we have to fetch a little bit more data, thus, we perform the
##                 ## query again and return it as GRanges.
##                 if(resTab == "Genes")
##                     res <- genes(edb, filter=TheFilter(input), return.type="GRanges")
##                 if(resTab == "Transcripts")
##                     res <- transcripts(edb,filter=TheFilter(input), return.type="GRanges")
##                 if(resTab == "Exons")
##                     res <- exons(edb,filter=TheFilter(input), return.type="GRanges")
##             }
##             rm(".ENS_TMP_RES", envir=globalenv())
##             stopApp(res)
##         }
##     })
## })



## ## ## Define server logic required to draw a histogram
## ## shinyServer(function(input, output) {

## ##     ## generate the select field for the package...
## ##     output$packages <- renderUI(
## ##         selectInput("package", "Select EnsDb package", as.list(epacks))
## ##     )

## ##     ## generating metadata info.
## ##     output$metadata_organism <- renderText({
## ##         if(!checkSelectedPackage(input))
## ##             return()
## ##         paste0("Organism: ", organism(edb))
## ##     })
## ##     output$metadata_ensembl <- renderText({
## ##         if(!checkSelectedPackage(input))
## ##             return()
## ##         md <- metadata(edb)
## ##         rownames(md) <- md$name
## ##         paste0("Ensembl version: ", md["ensembl_version", "value"])
## ##     })
## ##     output$metadata_genome <- renderText({
## ##         if(!checkSelectedPackage(input))
## ##             return()
## ##         md <- metadata(edb)
## ##         rownames(md) <- md$name
## ##         paste0("Genome build: ", md["genome_build", "value"])
## ##     })
## ##     ## output$genename <- renderText({
## ##     ##     if(!checkSelectedPackage(input))
## ##     ##         return()
## ##     ##     input$geneName
## ##     ## })
## ##     ## That's the actual queries for genes, transcripts and exons...
## ##     output$Genes <- renderDataTable({
## ##         if(!checkSelectedPackage(input))
## ##             return()
## ##         if(!is.na(input$geneName) & length(input$geneName) > 0 & input$geneName!=""){
## ##             res <- genes(edb, filter=TheFilter(input),
## ##                          return.type="data.frame")
## ##             return(res)
## ##         }
## ##     })
## ##     output$Transcripts <- renderDataTable({
## ##         if(!checkSelectedPackage(input))
## ##             return()
## ##         if(!is.na(input$geneName) & length(input$geneName) > 0 & input$geneName!=""){
## ##             res <- transcripts(edb, filter=TheFilter(input),
## ##                          return.type="data.frame")
## ##             return(res)
## ##         }
## ##     })
## ##     output$Exons <- renderDataTable({
## ##         if(!checkSelectedPackage(input))
## ##             return()
## ##         if(!is.na(input$geneName) & length(input$geneName) > 0 & input$geneName!=""){
## ##             res <- exons(edb, filter=TheFilter(input),
## ##                          return.type="data.frame")
## ##             return(res)
## ##         }
## ##     })
## ##     ## observe({
## ##     ##     if(input$Close > 0){
## ##     ##         stopApp("AAARGHHH")
## ##     ##     }
## ##     ## })
## ## })


