####============================================================
##  Running the shiny app.
##
####------------------------------------------------------------
visualizeXcmsSet <- function(...){
    if(requireNamespace("shiny", quietly=TRUE)){
        if(requireNamespace("shinyjs", quietly=TRUE)){
            message("Starting the xcmsSet visualization shiny app. Use Ctrl-C to stop.")
            shiny::runApp(appDir=system.file("shinyHappyPeople",
                                             package="xcmsExtensions"), ...)
        }else{
            stop("Package 'shinyjs' not installed! Please install it to run this app.")
        }
    }else{
        stop("Package 'shiny' not installed! Please install it to run this app.")
    }
}


