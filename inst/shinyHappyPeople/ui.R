####============================================================
##  shinyUI
##
##  The shiny UI
####------------------------------------------------------------
shinyUI(fluidPage(
    titlePanel("xcmsSet visualization"),

    shiny::fluidRow(
        shiny::column(3,
               wellPanel(
                   helpText(paste0("Select an xcmsSet object to visualize;",
                                   " hit the 'Update view' button to update the plots.")),
                   ## Drop down menu for xcmsSet object selection.
                   uiOutput("xcmsSetObjects"),
                   ## Range slider.
                   uiOutput("rtrange"),
                   uiOutput("mzrange"),
                   selectInput("nbin", label="Number of bins",
                               choices=c(20, 50, 100, 200, 400, 800, NULL), selected=100),
                   submitButton("Update view"),
                   actionButton("closeButton", "Return & close")
                   ## checkboxInput("chooseCol", "Choose colors"),
                   ## conditionalPanel(
                   ##     condition = "input.chooseCol == true",
                   ##     uiOutput("colorPicker")
                   ## )
                   ),
               wellPanel(
                   strong("Choose sample color"),
                   uiOutput("colorPicker")
               )
               ),
        shiny::column(9,
               wellPanel(
                   h4("Chromatogram"),
                   plotOutput("chromPlot")
               ),
               wellPanel(
                   h4("Spectrum"),
                   plotOutput("specPlot")
               )
               )
    )

