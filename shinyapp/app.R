library(shiny)
library(shinyBS)
library(R2easyR)
library(ggplot2)

## devmode shows a more comprehensive error message
# devmode(TRUE)

options(shiny.port = 8080)
options(shiny.host = "0.0.0.0")


workdir <- getwd()
#palettes
palettes <- r2easyR.palettes()


ui <- fluidPage(
  titlePanel("R2easyR"),
  sidebarLayout(
    sidebarPanel(
      bsCollapse(
        id = "collapsePanel",
        open = "Step 1",
        bsCollapsePanel("Step 1",
          h1("Step 1"),
          fileInput("file", "Choose CSV File",
                    multiple = FALSE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          tags$hr(),
          checkboxInput("header", "Header", TRUE),
          radioButtons("sep", "Separator",
                       choices = c(Comma = ",",
                                   Semicolon = ";",
                                   Tab = "\t"),
                       selected = ","),
          radioButtons("quote", "Quote",
                       choices = c(None = "",
                                   "Double Quote" = "\"",
                                   "Single Quote" = "'"),
                       selected = "\""),
          tags$hr(),
          radioButtons("disp", "Display",
                       choices = c(Head = "head",
                                   All = "all"),
                       selected = "head")
        ),
        bsCollapsePanel("Step 2",
                        h1("Step 2"),
                        tags$hr(),
                        ## dropdown menu for selecting the palette
                        selectInput("palette", "Select a palette",
                                    choices = names(palettes)),
                        ## plot name input box
                        textInput("plotname", "Enter the plot name", value = "plot"),

                        ## action button, only active when the plot name is not empty and palette is selected
                        actionButton("plot", "Plot", class = "btn-primary", icon = icon("bar-chart-o"))

        )
      )
    ),
    mainPanel(
      tableOutput("contents")
    )
  )
)

server <- function(input, output, session) {

  fileContent <- reactive({
    file <- input$file
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    
    validate(need(ext == "csv", "Please upload a csv file"))
    
    read.csv(file$datapath,
             header = input$header,
             sep = input$sep,
             quote = input$quote)
  })
  
  output$contents <- renderTable({
   req(fileContent())
    if(input$disp == "head") {
      head(fileContent())
    }
    else {
      fileContent()
    }
  },
  striped = FALSE,
  width = "auto")
  
  observeEvent(input$file, {
    updateCollapse(session, "collapsePanel", "Step 2")
  })

  observeEvent(input$plot, {
    req(input$plotname)
    req(input$palette)
    df <- fileContent()
    df <- r2easyR.color(df, palettes[[input$palette]], abs_reactivity_threshold = 0.1)
    r2easyR.write("demo", df, input$plotname, colors = "circles")
    r2easyR.stem_editor("demo.sto")
    name <- paste0(stringi::stri_rand_strings(1,10), ".svg")
    system(paste0("r2r --disable-usage-warning demo.r2r_meta www/",name))
    #load local file demo.svg plot svg plot
    output$contents <- renderUI({
      tags$img(src = name, width = "100%", height = "100%")
    })
  })
  
}

shinyApp(ui, server)
