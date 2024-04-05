# Load required libraries
library(shiny)
library(shinythemes)
library(DT)  # For displaying data tables

# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                
                navbarPage(
                  "Box-Cox Function App",
                  
                  # Tab 1: File Upload
                  tabPanel("Upload CSV",
                           sidebarPanel(
                              fileInput("file", "Input Patient Data (CSV)"),
                              fileInput("file", "Input AIM Variables (CSV)"),
                              downloadButton("download", "Download Transformed Data")
                           ),
                           mainPanel(
                             DTOutput("table")
                           )
                  ),
                  
                  # Tab 2: About
                  tabPanel("About",
                           tags$br(),
                           tags$p("This is a Shiny web application to transform data using the Box-cox function."),
                  )
                )
)

# Define server function
server <- function(input, output) {
  
  # Reactive function to read uploaded CSV file
  data <- reactive({
    req(input$file)  # Ensure file is uploaded
    read.csv(input$file$datapath, stringsAsFactors = FALSE)
  })
  
  # Render data table
  output$table <- renderDT({
    data()
  })
  
  # Function to prepare original CSV file for download
  output$download <- downloadHandler(
    filename = function() {
      # Construct filename based on original filename
      input$file$name
    },
    content = function(file) {
      # Ensure file is uploaded
      req(input$file)
      # Copy the uploaded CSV file to the download location
      file.copy(input$file$datapath, file)
    }
  )
}

# Create Shiny app object
shinyApp(ui = ui, server = server)
