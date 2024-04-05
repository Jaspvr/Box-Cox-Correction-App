# Load required libraries for calculation
#Note: tidyverse contains ggplot2, dplyr, tidyr, stringr
library(openxlsx)
library(IDPmisc)
library(glue)
library(tidyverse)
library(conflicted)

# Load required libraries for UI
library(shiny)
library(shinythemes)
library(DT)  # For displaying data tables

# Resolve conflicts
conflicts_prefer(dplyr::filter)

# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                
                navbarPage(
                  "Box-Cox Function App",
                  
                  # Tab 1: File Upload
                  tabPanel("Upload CSV",
                           sidebarPanel(
                              fileInput("patientData", "Input Patient Data (CSV)"),
                              fileInput("AIMVariables", "Input AIM Variables (CSV)"),
                              downloadButton("download", "Download Transformed Data")
                           ),
                           mainPanel(
                             DTOutput("table1"),
                             DTOutput("variablesTable")
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
  
  # # Reactive function to read uploaded CSV file
  # data <- reactive({
  #   req(input$file)  # Ensure file is uploaded
  #   read.csv(input$file$datapath, stringsAsFactors = FALSE)
  # })
  
  # Reactive function to read uploaded Patient Data CSV file
  all_data <- reactive({
    req(input$patientData)  # Ensure file is uploaded
    read.csv(input$patientData$datapath, stringsAsFactors = FALSE)
  })
  
  # Reactive function to read uploaded AIM Variables CSV file
  variables <- reactive({
    req(input$AIMVariables)  # Ensure file is uploaded
    vars <- read.csv(input$AIMVariables$datapath, stringsAsFactors = FALSE) %>% 
      pull()
    data.frame(variable = vars)
  })
  
  # Render data table for all_data
  output$table1 <- renderDT({
    req(all_data())
    all_data()
  })
  
  # Render data table for variables
  # Ensure this output ID is unique and matches an output in your UI
  output$variablesTable <- renderDT({
    req(variables())
    variables()
  })
  
  
  # ---------------------------- Box-Cox Calculation Start --------------------------------------------
  
  # Create last set of graphs without any parameters highlighted (variables and data read in above)
  stimulants<-c("DMSO","Fluzone", "COVID_WT", "COVID_BA4_5")
  
  
  
  
  
  
  
  
  
  
  
  
  # ---------------------------- Box-Cox Calculation End -------------------------------------------
  
  
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
