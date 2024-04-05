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
  # Reactive value to store the path of the last created file
  lastCreatedFile <- reactiveVal()
  
  # # Reactive function to read uploaded CSV file
  # data <- reactive({
  #   req(input$file)  # Ensure file is uploaded
  #   read.csv(input$file$datapath, stringsAsFactors = FALSE)
  # })
  
  # Reactive function to read uploaded Patient Data CSV file
  all_data_raw <- reactive({
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
    req(all_data_raw())
    all_data_raw()
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
  
  Neg_to_Zero<-function(x){
    ifelse((x<=0.005), 0.005, x)
  }
  
  # boxcox function (I added a default lambda value of 0.5)
  #default lambda value of 0.5
  bc <- function(x,l=0.5){if(l==0) return(log(x)) else return((x^l-1)/l)}
  # inverse boxcox function
  ibc <- function(x,l=0.5){if(l==0) return(exp(x)) else return((x*l+1)^(1/l))}
  
  #var<-"CD4_CD134pCD25p"
  
  # Assuming all_data is defined as before, we create a new reactive expression for the filtered data
  all_data_filtered <- reactive({
    req(all_data_raw())  # Ensure all_data is available
    data <- all_data_raw()  # Get the data frame
    
    # Now filter the data frame
    data %>% 
      filter(Stim != "Fluzone" & Stim != "Cytostim") %>%
      filter(DonorID != "UBC_004" & DonorID != "UBC_059") %>%
      filter(Timepoint != "V1" & Timepoint != "6M" & Timepoint != "1Y") %>%
      arrange(Stim, Timepoint, DonorID)
  })
  
  
  # all_data <- all_data %>%
  #   arrange(Stim, Timepoint, DonorID)
  # all_data
  
  
  observe({
    req(input$patientData, input$AIMVariables) # Ensure both files are uploaded
    allDataValue <- all_data_filtered() # Get the current value of all_data
    variableNames <- variables()$variable # Assuming variables() returns a dataframe with a column 'variable'
    
    for (var in variableNames) {
      
      subtracted_data<-allDataValue %>% 
        dplyr::select(DonorID:Stim,{{var}}) %>% 
        tidyr::pivot_wider(names_from = Stim, values_from = {{var}}) %>% 
        dplyr::mutate(Covid_WT_sum=COVID_WT-DMSO, Covid_BA4_5_sum=COVID_BA4_5-DMSO) %>% 
        tidyr::pivot_longer(DMSO:Covid_BA4_5_sum, names_to="Stimulant", values_to = {{var}}) %>%
        dplyr::mutate(across(DonorID:Stimulant,as.factor)) %>%
        dplyr::filter(Stimulant != "COVID_WT" & Stimulant != "COVID_BA4_5") %>%
        dplyr::mutate(Stimulant=fct_recode(Stimulant, COVID_WT="Covid_WT_sum", COVID_BA4_5="Covid_BA4_5_sum")) %>% 
        dplyr::mutate(Timepoint=fct_relevel(Timepoint,"VY","V2","V3"))
      
      #use lapply to apply the Neg to Zero function to all rows
      #subtracted_data[[var]] <- lapply(subtracted_data[[var]], Neg_to_Zero)
      #subtracted_data[[var]] <- unlist(subtracted_data[[var]]) 
      filePath <- glue::glue("PREVENT_Subtracted_{var}.csv")
      subtracted_data <- arrange(subtracted_data, Stimulant, Timepoint, DonorID)
      write.csv(subtracted_data, filePath)
      lastCreatedFile(filePath)
    }
    
    for (var in variableNames) {
      
      bcxsubtracted_data<-allDataValue %>% 
        dplyr::select(DonorID:Stim,{{var}}) %>% 
        tidyr::pivot_wider(names_from = Stim, values_from = {{var}}) %>% 
        dplyr::mutate(Covid_WT_sum=ibc(bc(COVID_WT)-bc(DMSO)), Covid_BA4_5_sum=ibc(bc(COVID_BA4_5)-bc(DMSO))) %>% 
        tidyr::pivot_longer(DMSO:Covid_BA4_5_sum, names_to="Stimulant", values_to = {{var}}) %>%
        dplyr::mutate(across(DonorID:Stimulant,as.factor)) %>%
        dplyr::filter(Stimulant != "COVID_WT" & Stimulant != "COVID_BA4_5") %>%
        dplyr::mutate(Stimulant=fct_recode(Stimulant, COVID_WT="Covid_WT_sum", COVID_BA4_5="Covid_BA4_5_sum")) %>% 
        dplyr::mutate(Timepoint=fct_relevel(Timepoint,"VY","V2","V3"))
      
      #use lapply to apply the Neg to Zero function to all rows
      #bcxsubtracted_data[[var]] <- lapply(bcxsubtracted_data[[var]], Neg_to_Zero)
      #bcxsubtracted_data[[var]] <- unlist(bcxsubtracted_data[[var]]) 
      filePath <- glue::glue("PREVENT_Boxcox_{var}.csv")
      bcxsubtracted_data <- arrange(bcxsubtracted_data, Stimulant, Timepoint, DonorID)
      write.csv(bcxsubtracted_data, filePath)
      lastCreatedFile(filePath)
    }
  })
  
  
  
  # ---------------------------- Box-Cox Calculation End -------------------------------------------
  
  
  # Function to prepare original CSV file for download
  output$download <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "transformed-data.csv", sep = "_") # Provide a meaningful default filename
    },
    content = function(file) {
      req(lastCreatedFile())  # Ensure there is a file to download
      file.copy(lastCreatedFile(), file)  # Copy the last created file to the download location
    }
  )
}

# Create Shiny app object
shinyApp(ui = ui, server = server)
