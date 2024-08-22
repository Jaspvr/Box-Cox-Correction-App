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
library(shinyalert)  # For displaying popup alerts

# Resolve conflicts
conflicts_prefer(dplyr::filter)

# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
    navbarPage(
      "Box-Cox Correction App",
      
      # Tab 1: File Upload, user specification based on their data
      tabPanel("Upload CSV",
         sidebarPanel(
           numericInput("lambda", "Input Lambda value", value = 0.5),
           textInput("unstimulated", "Unstimulated Parameter", value = "unstimulated"),
           textAreaInput("stimulants", "Stimulants (comma separated)", value = "CMV_protein, CMV_peptides, CytoStim, Infanrix, COVID_S_Ag"),
           fileInput("patientData", "Input Patient Data (CSV)"),
           selectInput("AIMVariables", "AIM Variables (from patient data headers)", choices = NULL, multiple = TRUE),
           selectInput("grouping_columns", "Grouping Columns (from patient data headers)", choices = NULL, multiple = TRUE),
           selectInput("stim_column", "Stimulant Column (from patient data headers)", choices = NULL),
           downloadButton("download", "Download Transformed Data")
         ),
         mainPanel(
           DTOutput("table1"),
         )
      ),
      
      # Tab 2: About section
      tabPanel("About",
         tags$br(),
         tags$p("Analysis of T cell activation-induced marker (AIM) assay data requires normalization of AIM+ cell frequencies to background AIM+ frequencies in an unstimulated control. 
Subtracting or dividing by the unstimulated control each have specific disadvantages and can amplify technical variability in the assay. 
The Box-Cox correction is an innovative method with features of both division and linear subtraction, allowing a more sophisticated correction for unstimulated AIM+ cell frequencies that better aligns with the mathematical properties of AIM datasets and reduces technical variability."), 
         tags$br(),
         tags$p("To take advantage of the Box-Cox correction, upload your full AIM dataset and the set of variables to be corrected. The Box-Cox Correction App will immediately return the corrected values which are then ready for data display or statistical analysis.")
      ),
    )
)

# Define server function
server <- function(input, output, session) {
  # Reactive value to store the path of the last created file
  lastCreatedFile <- reactiveVal()
  
  # Reactive function to read uploaded Patient Data CSV file
  all_data_raw <- reactive({
    req(input$patientData)  # Ensure file is uploaded
    tryCatch({
      read.csv(input$patientData$datapath, stringsAsFactors = FALSE)
    }, error = function(e) {
      shinyalert("Error", paste("Error reading patient data:", e$message), type = "error")
      return(NULL)
    })
  })
  
  # Update the selectInput choices based on the uploaded data
  observe({
    req(all_data_raw())
    updateSelectInput(session, "grouping_columns", choices = names(all_data_raw()))
    updateSelectInput(session, "AIMVariables", choices = names(all_data_raw()))
    updateSelectInput(session, "stim_column", choices = names(all_data_raw()))
  })
  
  # Render data table for all_data
  output$table1 <- renderDT({
    req(all_data_raw())
    all_data_raw()
  })
  
  
  # ---------------------------- Box-Cox Calculation Start --------------------------------------------
  
  Neg_to_Zero<-function(x){
    ifelse((x<=0.005), 0.005, x)
  }
  
  # Default lambda value of 0.5
  reactiveLambda <- reactive({
    input$lambda  # This will be a number because of numericInput
  })
  
  # Box Cox function
  bc <- function(x) {
    l <- reactiveLambda()
    if (l == 0) {
      return(log(x))
    } else {
      return((x^l - 1) / l)
    }
  }
  
  # Inverse Box Cox function
  ibc <- function(x) {
    l <- reactiveLambda()
    if (l == 0) {
      return(exp(x))
    } else {
      return((x * l + 1)^(1 / l))
    }
  }
  
  # Assuming all_data is defined as before, we create a new reactive expression for the filtered data
  all_data_filtered <- reactive({
    req(all_data_raw())  # Ensure all_data is available
    data <- all_data_raw()  # Get the data frame
    
    # These are the columns that are used to make groups based on rows having the same value for each of these columns
    grouping_columns <- input$grouping_columns
    
    tryCatch({
      data %>% 
        arrange(across(all_of(grouping_columns))) # Arrange by inputted identifier columns instead
    }, error = function(e) {
      shinyalert("Error", paste("Error grouping data:", e$message), type = "error")
      return(NULL)
    })
  })
  
  # Function to prepare original CSV file for download
  output$download <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "transformed-data.csv", sep = "_")  # Provide a meaningful default filename
    },
    content = function(file) {
      # Ensure Data and variables are inputted
      req(input$patientData, input$AIMVariables)
      
      allDataValue <- all_data_filtered()  # Get the current value of all_data
      variableNames <- input$AIMVariables
      stim_column <- input$stim_column
      if (is.null(allDataValue)) {
        shinyalert("Error", "Patient data is empty or not properly loaded", type = "error")
        return()
      }
      if (is.null(variableNames)) {
        shinyalert("Error", "Patient data is empty or not properly loaded", type = "error")
        return()
      }
      
      # Get user inputted values for stimulants, unstimulated parameter, and for the grouping columns
      stimulants <- strsplit(input$stimulants, ",\\s*")[[1]]
      unstimulated_parameter <- input$unstimulated
      grouping_columns <- input$grouping_columns
      # grouping_columns <- strsplit(input$grouping_columns, ",\\s*")[[1]]
      # print(grouping_columns)
      
      # Error handling:
      if (length(grouping_columns) == 0) {
        shinyalert("Error", "Grouping columns are not properly selected.", type = "error")
        return()
      }
      
      # Filter the variable names to only those present in the data
      # variableNames <- variableNames[variableNames %in% names(allDataValue)]
      # if (length(variableNames) == 0) {
      #   shinyalert("Error", "None of the selected AIM variables are present in the patient data.", type = "error")
      #   return()
      # }
      
      # Group by matching grouping columns
      tryCatch({
        grouped_data <- allDataValue %>%
          group_by(across(all_of(grouping_columns)))
      }, error = function(e) {
        shinyalert("Error", paste("Error grouping data:", e$message), type = "error")
        return(NULL)
      })
      
      # Function to apply Box-Cox transformation, subtract unstimulated control, and apply inverse Box-Cox transformation
      transform_and_subtract <- function(df, unstim_var) {
        unstimulated <- df %>% filter(!!sym(stim_column) == unstim_var)

        if (nrow(unstimulated) == 0) {
          return(df)
        }

        unstim_values <- unstimulated[variableNames]


        df %>%
          mutate(across(all_of(variableNames), ~ ibc(bc(.x) - bc(unstim_values[[cur_column()]]))))
      }

      # Apply the transformations within each group
      tryCatch({
        transformed_data <- grouped_data %>%
          group_modify(~ transform_and_subtract(.x, unstimulated_parameter)) %>%
          ungroup()

        temp_file <- tempfile(fileext = ".csv")
        write.csv(transformed_data, temp_file, row.names = FALSE)
        lastCreatedFile(temp_file)
        file.copy(lastCreatedFile(), file)  # Copy the last created file to the download location
      }, error = function(e) {
        shinyalert("Error", paste("Error during transformation:", e$message), type = "error")
        return(NULL)
      })
    }
  )
}

# Create Shiny app object
shinyApp(ui = ui, server = server)