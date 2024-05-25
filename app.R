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
    "Box-Cox Correction App",
    
    # Tab 1: File Upload
    tabPanel("Upload CSV",
     sidebarPanel(
       numericInput("lambda", "Input Lambda value", value = 0.5),
       textInput("unstimulated", "Unstimulated Parameter", value = "DMSO"),
       textAreaInput("stimulants", "Stimulants (comma separated)", value = "Fluzone, COVID_WT, COVID_BA4_5, Cytostim"),
       textAreaInput("grouping_columns", "Grouping Columns (comma separated)", value = "Timepoint, DonorID"),
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
     tags$p("Analysis of T cell activation-induced marker (AIM) assay data requires normalization of AIM+ cell frequencies to background AIM+ frequencies in an unstimulated control. 
     Subtracting or dividing by the unstimulated control each have specific disadvantages and can amplify technical variability in the assay. 
     The Box-Cox correction is an innovative method with features of both division and linear subtraction, allowing a more sophisticated correction for unstimulated AIM+ cell frequencies that better aligns with the mathematical properties of AIM datasets and reduces technical variability."), 
     tags$br(),
     tags$p("To take advantage of the Box-Cox correction, upload your full AIM dataset and the set of variables to be corrected. The Box-Cox Correction App will immediately return the corrected values which are then ready for data display or statistical analysis.")
    ),
  )
)

# Define server function
server <- function(input, output) {
  # Reactive value to store the path of the last created file
  lastCreatedFile <- reactiveVal()
  
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
  
  # These will be eventually user defined
  # stimulants<-c("Fluzone", "COVID_WT", "COVID_BA4_5", "Cytostim")
  # unstimulated_parameter<-"DMSO"
  
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
  
  # Stim is special, this is what we use to do operations on groups
  # grouping_columns <- c("Timepoint", "DonorID")
  
  # Assuming all_data is defined as before, we create a new reactive expression for the filtered data
  all_data_filtered <- reactive({
    req(all_data_raw())  # Ensure all_data is available
    data <- all_data_raw()  # Get the data frame
    
    grouping_columns <- strsplit(input$grouping_columns, ",\\s*")[[1]]
    
    # Now arrange the data
    data %>% 
      arrange(across(all_of(grouping_columns))) # Arrange by inputted identifier columns instead
  })
  
  observe({
    req(input$patientData, input$AIMVariables) # Ensure both files are uploaded
    allDataValue <- all_data_filtered() # Get the current value of all_data
    variableNames <- variables()$variable # Assuming variables() returns a dataframe with a column 'variable'
    
    # Get user inputted values for stimulants, unstimulated parameter, and for the grouping columns
    stimulants <- strsplit(input$stimulants, ",\\s*")[[1]]
    unstimulated_parameter <- input$unstimulated
    grouping_columns <- strsplit(input$grouping_columns, ",\\s*")[[1]]
    
    # Group by matching grouping columns
    grouped_data <- allDataValue %>%
      group_by(across(all_of(grouping_columns)))
    
    # Function to apply Box-Cox transformation, subtract unstimulated control, and apply inverse Box-Cox transformation
    transform_and_subtract <- function(df, unstim_var) {
      unstimulated <- df %>% filter(Stim == unstim_var)
      if (nrow(unstimulated) == 0) return(df)
      
      unstim_values <- unstimulated[variableNames]
      
      df %>%
        mutate(across(all_of(variableNames), ~ ibc(bc(.x) - bc(unstim_values[[cur_column()]]))))
    }
    
    # Apply the transformations within each group
    transformed_data <- grouped_data %>%
      group_modify(~ transform_and_subtract(.x, unstimulated_parameter)) %>%
      ungroup()
    
    # # Function to subtract unstimulated control values
    # subtract_unstimulated <- function(df, unstim_var) {
    #   unstimulated <- df %>% filter(Stim == unstim_var)
    #   if (nrow(unstimulated) == 0) return(df)
    #   
    #   unstim_values <- unstimulated[variableNames]
    #   
    #   df %>% 
    #     mutate(across(all_of(variableNames), ~ .x - unstim_values[[cur_column()]]))
    # }
    # 
    # # Apply the subtraction within each group
    # transformed_data <- grouped_data %>%
    #   group_modify(~ subtract_unstimulated(.x, unstimulated_parameter)) %>%
    #   ungroup()
    
    
    temp_file <- tempfile(fileext = ".csv")
    write.csv(transformed_data, temp_file, row.names = FALSE)
    # write.csv(grouped_data, temp_file, row.names = FALSE)
    lastCreatedFile(temp_file)
    
    
    
    
    
    
    # for (var in variableNames) {
    #   
    #   # Process data to subtract the unstimulated value from other stimulants
    #   processedData <- allDataValue %>%
    #     group_data() %>%  # Use the defined function for grouping
    #     mutate(
    #       ReferenceValue = var[Stimulant == unstimulated_parameter]  # Create a reference value column with unstimulated value
    #     ) %>%
    #     mutate(
    #       AdjustedValue = if_else(Stimulant == unstimulated_parameter, Value, Value - ReferenceValue)  # Adjust values
    #     ) %>%
    #     ungroup()  # Ungroup data for further processing if necessary
    #   
    #   # Optionally filter out the original value column if only adjusted values are needed
    #   processedData <- processedData %>%
    #     select(-Value) %>%
    #     rename(Value = AdjustedValue)
    #   
    #   # Merge processed data into combinedData or simply assign it if combinedData is not used elsewhere
    #   if (is.null(combinedData)) {
    #     combinedData <- processedData
    #   } else {
    #     combinedData <- bind_rows(combinedData, processedData)
    #   }
    # }
    
    
    
    # for (var in variableNames) {
    #   # Get the Unstimulated parameter value for the 
    #   
    #   for (stimulant in stimulants) {
    #     cat(var, stimulant, "\n")
    #   }
    # }
    
    # write.csv(combinedData, "PREVENT_Boxcox_Combined.csv")
    # lastCreatedFile("PREVENT_Boxcox_Combined.csv")
  })
  
  
  # observe({
  #   req(input$patientData, input$AIMVariables) # Ensure both files are uploaded
  #   allDataValue <- all_data_filtered() # Get the current value of all_data
  #   variableNames <- variables()$variable # Assuming variables() returns a dataframe with a column 'variable'
  #   
  #   combinedData <- NULL
  #   
  #   for (var in variableNames) {
  #     
  #     bcxsubtracted_data<-allDataValue %>% 
  #       dplyr::select(DonorID:Stim,{{var}}) %>% 
  #       tidyr::pivot_wider(names_from = Stim, values_from = {{var}}) %>% 
  #       dplyr::mutate(Covid_WT_sum=ibc(bc(COVID_WT)-bc(DMSO)), Covid_BA4_5_sum=ibc(bc(COVID_BA4_5)-bc(DMSO))) %>% 
  #       tidyr::pivot_longer(DMSO:Covid_BA4_5_sum, names_to="Stimulant", values_to = {{var}}) %>%
  #       dplyr::mutate(across(DonorID:Stimulant,as.factor)) %>%
  #       dplyr::filter(Stimulant != "COVID_WT" & Stimulant != "COVID_BA4_5") %>%
  #       dplyr::mutate(Stimulant=fct_recode(Stimulant, COVID_WT="Covid_WT_sum", COVID_BA4_5="Covid_BA4_5_sum")) %>% 
  #       dplyr::mutate(Timepoint=fct_relevel(Timepoint,"VY","V2","V3"))
  #     
  #     # Remove the COVID_BA4_5.x and COVID_WT.x columns
  #     bcxsubtracted_data <- bcxsubtracted_data %>% 
  #       dplyr::select(DonorID, Stimulant, Timepoint, {{var}})
  #     
  #     # First iteration: initialize combinedData with bcxsubtracted_data. Otherwise, merge the new data into combinedData by your unique identifiers
  #     if (is.null(combinedData)) {
  #       combinedData <- bcxsubtracted_data
  #     } else {
  #       combinedData <- combinedData %>%
  #         dplyr::left_join(bcxsubtracted_data, by = c("DonorID", "Stimulant", "Timepoint"))
  #     }
  #   }
  #   write.csv(combinedData, "PREVENT_Boxcox_Combined.csv")
  #   lastCreatedFile("PREVENT_Boxcox_Combined.csv")
  # })
  
  
  
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


