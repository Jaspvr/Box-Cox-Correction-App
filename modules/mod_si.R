siUI <- function(id, title = "Advanced Scale-Independent Modified Stimulation Index (mSI)") {
  ns <- NS(id)
  tabPanel(
    title,
    sidebarPanel(
      numericInput(ns("lambda"),     "Lambda (λ) value",            value = 0.5, step = 0.05, min = 0, max = 1),
      checkboxInput(ns("corrected"), "Apply F1(λ) correction", value = TRUE),
      numericInput(ns("theta_H"),    "Theta (θ) – leave NA for automatic", value = NA, step = 0.05),
      numericInput(ns("offset_e"),   "Epsilon",             value = 1e-3),
      numericInput(ns("oob"),        "Replacement value for undefined SI",  value = 1e-3, min = 0),
      textInput(   ns("unstimulated"), "Unstimulated Parameter", value = "DMSO"),
      textAreaInput(ns("stimulants"),  "Stimulants (comma separated)",
                    "SARSCoV2_Spike"),
      fileInput(   ns("patientData"),  "Upload CSV)"),
      selectInput( ns("AIMVariables"),     "AIM Variables",      choices = NULL, multiple = TRUE),
      selectInput( ns("grouping_columns"), "Grouping Columns",    choices = NULL, multiple = TRUE),
      selectInput( ns("stim_column"),      "Stimulant Column",    choices = NULL),
      
      # Restore default placeholders and clear all placeholders
      tags$hr(),
      actionButton(ns("clear"),  "Clear inputs"),
      actionButton(ns("reset"),  "Restore defaults"),
      tags$hr(),
      
      downloadButton(ns("download"),   "Download Transformed Data")
    ),
    mainPanel(
      DTOutput(ns("table1"))
    )
  )
}


siServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    SI <- function(x1, x0, L = 0.5, H = NULL, s = 1, e = 0, OOB.V = 1e-3, corrected = TRUE) {
      F1 <- function(L) {
        sign0 <- function(v) { s <- sign(v); s[s == 0] <- 1; s }
        L0  <- abs(L - 0.5) + 0.5
        f1  <- L0 ^ ((1 - L0) / L0)
        f1  <- sign0(L - 0.5) * f1 + (-sign0(L - 0.5) + 1) / 2
        f1
      }
      
      d <- dim(x1)
      if (!is.null(d) && any(dim(x0) != d)) stop("x0 and x1 not the same dimension")
      
      # scaling
      x1 <- as.vector(unlist((x1 - e) / s))
      x0 <- as.vector(unlist((x0 - e) / s))
      if (length(x0) != length(x1)) stop("x0 and x1 not the same length")
      
      if (!is.null(H) && corrected) { corrected <- FALSE; warning("corrected set to FALSE (because H supplied)") }
      if (corrected) H <- F1(L)
      if (is.null(H)) H <- 0
      
      res <- rep(NA_real_, length(x1))
      if (L > 0 && L <= 1) {
        OOB   <- which((x1^L + 1 - H)^(1/L) <  x0)
        N.OOB <- which((x1^L + 1 - H)^(1/L) >= x0)
        res[N.OOB] <- (x1[N.OOB]^L - x0[N.OOB]^L + 1 - H)^(1/L)
        res[OOB]   <- OOB.V
      }
      if (L == 0) res <- x1 / x0
      
      dim(res) <- d
      res
    }
    
    lastCreatedFile <- reactiveVal(NULL)
    
    defaults <- list(
      lambda       = 0.5,
      corrected    = TRUE,
      theta_H      = NA,
      offset_e     = 1e-3,
      oob          = 1e-3,
      unstimulated = "DMSO",
      stimulants   = "SARSCoV2_Spike"
    )
    
    # Clear inputs
    observeEvent(input$clear, {
      updateNumericInput(session, "lambda",     value = NA)
      updateCheckboxInput(session, "corrected", value = FALSE)
      updateNumericInput(session, "theta_H",    value = NA)
      updateNumericInput(session, "offset_e",   value = NA)
      updateNumericInput(session, "oob",        value = NA)
      
      updateTextInput( session, "unstimulated", value = "")
      updateTextAreaInput(session, "stimulants", value = "")
      
      updateSelectInput(session, "AIMVariables",      selected = character(0))
      updateSelectInput(session, "grouping_columns",  selected = character(0))
      updateSelectInput(session, "stim_column",       selected = character(0))
    })
    
    # Restore default values for inputs
    observeEvent(input$reset, {
      updateNumericInput(session, "lambda",     value = defaults$lambda)
      updateCheckboxInput(session, "corrected", value = defaults$corrected)
      updateNumericInput(session, "theta_H",    value = defaults$theta_H)
      updateNumericInput(session, "offset_e",   value = defaults$offset_e)
      updateNumericInput(session, "oob",        value = defaults$oob)
      
      updateTextInput( session, "unstimulated", value = defaults$unstimulated)
      updateTextAreaInput(session, "stimulants", value = defaults$stimulants)
    })
    
    # raw upload
    all_data_raw <- reactive({
      req(input$patientData)
      tryCatch(
        read.csv(input$patientData$datapath, stringsAsFactors = FALSE),
        error = function(e) { shinyalert("Error", paste("Error reading patient data:", e$message), type = "error"); NULL }
      )
    })
    
    # update selectors whenever the file changes
    observe({
      req(all_data_raw())
      updateSelectInput(session, "grouping_columns", choices = names(all_data_raw()))
      updateSelectInput(session, "AIMVariables",     choices = names(all_data_raw()))
      updateSelectInput(session, "stim_column",       choices = names(all_data_raw()))
    })
    
    # preview table
    output$table1 <- renderDT({ req(all_data_raw()); all_data_raw() })
    
    # arranged version
    all_data_filtered <- reactive({
      req(all_data_raw())
      gc <- input$grouping_columns
      tryCatch(
        all_data_raw() %>% arrange(across(all_of(gc))),
        error = function(e) { shinyalert("Error", paste("Error arranging data:", e$message), type = "error"); NULL }
      )
    })
    
    output$download <- downloadHandler(
      filename = function() paste(Sys.Date(), "si-transformed.csv", sep = "_"),
      content  = function(file) {
        req(input$patientData, input$AIMVariables)
        
        df            <- all_data_filtered(); if (is.null(df)) return()
        variableNames <- input$AIMVariables
        stim_column   <- input$stim_column
        
        stimulants    <- strsplit(input$stimulants, ",\\s*")[[1]]
        unstim_param  <- input$unstimulated
        grouping_cols <- input$grouping_columns
        
        if (length(grouping_cols) == 0) { shinyalert("Error", "Grouping columns not selected.", type = "error"); return() }
        
        # parameters for SI (user‑supplied)
        lambda_val <- input$lambda
        corrected <- isTRUE(input$corrected)
        oob_value <- input$oob
        h_val <- if (is.na(input$theta_H)) NULL else input$theta_H
        e_val <- input$offset_e
        
        s_val <- 1    # fixed scaling factor
        
        unstim_baseline <- 0.25
        
        # Transformation function applied per-group
        transform_si <- function(chunk, unstim) {
          unstim_row <- chunk %>% filter(!!sym(stim_column) == unstim)
          if (nrow(unstim_row) == 0) return(chunk)
          
          unstim_vals <- unstim_row[variableNames]
          current_stim <- chunk[[stim_column]]
          
          chunk <- chunk %>% mutate(
            across(all_of(variableNames),
                   ~ {
                     col_name <- cur_column()
                     mapply(function(val, stim_label) {
                       if (stim_label == unstim) return(unstim_baseline)
                       SI(
                         x1        = val,
                         x0        = unstim_vals[[col_name]],
                         L         = lambda_val,
                         H         = h_val,
                         s         = s_val,
                         e         = e_val,
                         OOB.V     = oob_value,
                         corrected = corrected
                       )
                     }, .x, current_stim)
                   })
          )
          chunk
        }
        
        # apply transformation
        transformed <- tryCatch(
          df %>%
            group_by(across(all_of(grouping_cols))) %>%
            group_modify(~ transform_si(.x, unstim_param)) %>%
            ungroup(),
          error = function(e) { shinyalert("Error", paste("Transformation error:", e$message), type = "error"); NULL }
        )
        
        tmp <- tempfile(fileext = ".csv")
        write.csv(transformed, tmp, row.names = FALSE)
        lastCreatedFile(tmp)
        file.copy(tmp, file)
      }
    )
  })
}
