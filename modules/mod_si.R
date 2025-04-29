## ── modules/mod_si.R ────────────────────────────────────────────────────
siUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "Stimulation-Index",
    sidebarLayout(
      sidebarPanel(
        numericInput(ns("lambda"),     "Lambda (L)",        0.5, step = .1),
        textInput(  ns("unstim"),      "Unstimulated label", "unstimulated"),
        textAreaInput(ns("stims"),     "Stim labels (comma separated)",
                      "CMV_peptides, Infanrix, CytoStim"),
        fileInput(  ns("file"),        "Patient data (CSV)"),
        selectInput(ns("vars"),  "AIM variables", choices = NULL, multiple = TRUE),
        selectInput(ns("group"), "Grouping columns", choices = NULL, multiple = TRUE),
        selectInput(ns("stimcol"), "Stimulant column", choices = NULL),
        downloadButton(ns("dl"), "Download SI data")
      ),
      mainPanel(
        DTOutput(ns("preview"))
      )
    )
  )
}

siServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    ## ------------ load data ------------------------------------------------
    rawDat <- reactive({
      req(input$file)
      tryCatch(
        read.csv(input$file$datapath, stringsAsFactors = FALSE),
        error = function(e) { shinyalert("Read error", e$message, "error"); NULL }
      )
    })
    
    ## ---- update selectInput choices once data available -------------------
    observe({
      req(rawDat())
      cols <- names(rawDat())
      updateSelectInput(session, "vars",     choices = cols)
      updateSelectInput(session, "group",    choices = cols)
      updateSelectInput(session, "stimcol",  choices = cols)
    })
    
    ## ---- preview raw data --------------------------------------------------
    output$preview <- renderDT({ req(rawDat()); datatable(rawDat()) })
    
    ## ---- grouped & ordered copy -------------------------------------------
    grouped <- reactive({
      req(rawDat(), input$group)
      rawDat() |> arrange(across(all_of(input$group)))
    })
    
    ## ---- SI transformation -------------------------------------------------
    siData <- reactive({
      req(grouped(), input$vars, input$stimcol)
      unstim_label <- input$unstim
      stimlabels   <- strsplit(input$stims, ",\\s*")[[1]]
      L            <- input$lambda
      
      # Work per grouping block
      grouped() |>
        group_by(across(all_of(input$group))) |>
        group_modify(function(df, key) {
          
          # cache unstim row(s)
          unstim <- df |> filter(.data[[input$stimcol]] == unstim_label)
          if (nrow(unstim) == 0) return(df)   # skip if no control row
          
          # iterate through each stim row and each variable
          df |> mutate(across(
            all_of(input$vars),
            ~ {
              stim_val <- .
              control  <- unstim[[cur_column()]]
              SI(x1 = stim_val, x0 = control, L = L)
            }
          ))
        }) |>
        ungroup()
    })
    
    ## ---- download handler --------------------------------------------------
    output$dl <- downloadHandler(
      filename = function() paste(Sys.Date(), "SI_data.csv", sep = "_"),
      content  = function(f) {
        write.csv(siData(), f, row.names = FALSE)
      }
    )
  })
}
