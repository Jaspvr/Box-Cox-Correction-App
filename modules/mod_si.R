## ── modules/mod_si.R ────────────────────────────────────────────────────
siUI <- function(id, title = "Stimulation-Index") {
  ns <- NS(id)
  tabPanel(
    title,
    sidebarPanel(
      ## ------------- SI parameters --------------------------------------
      numericInput(ns("L"),    "Lambda (L)", value = 0.5, step = .1),
      textInput(   ns("H"),    "H (blank = auto)", value = ""),
      numericInput(ns("s"),    "Scale  s", value = 1),
      numericInput(ns("e"),    "Shift  e", value = 0),
      numericInput(ns("oob"),  "OOB.V",    value = 1e-3, min = 0),
      checkboxInput(ns("corr"),"Use automatic correction (F1(L))", value = TRUE),
      hr(),
      ## ------------- data options ---------------------------------------
      textInput(   ns("unstim"),  "Unstim label", value = "unstimulated"),
      fileInput(   ns("file"),    "Patient data (CSV)"),
      selectInput( ns("vars"),    "AIM variables", choices = NULL, multiple = TRUE),
      selectInput( ns("groups"),  "Grouping columns", choices = NULL, multiple = TRUE),
      selectInput( ns("stimCol"), "Stimulant column", choices = NULL),
      downloadButton(ns("dl"),    "Download SI data")
    ),
    mainPanel(
      DTOutput(ns("preview"))
    )
  )
}

siServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    ## ------------- updated SI() ----------------------------------------
    SI <- function(x1,x0,L=0.5,H=NULL,s=1,e=0,OOB.V=1E-3,corrected=TRUE){
      F1 <- function(L){
        sign0 <- function(v){ s <- sign(v); s[s==0] <- 1; s }
        L0 <- abs(L-0.5)+0.5
        f1 <- L0^((1-L0)/L0)
        f1 <- sign0(L-0.5)*f1 + (-sign0(L-0.5)+1)/2
        f1
      }
      d <- dim(x1)
      if(!is.null(d)&&any(dim(x0)!=d)) stop("x0 and x1 not same dimension")
      x1 <- as.vector(unlist((x1-e)/s))
      x0 <- as.vector(unlist((x0-e)/s))
      if(length(x0)!=length(x1)) stop("x0 and x1 not same length")
      if(!is.null(H)&&corrected) corrected <- FALSE
      if(corrected) H <- F1(L)
      if(is.null(H)) H <- 0
      res <- x1*NA
      if(L>0 & L<=1){
        oob <- which((x1^L+1-H)^(1/L) <  x0)
        ok  <- which((x1^L+1-H)^(1/L) >= x0)
        res[ok]  <- (x1[ok]^L - x0[ok]^L + 1 - H)^(1/L)
        res[oob] <- OOB.V
      }
      if(L==0) res <- x1/x0
      dim(res) <- d
      res
    }
    
    ## ------------- load CSV -------------------------------------------
    rawDat <- reactive({
      req(input$file)
      read.csv(input$file$datapath, stringsAsFactors = FALSE)
    })
    
    ## fill selectors
    observe({
      req(rawDat())
      cols <- names(rawDat())
      updateSelectInput(session, "vars",    choices = cols)
      updateSelectInput(session, "groups",  choices = cols)
      updateSelectInput(session, "stimCol", choices = cols)
    })
    
    output$preview <- renderDT({ req(rawDat()); rawDat() })
    
    ## ------------- compute SI ------------------------------------------
    siData <- reactive({
      req(rawDat(), input$vars, input$groups, input$stimCol)
      
      L   <- input$L
      s_v <- input$s
      e_v <- input$e
      oob <- input$oob
      H_v <- suppressWarnings(as.numeric(input$H))
      if(input$H=="" || is.na(H_v)) H_v <- NULL
      
      rawDat() |>
        arrange(across(all_of(input$groups))) |>
        group_by(across(all_of(input$groups))) |>
        group_modify(function(df, key){
          
          ctl <- df |> filter(.data[[input$stimCol]] == input$unstim)
          if (nrow(ctl)==0) return(df)
          
          ctl_vals <- ctl[input$vars]
          
          df |> mutate(across(
            all_of(input$vars),
            ~ {
              stim_val <- .
              ctl_val  <- ctl_vals[[cur_column()]]
              SI(stim_val, ctl_val, L=L, H=H_v, s=s_v, e=e_v,
                 OOB.V=oob, corrected=input$corr) %>%
                round(3)      # match your example’s rounding
            }
          ))
        }) |>
        ungroup()
    })
    
    ## ------------- download --------------------------------------------
    output$dl <- downloadHandler(
      filename = function()
        paste(Sys.Date(), "SI_data.csv", sep = "_"),
      content  = function(file) write.csv(siData(), file, row.names = FALSE)
    )
  })
}
