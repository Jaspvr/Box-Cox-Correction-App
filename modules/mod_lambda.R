library(shiny)
library(MASS)
library(tidyverse)
library(ggplot2)
library(patchwork)


# Functions from given code
bc <- function(x,L=0.5){if(L==0) return(log(x)) else return((x^L-1)/L)}
ibc <- function(x,L=0.5){if(L==0) return(exp(x)) else return((x*L+1)^(1/L))}

test.lambda.lik <- function(x1, x0, LL=NULL,f1=NULL, method="rlm", DAT=NULL, plotit=TRUE){
  # "best" lambda here will have slope = 0 on this model f(x1)-f(x0) ~ f(x0) 
  # (equivalent to slope of 1 on f(x1)~f(x0) )
  # We use the density distribution of the slope estimate as the posterior distribution
  # under a non informative prior
  # the distribution is non-central t we can approximate by normal
  # normal would put more weight on the estimated coef
  # regression is sensitive to extreme values and distribution should be considered an approximation
  # rlm is less sensitive to outliers
  # we assume independence of the samples
  # Spearman rho's is a more robust method, but we do not compute a posterior
  # here we simply return the rho value, which should be ~0 at best lambda
  
  # Note: negatives estimate for rho correspond to "over-corrected"
  # if x1 and x0 are not related, we expect negatives rho
  # if x1 and x0 are unrelated, the maximum for post.L is expected to be close to 1
  # when there is no lambda with beta=0, there is no reliable estimate of lambda 
  
  
  # This version also incorporates the likelihood with a fixed beta
  # data is corrected for the Jacobian of the transformation
  
  # f1 accounts for covariates. This changes the interpretation of the results
  
  if(is.null(LL)) LL <- seq(0,1,0.01)
  nona <- !is.na(x0)&!is.na(x1)
  x0 <- x0[nona]
  x1 <- x1[nona]
  n <- length(nona)
  
  DAT <- DAT[nona,]
  
  if(is.null(f1)){
    fb1 <- as.formula("yy~x")
    fb2 <- as.formula("y~offset(x)")
  }
  if(!is.null(f1)){
    #consider covariate from f1
    fb1 <- as.formula(paste0("yy~x+", f1) )
    fb2 <- as.formula(paste0("y~offset(x)+", f1) )
  }
  
  post.b <- c()
  rho <- c()
  beta <- c()
  L.Lik <- c()
  
  for(i in 1:length(LL)){
    JJ <- (LL[i]-1)*sum(log(x1)) # log Jacobian
    K <- (exp(JJ))^(1/n) #data rescaling factor to account for Jacobian
    
    y <- bc(x1, LL[i])
    x <- bc(x0, LL[i])
    
    #We rescale both x and y to keep the relationship
    y <- y/K
    x <- x/K
    yy <- y-x
    
    if(method=="lm"){
      Mb <- lm(fb1,data=DAT)
      M1 <- lm(fb2,data=DAT)
    }
    if(method=="rlm"){
      Mb <- rlm(fb1,data=DAT)
      M1 <- rlm(fb2,data=DAT)
    }
    beta[i] <- Mb$coef[2]
    resid1 <- resid(M1)
    rho[i] <- cor( resid1,x, method="spearman", use="complete.obs")
    #dfs <- Mb$df.residual
    dfs <- n - length(Mb$coefficients) # for rlm compatibility
    post.b[i] <- dt((0-Mb$coef[2])/summary(Mb)$coef[2,2],dfs)
    #post.b[i] <- dt((1-Mb$coef[2])/summary(Mb)$coef[2,2],dfs)
    
    Sig2 <- sum(resid1^2)/n
    tmp <- anova(M1)
    #J <- (LL[i]-1)*sum(log(x1)) # log Jacobian if not corrected data
    J <- 0 # if jacobian already included in data transformation
    L.Lik[i] <- (-n/2)*log(Sig2)+J
    #L.Lik[i] <- logLik(M1)+J   # equivalent
    
  }
  
  # under non informative prior for L, beside being 0-1 (here being LL)
  # we get the posterior P(L | x1,x0,b=0)
  # post.L <- post.b/sum(post.b)
  post.L <- post.b/sum(post.b, na.rm=TRUE)
  Lik <- exp(L.Lik-median(L.Lik))
  Lik <- Lik/sum(Lik)
  
  res <- data.frame(LL=LL,post.L=post.L, beta, rho, L.Lik, Lik)
  
  if(plotit){
    tmp<- plotres(res)
    print(tmp)
  }
  
  return(res)
}

plotres <- function(res, line_w = 1.4) {
  
  # indices of the ‘best’ λ for each diagnostic
  pick <- list(
    beta = which.min(abs(res$beta)),
    post = which.max(res$post.L),
    rho  = which.min(abs(res$rho)),
    ll   = which.max(res$L.Lik)
  )
  
  base_thm <- theme_minimal(base_size = 13) +
    theme(panel.grid  = element_line(linewidth = .25),
          plot.margin = margin(6, 30, 6, 6))          # space on the right for λ
  
  add_lab <- function(col, idx, digits = 2) {
    annotate(
      "text",
      x      = res$LL[idx],
      y      = Inf,          # very top edge of panel
      vjust  = -0.4,         # -ve moves text *above* the edge
      colour = col,
      hjust  = 0.5,
      label  = bquote(lambda==.(round(res$LL[idx], digits))),
      size   = 4.3
    )
  }
  
  make_panel <- function(y, ylab, col, idx, ylim = NULL) {
    ggplot(res, aes(LL, .data[[y]])) +
      { if (!is.null(ylim)) scale_y_continuous(limits = ylim) } +
      geom_line(linewidth = line_w, colour = col) +
      geom_vline(xintercept = res$LL[idx], colour = col,
                 linewidth = 1, linetype  = "dashed") +
      add_lab(col, idx) +
      labs(x = expression(lambda), y = ylab) +
      base_thm +
      coord_cartesian(clip = "off")                 # let λ text sit in margin
  }
  
  p_beta <- make_panel("beta",   expression(beta),   "#1f77b4", pick$beta)
  p_post <- make_panel("post.L", "Posterior",        "#E69F00", pick$post)
  p_rho  <- make_panel("rho",    expression(rho),    "#D62728", pick$rho,
                       ylim = c(-0.5, 0.5))
  p_ll   <- make_panel("L.Lik",  "Log–lik.",         "#8E44AD", pick$ll)
  
  (p_beta | p_post) / (p_rho | p_ll)                 # patchwork grid
}

# UI and Server
lambdaUI <- function(id, title = "Lambda (λ) Estimation Tools") {
  ns <- NS(id)

  tabPanel(
    title,
    sidebarPanel(
      fileInput(ns("csv"), "Upload CSV", accept = ".csv"),
      textInput(ns("stim_col"),   "Stimulant Column Name",  "Stim"),
      textInput(ns("stim_lvl"),   "Stimulant",            "SARSCoV2_Spike"),
      textInput(ns("unstim_lvl"), "Unstimulated Parameter",          "DMSO"),
      numericInput(ns("eps"), "Epsilon", value = 0.001, step = 0.0005),
      uiOutput(ns("aim_selector")),
      
      # Restore default placeholders and clear all placeholders
      tags$hr(),
      actionButton(ns("clear"),  "Clear inputs"),
      actionButton(ns("reset"),  "Restore defaults"),
      tags$hr(),
      
      actionButton(ns("run"), "Run")
    ),
    mainPanel(
      DT::DTOutput(ns("preview_table"), height = "250px"),   # new
      tags$hr(),
      plotOutput(ns("lambda_plot"), height = 600)
    )
  )
}


lambdaServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    defaults <- list(
      stim_col  = "Stim",
      stim_lvl  = "SARSCoV2_Spike",
      unstim_lvl = "DMSO",
      eps       = 0.001
    )
    
    # Read the csv
    data_raw <- reactive({
      req(input$csv)
      tryCatch(
        read.csv(input$csv$datapath,
                 header = TRUE, sep = ",",
                 fileEncoding = "UTF-8-BOM",
                 check.names = FALSE,
                 stringsAsFactors = FALSE),
        error = function(e) {
          shinyalert("Error",
                     paste("Error reading CSV:", e$message),
                     type = "error")
          NULL
        }
      )
    })
    
    # Preview table
    output$preview_table <- renderDT({
      req(data_raw())
      datatable(
        data_raw(),
        options = list(pageLength = 5, scrollX = TRUE)
      )
    })
    
    # Populate AIM selector
    output$aim_selector <- renderUI({
      req(data_raw())
      selectInput(
        session$ns("aim_var"),
        "AIM Variable (column)",
        choices  = names(data_raw()),
        selected = names(data_raw())[1]
      )
    })
    
    # Clear inputs
    observeEvent(input$clear, {
      updateTextInput( session, "stim_col",   value = "")
      updateTextInput( session, "stim_lvl",   value = "")
      updateTextInput( session, "unstim_lvl", value = "")
      updateNumericInput(session, "eps",      value = NA)
    })
    
    # Reset inputs to defaults
    observeEvent(input$reset, {
      updateTextInput( session, "stim_col",   value = defaults$stim_col)
      updateTextInput( session, "stim_lvl",   value = defaults$stim_lvl)
      updateTextInput( session, "unstim_lvl", value = defaults$unstim_lvl)
      updateNumericInput(session, "eps",      value = defaults$eps)
    })

    # populate AIM variable selector once the CSV is in
    output$aim_selector <- renderUI({
      req(input$csv)
      df <- read.csv(input$csv$datapath, nrows = 1, stringsAsFactors = FALSE)
      selectInput(session$ns("aim_var"), "AIM Variable (column)",
                  choices = names(df), selected = names(df)[1])
    })

    # run on "Run"
    observeEvent(input$run, {
      req(data_raw(), input$aim_var, input$stim_col,
          nzchar(input$stim_lvl), nzchar(input$unstim_lvl))

      # User input
      epsilon <- input$eps
      data    <- data_raw()

      Stim   <- data %>% filter(.data[[input$stim_col]] == input$stim_lvl)
      Unstim <- data %>% filter(.data[[input$stim_col]] == input$unstim_lvl)

      Stim   <- Stim[[input$aim_var]]   + epsilon
      Unstim <- Unstim[[input$aim_var]] + epsilon

      # open a fresh graphics device inside renderPlot
      output$lambda_plot <- renderPlot({
        test.lambda.lik(x1 = Stim, x0 = Unstim)  # this ultimately prints the plots
      })
    })
  })
}
