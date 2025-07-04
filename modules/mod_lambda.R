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

plotres <- function(res, tol = 0.1) {
  
  pick_idx <- function(v, want_zero = FALSE) {
    rng <- range(v, na.rm = TRUE)
    
    ## a) for Beta & rho we want the value *near 0*
    if (want_zero) {
      if (min(abs(v)) > tol) return(NA_integer_)        # never close to zero
      idx <- which.min(abs(v))
    } else {
      ## b) for post.L & L.Lik we want the maximum – but not at an edge
      idx <- which.max(v)
      if (idx %in% c(1, length(v))) return(NA_integer_) # max at boundary
    }
    idx
  }

  # Pick the “best” λ for each curve
  best <- list(
    beta = pick_idx(res$beta,  want_zero = TRUE),
    post = pick_idx(res$post.L),
    rho  = pick_idx(res$rho,   want_zero = TRUE),
    ll   = pick_idx(res$L.Lik)
  )
  
  # centre coordinates for the fallback message in each panel
  centre_x <- mean(range(res$LL))
  centre_y <- list(
    beta = mean(range(res$beta)),
    post = mean(range(res$post.L)),
    rho  = 0,                                # we scale Rho to −0.5 ↔ 0.5 later
    ll   = mean(range(res$L.Lik))
  )

  # Helper that writes λ value just above the panel
  add_lab <- function(col, idx, digits = 2) {
    if (is.na(idx)) return(annotate("blank", 0, 0))
    annotate(
      "text",
      x      = res$LL[idx],
      y      = Inf,         # top edge
      vjust  = -0.2,
      hjust  = .5,
      colour = col,
      label  = bquote(lambda==.(round(res$LL[idx], digits))),
      size   = 4
    )
  }
  
  safe_geoms <- function(best_idx, col, y0 = NULL) {
    if (is.na(best_idx)) return(NULL)
    
    list(
      if (!is.null(y0))
        geom_hline(yintercept = y0, colour = "grey60"),
      geom_line(linewidth = 1.3, colour = col),
      geom_vline(
        xintercept = res$LL[best_idx],
        colour = col, linetype = "dashed", linewidth = 1
      )
    )
  }

  add_warn_panel <- function(txt = "Optimal λ not defined",
                             x_npc = .5, y_npc = .6,
                             txt_size = 5.5) {

    # one grob that contains the white box + text
    g <- grid::grobTree(
      grid::roundrectGrob(
        x = unit(x_npc, "npc"), y = unit(y_npc, "npc"),
        r = unit(4, "pt"),
        gp = grid::gpar(fill = "white", col = NA)
      ),
      grid::textGrob(
        label = txt,
        x = unit(x_npc, "npc"), y = unit(y_npc, "npc"),
        gp = grid::gpar(
          col = "grey20",
          fontsize = txt_size * ggplot2::.pt,
          fontface = "italic"
        )
      )
    )

    # add that grob as a single annotation layer
    ggplot2::annotation_custom(g)
  }

  base_thm <- theme_minimal(base_size = 13) +
    theme(
      panel.grid     = element_line(linewidth = .25),
      panel.spacing  = unit(2.5, "lines"),
      plot.margin    = margin(t = 20, r = 18, b = 14, l = 18),
      plot.title.position = "plot",
      plot.title          = element_text(size = rel(.9), hjust = 0.5, margin = margin(b = 16)),
    )

  # Four ggplots
  p_beta <- ggplot(res, aes(LL, beta)) +
    safe_geoms(best$beta, "#1f77b4", y0 = 0) +
    add_lab("#1f77b4", best$beta) +
    
    { if (is.na(best$beta))
      add_warn_panel()
    } +
    
    labs(
      title = "Beta (β) Coefficient Method\n(Linear Regression)",
      y     = "Beta (β) Coefficient",
      x     = "Lambda (λ)"
    ) + base_thm + coord_cartesian(clip = "off")

  p_post <- ggplot(res, aes(LL, post.L)) +
    safe_geoms(best$post, "#E69F00", y0 = 0) +
    add_lab("#E69F00", best$post) +
    { if (is.na(best$post))
      add_warn_panel()
    } +
    labs(
      title = "Posterior Probability Distribution Method",
      y     = "Posterior P(λ)",
      x     = "Lambda (λ)"
    ) + base_thm + coord_cartesian(clip = "off")

  p_rho <- ggplot(res, aes(LL, rho)) +
    safe_geoms(best$rho, "#D62728", y0 = 0) +
    add_lab("#D62728", best$rho) +
    { if (is.na(best$rho))
      add_warn_panel()
    } +
    scale_y_continuous(limits = c(-0.5, 0.5)) +
    labs(
      title = "Spearman Correlation Method",
      y     = "Spearman rho (ρ)",
      x     = "Lambda (λ)"
    ) + base_thm + coord_cartesian(clip = "off")

  p_ll <- ggplot(res, aes(LL, L.Lik)) +
    safe_geoms(best$ll, "#8E44AD", y0 = 0) +
    add_lab("#8E44AD", best$ll) +
    { if (is.na(best$ll))
      add_warn_panel()
    } +
    labs(
      title = "Likelihood Function Method",
      y     = "Likelihood L(λ)",
      x     = "Lambda (λ)"
    ) + base_thm + coord_cartesian(clip = "off")

  # Assemble the 2×2 grid (patchwork)
  (p_beta | p_post) / (p_rho | p_ll)
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
