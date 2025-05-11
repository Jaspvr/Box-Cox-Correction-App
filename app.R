# Libraries and modules
source("global.R",      local = FALSE)
source("modules/mod_boxcox.R", local = TRUE)
source("modules/mod_si.R",     local = TRUE)
source("modules/mod_lambda.R",     local = TRUE)

# Define UI
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  navbarPage(
    "Box-Cox Correction App",
    boxcoxUI("boxcox"),
    siUI("si"),
    lambdaUI("lambda"),
    tabPanel("About",
             tags$br(),
             tags$p("Analysis of T cell activation-induced marker (AIM) assay data requires normalization of AIM+ cell frequencies to background AIM+ frequencies in an unstimulated control. 
Subtracting or dividing by the unstimulated control each have specific disadvantages and can amplify technical variability in the assay. 
The Box-Cox correction is an innovative method with features of both division and linear subtraction, allowing a more sophisticated correction for unstimulated AIM+ cell frequencies that better aligns with the mathematical properties of AIM datasets and reduces technical variability."),
             tags$br(),
             tags$p("To take advantage of the Box-Cox correction, upload your full AIM dataset and the set of variables to be corrected. The Box-Cox Correction App will immediately return the corrected values which are then ready for data display or statistical analysis.")
    )
  )
)

server <- function(input, output, session) {
  boxcoxServer("boxcox")       # <<— call module server
  siServer("si")
  lambdaServer("lambda")
}

# Create Shiny app object
shinyApp(ui = ui, server = server)