# Libraries and modules
source("global.R",      local = FALSE)
source('modules/mod_about.R',  local = TRUE)
source("modules/mod_boxcox.R", local = TRUE)
source("modules/mod_si.R",     local = TRUE)
source("modules/mod_lambda.R",     local = TRUE)

# Define UI
ui <- fluidPage(
  theme = shinytheme("cerulean"),
    navbarPage(
    "Box-Cox Correction App",
    aboutUI("about"),
    lambdaUI("lambda"),
    boxcoxUI("boxcox"),
    siUI("si")
  )
)

server <- function(input, output, session) {
  aboutServer('about') 
  boxcoxServer("boxcox")
  siServer("si")
  lambdaServer("lambda")
}

# Create Shiny app object
shinyApp(ui = ui, server = server)