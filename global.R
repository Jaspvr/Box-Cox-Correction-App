## ── global.R ─────────────────────────────────────────────────────────────
library(shiny)
library(shinythemes)
library(shinyalert)
library(DT)

library(tidyverse)      # ggplot2, dplyr, tidyr, readr, stringr …
library(openxlsx)
library(IDPmisc)
library(glue)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(shiny::hr)   
