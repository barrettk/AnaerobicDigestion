
suppressPackageStartupMessages({
  library(rsconnect)
  library(shiny)
  library(shinyAce,quietly = T)
  library(shinydashboard)
  #library(shiny.semantic)
  library(shinyWidgets,quietly = T)
  library(rmarkdown)
  library(markdown)
  library(knitr)
  library(webshot)
  library(kableExtra)
  library(polynom,quietly = T)
  library(plyr,quietly = T)
  library(dplyr,quietly = T)
  library(ggplot2,quietly = T)
  library(tidyverse,quietly = T)
  library(sn,quietly = T)
  library(rlang)
  library(wrapr)
  library(gridExtra,quietly = T)
  library(extrafont,quietly = T)
  library(reshape2,quietly = T)
  library(mrgsolve)
  
  
})

options(mrgsolve_mread_quiet=TRUE)


# tell shiny to log all reactivity
# library(reactlog)
# options(shiny.reactlog = TRUE)

# Compile the model
mod <- mrgsolve::mread("AnaerobicDigestionShinyV3.cpp")
.model <- paste(mrgsolve:::code(mod), collapse='\n')
# assign("mod",mod,envir = globalenv())
assign(".model",.model,envir = globalenv())

source("Bacteria_kinetics.R",local=globalenv())

SensChoices <- c("Acidogen Conc.","Acetogen Conc.","Methanogen Conc.","Bacteroid Conc.","Bacteria Decay Rate","Head Space Ratio","Temperature","Number of Wells","WT % of Guar Gum","WT % of PEG-400","WT % of Methanol","WT % of Isopropanol")
