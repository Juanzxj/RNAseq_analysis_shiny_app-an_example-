# This file will serve as the main entry point of your Shiny application. It will define the overall UI layout, 
#load the necessary modules, and call the server functions.

library(shiny)
library(shinydashboard)

# Source the modular and helper functions
source("modular.R")
source("helper.R")

ui <- dashboardPage(
  dashboardHeader(title = "RNASeq Analysis App"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Upload", tabName = "data_upload", icon = icon("upload")),
      menuItem("Volcano Plot", tabName = "volcano_plot", icon = icon("fire")),
      menuItem("GO Enrichment", tabName = "go_enrichment", icon = icon("list")),
      menuItem("KEGG Enrichment", tabName = "kegg_enrichment", icon = icon("list"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "data_upload", data_upload_UI("data_upload")),
      tabItem(tabName = "volcano_plot", volcano_plot_UI("volcano_plot")),
      tabItem(tabName = "go_enrichment", go_enrichment_UI("go_enrichment")),
      tabItem(tabName = "kegg_enrichment", kegg_enrichment_UI("kegg_enrichment"))
    )
  )
)

# Main server function
server <- function(input, output, session) {
  
  # Call the data upload module and get the data
  data <- data_upload_Server("data_upload")
  
  # Pass the data to other modules
  volcano_plot_Server("volcano_plot", data)
  go_enrichment_Server("go_enrichment", data)
  kegg_enrichment_Server("kegg_enrichment", data)
}



shinyApp(ui, server)

