library(shiny)
library(shinydashboard)
library(ggplot2)
library(readxl)
library(dplyr)
library(ggrepel)

source("functions_module.R")

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
      tabItem(tabName = "data_upload",
              fluidRow(
                box(title = "Data Status", status = "primary", solidHeader = TRUE, width = 12,
                    uiOutput("upload_status")
                )
              )
      ),
      tabItem(tabName = "volcano_plot",
              fluidRow(
                box(title = "Volcano Plot", status = "primary", solidHeader = TRUE, width = 12,
                    numericInput("log2fc_threshold_volcano", "Log2 Fold Change Threshold", value = 1, min = 0),
                    numericInput("padj_threshold_volcano", "Adjusted p-value Threshold", value = 0.05, min = 0, max = 1, step = 0.01),
                    numericInput("num_significant_genes", "Number of Top Significant Genes", value = 10, min = 1),
                    actionButton("update_volcano", "Update Volcano Plot"),
                    plotOutput("volcano_plot"),
                    numericInput("volcano_plot_width", "Plot Width (in inches)", value = 10, min = 1),
                    numericInput("volcano_plot_height", "Plot Height (in inches)", value = 8, min = 1),
                    selectInput("volcano_plot_type", "File Type", choices = c("png", "pdf", "svg")),
                    downloadButton("output_volcano", "Output Volcano Plot")
                )
              )
      ),
      tabItem(tabName = "go_enrichment",
              fluidRow(
                box(title = "GO Enrichment", status = "primary", solidHeader = TRUE, width = 12,
                    numericInput("log2fc_threshold", "Log2 Fold Change Threshold", value = 1, min = 0),
                    numericInput("padj_threshold", "Adjusted p-value Threshold", value = 0.05, min = 0, max = 1, step = 0.01),
                    selectInput("go_ontology", "Select Ontology", choices = c("BP", "MF", "CC")),
                    actionButton("run_go_enrichment", "Run GO Enrichment"),
                    plotOutput("go_plot", height = "600px"),
                    numericInput("go_plot_width", "Plot Width (in inches)", value = 10, min = 1),
                    numericInput("go_plot_height", "Plot Height (in inches)", value = 8, min = 1),
                    selectInput("go_plot_type", "File Type", choices = c("png", "pdf", "svg")),
                    downloadButton("output_go_plot", "Output GO Plot"),
                    downloadButton("download_go_terms", "Download GO Terms")
                )
              )
      ),
      tabItem(tabName = "kegg_enrichment",
              fluidRow(
                box(title = "KEGG Enrichment", status = "primary", solidHeader = TRUE, width = 12,
                    numericInput("log2fc_threshold_kegg", "Log2 Fold Change Threshold", value = 1, min = 0),
                    numericInput("padj_threshold_kegg", "Adjusted p-value Threshold", value = 0.05, min = 0, max = 1, step = 0.01),
                    actionButton("run_kegg_enrichment", "Run KEGG Enrichment"),
                    plotOutput("kegg_plot", height = "600px"),
                    numericInput("kegg_plot_width", "Plot Width (in inches)", value = 10, min = 1),
                    numericInput("kegg_plot_height", "Plot Height (in inches)", value = 8, min = 1),
                    selectInput("kegg_plot_type", "File Type", choices = c("png", "pdf", "svg")),
                    downloadButton("output_kegg_plot", "Output KEGG Plot"),
                    downloadButton("download_kegg_terms", "Download KEGG Terms")
                )
              )
      )
    )
  )
)

server <- function(input, output, session) {
  
  data <- reactive({
    readRDS("BB_Y1 vs M_DEGs_Unique.rds")
  })
  
  output$upload_status <- renderUI({
    if (is.null(data())) {
      h4("No data loaded")
    } else {
      h4("Data successfully loaded")
    }
  })
  
  filtered_data <- reactive({
    req(data())
    df <- data()
    df <- df %>%
      dplyr::filter((log2FoldChange >= input$log2fc_threshold_volcano | log2FoldChange <= -input$log2fc_threshold_volcano) & padj <= input$padj_threshold_volcano)
    return(df)
  })
  
  output$volcano_plot <- renderPlot({
    input$update_volcano
    isolate({
      req(filtered_data())
      plot_volcano(filtered_data(), input$num_significant_genes)
    })
  })
  
  output$output_volcano <- downloadHandler(
    filename = function() {
      paste("volcano_plot", ".", input$volcano_plot_type, sep = "")
    },
    content = function(file) {
      ggplot2::ggsave(file, plot = plot_volcano(filtered_data(), input$num_significant_genes), device = input$volcano_plot_type, width = input$volcano_plot_width, height = input$volcano_plot_height)
    }
  )
  
  go_result <- eventReactive(input$run_go_enrichment, {
    req(filtered_data())
    if (!"Gene.name" %in% colnames(filtered_data())) {
      stop("Error: 'Gene.name' column not found in the data")
    }
    DEGs <- filtered_data()$Gene.name
    perform_go_enrichment(DEGs, input$go_ontology) %>% 
      head(5)
  })
  
  output$go_plot <- renderPlot({
    req(go_result())
    generate_barplot(as.data.frame(go_result()), "Top 5 GO Terms")
  })
  
  output$output_go_plot <- downloadHandler(
    filename = function() {
      paste("go_enrichment", ".", input$go_plot_type, sep = "")
    },
    content = function(file) {
      ggplot2::ggsave(file, plot = generate_barplot(as.data.frame(go_result()), "Top 5 GO Terms"), device = input$go_plot_type, width = input$go_plot_width, height = input$go_plot_height)
    }
  )
  
  output$download_go_terms <- downloadHandler(
    filename = function() { "go_terms.xlsx" },
    content = function(file) {
      req(filtered_data())
      DEGs <- filtered_data()$Gene.name
      go_result <- perform_go_enrichment(DEGs, input$go_ontology) %>% 
        head(5)
      writexl::write_xlsx(as.data.frame(go_result), file)
    }
  )
  
  kegg_result <- eventReactive(input$run_kegg_enrichment, {
    req(filtered_data())
    if (!"Gene.name" %in% colnames(filtered_data())) {
      stop("Error: 'Gene.name' column not found in the data")
    }
    DEGs <- filtered_data()$Gene.name
    perform_kegg_enrichment(DEGs) %>% 
      head(10)
  })
  
  output$kegg_plot <- renderPlot({
    req(kegg_result())
    generate_barplot(as.data.frame(kegg_result()), "Top 10 KEGG Enrichment Results")
  })
  
  output$output_kegg_plot <- downloadHandler(
    filename = function() {
      paste("kegg_enrichment", ".", input$kegg_plot_type, sep = "")
    },
    content = function(file) {
      ggplot2::ggsave(file, plot = generate_barplot(as.data.frame(kegg_result()), "Top 10 KEGG Enrichment Results"), device = input$kegg_plot_type, width = input$kegg_plot_width, height = input$kegg_plot_height)
    }
  )
  
  output$download_kegg_terms <- downloadHandler(
    filename = function() { "kegg_terms.xlsx" },
    content = function(file) {
      req(filtered_data())
      DEGs <- filtered_data()$Gene.name
      kegg_result <- perform_kegg_enrichment(DEGs) %>% 
        head(10)
      writexl::write_xlsx(as.data.frame(kegg_result), file)
    }
  )
}

shinyApp(ui, server)
