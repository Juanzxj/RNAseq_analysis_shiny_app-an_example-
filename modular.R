# This file will define the UI and server modules for each feature in your app. Each module will 
# encapsulate a specific part of your application, such as the Volcano Plot or GO Enrichment.


# Volcano Plot UI Module
volcano_plot_UI <- function(id) {
  ns <- NS(id)  # Namespace function for the module
  tagList(
    fluidRow(
      box(
        title = "Volcano Plot", 
        status = "primary", 
        solidHeader = TRUE, 
        width = 12,
        numericInput(ns("log2fc_threshold_volcano"), "Log2 Fold Change Threshold", value = 1, min = 0),
        numericInput(ns("padj_threshold_volcano"), "Adjusted p-value Threshold", value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput(ns("num_significant_genes"), "Number of Top Significant Genes", value = 10, min = 1),
        actionButton(ns("update_volcano"), "Update Volcano Plot"),
        plotOutput(ns("volcano_plot")),
        numericInput(ns("volcano_plot_width"), "Plot Width (in inches)", value = 10, min = 1),
        numericInput(ns("volcano_plot_height"), "Plot Height (in inches)", value = 8, min = 1),
        selectInput(ns("volcano_plot_type"), "File Type", choices = c("png", "pdf", "svg")),
        downloadButton(ns("output_volcano"), "Output Volcano Plot")
      )
    )
  )
}


# Volcano Plot Server Module
volcano_plot_Server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive expression to filter data based on input thresholds
    filtered_data <- reactive({
      req(data())
      df <- data()
      df <- df %>%
        dplyr::filter((log2FoldChange >= input$log2fc_threshold_volcano | 
                         log2FoldChange <= -input$log2fc_threshold_volcano) & 
                        padj <= input$padj_threshold_volcano)
      return(df)
    })
    
    # Render the Volcano plot
    output$volcano_plot <- renderPlot({
      input$update_volcano
      isolate({
        req(filtered_data())
        plot_volcano(filtered_data(), input$num_significant_genes)
      })
    })
    
    # Download handler for the Volcano plot
    output$output_volcano <- downloadHandler(
      filename = function() {
        paste("volcano_plot", ".", input$volcano_plot_type, sep = "")
      },
      content = function(file) {
        ggplot2::ggsave(file, 
                        plot = plot_volcano(filtered_data(), input$num_significant_genes), 
                        device = input$volcano_plot_type, 
                        width = input$volcano_plot_width, 
                        height = input$volcano_plot_height)
      }
    )
  })
}


# GO Enrichment UI Module
go_enrichment_UI <- function(id) {
  ns <- NS(id)  # Namespace function for the module
  tagList(
    fluidRow(
      box(
        title = "GO Enrichment", 
        status = "primary", 
        solidHeader = TRUE, 
        width = 12,
        numericInput(ns("log2fc_threshold"), "Log2 Fold Change Threshold", value = 1, min = 0),
        numericInput(ns("padj_threshold"), "Adjusted p-value Threshold", value = 0.05, min = 0, max = 1, step = 0.01),
        selectInput(ns("go_ontology"), "Select Ontology", choices = c("BP", "MF", "CC")),
        actionButton(ns("run_go_enrichment"), "Run GO Enrichment"),
        plotOutput(ns("go_plot"), height = "600px"),
        numericInput(ns("go_plot_width"), "Plot Width (in inches)", value = 10, min = 1),
        numericInput(ns("go_plot_height"), "Plot Height (in inches)", value = 8, min = 1),
        selectInput(ns("go_plot_type"), "File Type", choices = c("png", "pdf", "svg")),
        downloadButton(ns("output_go_plot"), "Output GO Plot"),
        downloadButton(ns("download_go_terms"), "Download GO Terms")
      )
    )
  )
}


# GO Enrichment Server Module
go_enrichment_Server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive expression to perform GO enrichment analysis
    go_result <- eventReactive(input$run_go_enrichment, {
      req(data())
      filtered_data <- data()
      if (!"Gene.name" %in% colnames(filtered_data)) {
        stop("Error: 'Gene.name' column not found in the data")
      }
      DEGs <- filtered_data$Gene.name
      perform_go_enrichment(DEGs, input$go_ontology) %>% 
        head(5)
    })
    
    # Render the GO enrichment bar plot
    output$go_plot <- renderPlot({
      req(go_result())
      generate_barplot(as.data.frame(go_result()), "Top 5 GO Terms")
    })
    
    # Download handler for the GO enrichment plot
    output$output_go_plot <- downloadHandler(
      filename = function() {
        paste("go_enrichment", ".", input$go_plot_type, sep = "")
      },
      content = function(file) {
        ggplot2::ggsave(file, 
                        plot = generate_barplot(as.data.frame(go_result()), "Top 5 GO Terms"), 
                        device = input$go_plot_type, 
                        width = input$go_plot_width, 
                        height = input$go_plot_height)
      }
    )
    
    # Download handler for the GO terms as an Excel file
    output$download_go_terms <- downloadHandler(
      filename = function() { "go_terms.xlsx" },
      content = function(file) {
        req(data())
        filtered_data <- data()
        DEGs <- filtered_data$Gene.name
        go_result <- perform_go_enrichment(DEGs, input$go_ontology) %>% 
          head(5)
        writexl::write_xlsx(as.data.frame(go_result), file)
      }
    )
  })
}


# KEGG Enrichment UI Module
kegg_enrichment_UI <- function(id) {
  ns <- NS(id)  # Namespace function for the module
  tagList(
    fluidRow(
      box(
        title = "KEGG Enrichment", 
        status = "primary", 
        solidHeader = TRUE, 
        width = 12,
        numericInput(ns("log2fc_threshold_kegg"), "Log2 Fold Change Threshold", value = 1, min = 0),
        numericInput(ns("padj_threshold_kegg"), "Adjusted p-value Threshold", value = 0.05, min = 0, max = 1, step = 0.01),
        actionButton(ns("run_kegg_enrichment"), "Run KEGG Enrichment"),
        plotOutput(ns("kegg_plot"), height = "600px"),
        numericInput(ns("kegg_plot_width"), "Plot Width (in inches)", value = 10, min = 1),
        numericInput(ns("kegg_plot_height"), "Plot Height (in inches)", value = 8, min = 1),
        selectInput(ns("kegg_plot_type"), "File Type", choices = c("png", "pdf", "svg")),
        downloadButton(ns("output_kegg_plot"), "Output KEGG Plot"),
        downloadButton(ns("download_kegg_terms"), "Download KEGG Terms")
      )
    )
  )
}


# KEGG Enrichment Server Module
kegg_enrichment_Server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive expression to perform KEGG enrichment analysis
    kegg_result <- eventReactive(input$run_kegg_enrichment, {
      req(data())
      filtered_data <- data()
      if (!"Gene.name" %in% colnames(filtered_data)) {
        stop("Error: 'Gene.name' column not found in the data")
      }
      DEGs <- filtered_data$Gene.name
      perform_kegg_enrichment(DEGs) %>% 
        head(10)
    })
    
    # Render the KEGG enrichment bar plot
    output$kegg_plot <- renderPlot({
      req(kegg_result())
      generate_barplot(as.data.frame(kegg_result()), "Top 10 KEGG Enrichment Results")
    })
    
    # Download handler for the KEGG enrichment plot
    output$output_kegg_plot <- downloadHandler(
      filename = function() {
        paste("kegg_enrichment", ".", input$kegg_plot_type, sep = "")
      },
      content = function(file) {
        ggplot2::ggsave(file, 
                        plot = generate_barplot(as.data.frame(kegg_result()), "Top 10 KEGG Enrichment Results"), 
                        device = input$kegg_plot_type, 
                        width = input$kegg_plot_width, 
                        height = input$kegg_plot_height)
      }
    )
    
    # Download handler for the KEGG terms as an Excel file
    output$download_kegg_terms <- downloadHandler(
      filename = function() { "kegg_terms.xlsx" },
      content = function(file) {
        req(data())
        filtered_data <- data()
        DEGs <- filtered_data$Gene.name
        kegg_result <- perform_kegg_enrichment(DEGs) %>% 
          head(10)
        writexl::write_xlsx(as.data.frame(kegg_result), file)
      }
    )
  })
}


# Data Upload UI Module
data_upload_UI <- function(id) {
  ns <- NS(id)  # Namespace function for the module
  tagList(
    fluidRow(
      box(
        title = "Data Status", 
        status = "primary", 
        solidHeader = TRUE, 
        width = 12,
        uiOutput(ns("upload_status"))  # Use the ns() function for the UI output ID
      )
    )
  )
}


# Data Upload Server Module
data_upload_Server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive expression to load the data
    data <- reactive({
      readRDS("BB_Y1 vs M_DEGs_Unique.rds")
    })
    
    # Render the upload status message
    output$upload_status <- renderUI({
      if (is.null(data())) {
        h4("No data loaded")
      } else {
        h4("Data successfully loaded")
      }
    })
    
    # Return the data reactive expression so it can be used elsewhere
    return(data)
  })
}


