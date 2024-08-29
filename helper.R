# This file will contain utility functions that are used across different modules in your Shiny app, such as the functions for 
# generating plots or performing enrichment analysis.

library(clusterProfiler)
library(org.Rn.eg.db)
library(ggplot2)
library(writexl)
library(stringr)
library(tibble)
library(tidyr)
library(dplyr)
library(ggrepel)

# Helper function to generate a volcano plot
plot_volcano <- function(res, num_genes) {
  res <- as.data.frame(res) %>% na.omit()
  res$color <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) >= 1.0,
                      ifelse(res$log2FoldChange >= 1.0, 'Up', 'Down'), 'Not')
  color <- c(Up = "#FD9145", Down = "#03AFD8", Not = "#2F343F")
  p <- ggplot(res, aes(log2FoldChange, -log10(padj), color = color)) +
    geom_point(size = 4, alpha = 0.5) +
    scale_color_manual(values = color) +
    labs(x = "log2 (fold change)", y = "-log10 (q-value)") +
    geom_hline(yintercept = -log10(0.05), lty = 4, col = "grey", lwd = 0.6) +
    geom_vline(xintercept = c(-1.0, 1.0), lty = 4, col = "grey", lwd = 0.6) +
    theme_bw() + theme(axis.title = element_text(size = 16),
                       axis.text = element_text(size = 14),
                       plot.title = element_text(size = 15, hjust = 0.5),
                       legend.title = element_text(size = 15),
                       legend.text = element_text(size = 13)) +
    ggtitle("Volcano Plot")
  
  # Highlight top significant genes
  top_genes <- res %>% arrange(padj) %>% head(num_genes)
  p <- p + geom_text_repel(data = top_genes, aes(label = Gene.name), size = 4)
  
  return(p)
}


# Helper function to perform GO enrichment analysis
perform_go_enrichment <- function(genes, ont) {
  # Perform GO enrichment analysis
  ego <- enrichGO(
    gene         = genes,
    OrgDb        = org.Rn.eg.db,  # This should match the organism you're working with
    keyType      = 'SYMBOL',      # Gene identifier type, e.g., SYMBOL, ENTREZID
    ont          = ont,           # Ontology: BP (Biological Process), MF (Molecular Function), or CC (Cellular Component)
    pAdjustMethod = "BH",         # Method for adjusting p-values, e.g., Benjamini-Hochberg
    qvalueCutoff = 0.05,          # Q-value cutoff for significant terms
    readable     = TRUE           # Convert gene identifiers to more readable names
  )
  
  # Return the top 25 GO terms
  top_ego <- head(ego, 10)
  return(top_ego)
}


# Helper function to perform KEGG enrichment analysis
perform_kegg_enrichment <- function(genes) {
  # Convert gene symbols to Entrez IDs
  genes_entrez <- mapIds(
    org.Rn.eg.db,      # Annotation database for Rattus norvegicus
    keys = genes,      # Gene symbols to convert
    keytype = "SYMBOL",# Input key type (gene symbols)
    column = "ENTREZID", # Output key type (Entrez IDs)
    multiVals = "first"  # Handle multiple mappings (take the first match)
  )
  
  # Perform KEGG enrichment analysis
  ekk <- enrichKEGG(
    gene = genes_entrez,  # List of Entrez IDs
    organism = 'rno',     # KEGG organism code for Rattus norvegicus
    keyType = 'kegg',     # Key type for KEGG (typically 'kegg' for pathway IDs)
    pAdjustMethod = "BH", # Method for adjusting p-values (Benjamini-Hochberg)
    qvalueCutoff = 0.05   # Q-value cutoff for significant terms
  )
  
  # Return the top 25 KEGG pathways
  top_ekk <- head(ekk, 10)
  return(top_ekk)
}


# Helper function to generate a bar plot
generate_barplot <- function(enrichment_df, title) {
  # Combine ID and Description for readability
  enrichment_df$ID_Description <- paste(enrichment_df$ID, ": ", enrichment_df$Description)
  
  # Wrap text for better display on the plot
  enrichment_df$ID_Description <- stringr::str_wrap(enrichment_df$ID_Description, width = 70)
  
  # Generate the bar plot
  ggplot(enrichment_df, aes(x = reorder(ID_Description, -p.adjust), y = Count, fill = -log10(p.adjust))) +
    geom_bar(stat = "identity", width = 0.8) +
    coord_flip() +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      axis.text.x = element_text(size = 15),   # Adjusted for readability
      axis.text.y = element_text(size = 15),
      legend.text = element_text(size = 15)
    ) +
    xlab("") +
    ylab("Number of Enriched Genes") +
    ggtitle(title)
}

