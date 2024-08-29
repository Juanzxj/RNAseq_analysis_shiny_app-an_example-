#functions_module

library(clusterProfiler)
library(org.Rn.eg.db)
library(ggplot2)
library(dplyr)
library(writexl)
library(stringr)
library(tibble)
library(tidyr)

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
                       legend.title = element_text(size = 15),  # Increased legend title size
                       legend.text = element_text(size = 13)) +
    ggtitle("Volcano Plot")
  
  # Highlight top significant genes
  top_genes <- res %>% arrange(padj) %>% head(num_genes)
  p <- p + geom_text_repel(data = top_genes, aes(label = Gene.name), size = 4)
  
  return(p)
}

perform_go_enrichment <- function(genes, ont) {
  ego <- enrichGO(
    gene         = genes,
    OrgDb        = org.Rn.eg.db,
    keyType      = 'SYMBOL',
    ont          = ont,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable     = TRUE
  )
  
  top_ego <- head(ego, 25)
  return(top_ego)
}

perform_kegg_enrichment <- function(genes) {
  genes_entrez <- mapIds(org.Rn.eg.db, 
                         keys = genes, 
                         keytype = "SYMBOL", 
                         column = "ENTREZID",
                         multiVals = "first")
  
  ekk <- enrichKEGG(
    gene = genes_entrez,
    organism = 'rno',
    keyType = 'kegg',
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  
  top_ekk <- head(ekk, 25)
  return(top_ekk)
}

generate_barplot <- function(enrichment_df, title) {
  enrichment_df$ID_Description <- paste(enrichment_df$ID, ": ", enrichment_df$Description)
  enrichment_df$ID_Description <- str_wrap(enrichment_df$ID_Description, width=70)
  ggplot(enrichment_df, aes(x = reorder(ID_Description, -p.adjust), y = Count, fill = -log10(p.adjust))) +
    geom_bar(stat="identity", width = 0.8) +
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
    xlab(" ") +
    ylab("Number of Enriched Genes") +
    ggtitle(title)
}
