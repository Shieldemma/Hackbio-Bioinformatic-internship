## Data Loading and Preprocessing

# Load necessary libraries
library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(gplots)
library(viridis)
library(org.Hs.eg.db)
library(AnnotationDbi)
# Load the data
gene_data <- read_csv("https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv") %>%
  column_to_rownames(var = "...1")

# Convert Ensembl IDs to gene symbols
ensembl_ids <- rownames(gene_data)
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
gene_symbols <- na.omit(gene_symbols)

# Update gene_data with gene symbols
gene_data <- gene_data[names(gene_symbols), ]
rownames(gene_data) <- gene_symbols

# Create sample info dataframe
sample_info <- data.frame(
  SampleID = colnames(gene_data),
  tissue_type = sapply(colnames(gene_data), function(barcode) {
    sample_type <- substr(barcode, 14, 15)
    switch(sample_type,
           "01" = "Solid Tissue Normal",
           "02" = "Primary Tumor",
           "06" = "Metastatic",
           "11" = "Solid Tissue Normal",
           "12" = "Primary Blood Derived Cancer",
           "14" = "Primary Blood Derived Cancer",
           "Unknown")
  }),
  row.names = colnames(gene_data)
)

## DESeq2 Analysis

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(gene_data),
                              colData = sample_info,
                              design = ~ tissue_type)

# Pre-filtering
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Set factor levels
dds$tissue_type <- relevel(dds$tissue_type, ref = "Solid Tissue Normal")

# Run DESeq
dds <- DESeq(dds)
res <- results(dds)

# Filter results
res_filtered <- res %>%
  as.data.frame() %>%
  filter(!is.na(padj), !is.na(log2FoldChange), padj < 0.05, abs(log2FoldChange) > 2) %>%
  arrange(desc(abs(log2FoldChange)))

# Get top 10 up and down regulated genes
top_genes <- bind_rows(
  head(res_filtered, 10) %>% mutate(regulation = "Up"),
  tail(res_filtered, 10) %>% mutate(regulation = "Down")
)

## Visualization

# Heatmap function
create_heatmap <- function(data, filename, title) {
  pdf(file = filename, width = 15, height = 12)
  heatmap.2(as.matrix(data),
            scale = "row",
            col = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                                     "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(100),
            trace = "none",
            dendrogram = "both",
            main = title,
            xlab = "Samples",
            ylab = "Genes",
            margins = c(12, 8),
            cexRow = 0.8,
            cexCol = 1,
            srtCol = 45,
            adjCol = c(1,1),
            keysize = 1.5,
            key.title = "Expression",
            key.xlab = "Z-score",
            density.info = "none",
            colsep = 1:ncol(data),
            sepcolor = "white",
            sepwidth = c(0.01, 0.01),
            lhei = c(2, 10),
            lwid = c(2, 10))
  dev.off()
}

# Create heatmaps
create_heatmap(gene_data, "heatmap_all_genes.pdf", "Gene Expression Heatmap - All Genes")
create_heatmap(gene_data[rownames(top_genes), ], "heatmap_top_genes.pdf", "Gene Expression Heatmap - Top Differentially Expressed Genes")

# Boxplot function
create_boxplot <- function(data, genes, filename, title) {
  plot_data <- data %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    filter(gene %in% genes) %>%
    pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
    left_join(sample_info, by = c("sample" = "SampleID")) %>%
    left_join(res_filtered %>% rownames_to_column("gene") %>% select(gene, log2FoldChange, padj),
              by = "gene")
  
  ggplot(plot_data, aes(x = tissue_type, y = log2(count + 1), fill = tissue_type)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ gene, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = c("Solid Tissue Normal" = "lightblue", "Primary Tumor" = "darkblue")) +
    labs(y = "log2CPM", x = NULL, title = title) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "none") +
    geom_text(data = plot_data %>% group_by(gene) %>% slice(1),
              aes(x = tissue_type, y = Inf, label = sprintf("q = %.2e", padj)),
              vjust = 1.5, size = 3, inherit.aes = FALSE)
  
  ggsave(filename, width = 15, height = 10, units = "in")
}

# Create boxplots
normalized_counts <- counts(dds, normalized = TRUE)
create_boxplot(normalized_counts[rownames(top_genes), ], rownames(top_genes), "boxplot_top_genes.pdf", "Expression of Top Differentially Expressed Genes")


#==============================================================================
#----------- 04# Enrichment Analaysis for significant genes -------------------
#------------------------------------------------------------------------------
#-=============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# Assuming 'res' is your DESeq2 results object
significant_genes <- res_filtered
gene_list <- rownames(significant_genes)

# Convert gene symbols to Entrez IDs
gene_symbols <- rownames(significant_genes)
gene_list <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Filter out any NA values that may result from conversion
gene_list <- gene_list[!is.na(gene_list$ENTREZID), "ENTREZID"]

# Perform GO enrichment analysis
Go.BP <- enrichGO(gene          = gene_list,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "BP", # Biological Process
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2)

go_data <- as.data.frame(Go.BP)
go_data$Count <- as.numeric(as.character(go_data$Count))
go_data$p.adjust <- as.numeric(as.character(go_data$p.adjust))

barplot_BP <- ggplot(go_data, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red", high = "blue", 
                      trans = "log10",  # Use log scale for better color distribution
                      name = "Adjusted\np-value") +
  labs(x = "Count", y = NULL, title = "GO Biological Process Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(go_data$Count) * 1.1))

# Save the plot as a PNG file
ggsave("GO_Biological_Process_Enrichment.png", plot = barplot_BP, width = 10, height = 8, dpi = 300)

# ========================================== GO with MF ================================
# GO enrichment analysis
Go.MF <- enrichGO(gene          = gene_list,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "MF", # Molecular Function 
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2)
Go.MF <- as.data.frame(Go.MF)

Go.MF$Count <- as.numeric(as.character(Go.MF$Count))
Go.MF$p.adjust <- as.numeric(as.character(Go.MF$p.adjust))

barplot_MF <- ggplot(Go.MF, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red", high = "blue", 
                      trans = "log10",  # Use log scale for better color distribution
                      name = "Adjusted\np-value") +
  labs(x = "Count", y = NULL, title = "GO Molecular Function Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(go_data$Count) * 1.1))

# Save the plot as a PNG file
ggsave("GO_Molecular_Function_Enrichment.png", plot = barplot_MF, width = 10, height = 8, dpi = 300)


#======================================= GO with Cellular Componenet ============================
# GO enrichment analysis
Go.CC <-enrichGO(gene          = gene_list,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = "ENTREZID",
                 ont           = "CC", # Molecular Function 
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2)
Go.CC<- as.data.frame(Go.CC@result)

Go.CC$Count <- as.numeric(as.character(Go.CC$Count))
Go.CC$p.adjust <- as.numeric(as.character(Go.CC$p.adjust))

barplot_CC <- ggplot(Go.CC, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red", high = "blue", 
                      trans = "log10",  # Use log scale for better color distribution
                      name = "Adjusted\np-value") +
  labs(x = "Count", y = NULL, title = "GO Cellular Component Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(go_data$Count) * 1.1))

# Save the plot as a PNG file
ggsave("GO_Cellular_Componenet_Enrichment.png", plot = barplot_CC, width = 10, height = 8, dpi = 300)
# ============================================================================================
# KEGG pathway enrichment analysis
# ============================================================================================
ekegg <- enrichKEGG(gene         = gene_list,
                    organism     = 'hsa', # Human
                    pvalueCutoff = 0.05)

if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
  stop("The GO enrichment analysis did not return any results.")
}

if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
  stop("The KEGG enrichment analysis did not return any results.")
}


ekegg@result
# Bar plot for GO enrichment
kegg_plot <- barplot(ekegg, showCategory = 5)
barplot()

ggsave("KEGG_Enrichment.png", plot = kegg_plot, width = 16, height = 20, dpi = 300)
# Dot plot for KEGG enrichment

dp <- dotplot(ekegg, showCategory = 5)
# save dot plot 
ggsave("kEGG_Enrichment_dotplot.png", plot = dp, width = 10, height = 8, dpi = 300)


# Load necessary libraries
library(ggplot2)
library(forcats)

# Assuming 'ekegg' is your enrichKEGG result object
# Convert the enrichKEGG result to a data frame
ekegg_df <- as.data.frame(ekegg@result)

# Assuming 'Count' and 'GeneRatio' are columns in your data frame
ekegg_df$RichFactor <- ekegg_df$Count / as.numeric(sapply(strsplit(ekegg_df$GeneRatio, "/"), `[`, 2))

# Create a lollipop plot using the calculated RichFactor
lp <- ggplot(ekegg_df[1:5, ], aes(x = RichFactor, y = fct_reorder(Description, RichFactor))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_gradientn(colours = c("#f7ca64", "#46bac2", "#7e62a3"), trans = "log10") +
  scale_size_continuous(range = c(2, 10)) +
  theme_minimal() +
  xlab("Rich Factor") +
  ylab(NULL) +
  ggtitle("KEGG Enrichment Lollipop Plot")

ggsave("kEGG_Enrichment_Lollipop.png", plot = lp, width = 10, height = 8, dpi = 300)
