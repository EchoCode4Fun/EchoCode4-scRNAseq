# Add metadata to indicate the condition
control_neutrophils$condition <- "control"
neutrophil_cluster$condition <- "Mutifc"

# Merge the control and HNifc neutrophils into one combined object
neutrophils_combined <- merge(control_neutrophils, y = neutrophil_cluster, add.cell.ids = c("control", "Mutifc"))

# Normalize the data
neutrophils_combined <- NormalizeData(neutrophils_combined)

# Find variable features
neutrophils_combined <- FindVariableFeatures(neutrophils_combined)

# Scale the data
neutrophils_combined <- ScaleData(neutrophils_combined)

# Perform PCA
neutrophils_combined <- RunPCA(neutrophils_combined)

# Perform UMAP for visualization
neutrophils_combined <- RunUMAP(neutrophils_combined, dims = 1:20)

# Find clusters
neutrophils_combined <- FindNeighbors(neutrophils_combined, dims = 1:20)
neutrophils_combined <- FindClusters(neutrophils_combined)

# Join data layers
neutrophils_combined <- JoinLayers(neutrophils_combined)

# UMAP plot colored by condition
DimPlot(neutrophils_combined, reduction = "umap", group.by = "condition")

# Perform DE analysis on neutrophils
de_genes_neutrophils <- FindMarkers(neutrophils_combined, ident.1 = "control", ident.2 = "Mutifc", group.by = "condition")

# View the top differentially expressed genes
head(de_genes_neutrophils)

## now for visualization 
# Create a volcano plot
EnhancedVolcano(de_genes_neutrophils,
                lab = rownames(de_genes_neutrophils),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Differential Expression of Neutrophils',
                pCutoff = 0.05,
                FCcutoff = 0.25,
                pointSize = 3.0,
                labSize = 3.5)

# Select top 20 DE genes based on adjusted p-value
top_genes <- head(rownames(de_genes_neutrophils[order(de_genes_neutrophils$p_val_adj), ]), 20)

# Extract normalized expression data for the top DE genes
expr_data <- FetchData(neutrophils_combined, vars = top_genes)

# Add condition metadata
expr_data$condition <- neutrophils_combined$condition

# Create a heatmap
pheatmap(t(expr_data[, -ncol(expr_data)]), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         annotation_col = data.frame(condition = expr_data$condition))

###Perform GOE analysis on neutrophils 


# Extract DE genes
de_genes <- rownames(de_genes_neutrophils[de_genes_neutrophils$p_val_adj < 0.05, ])

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(de_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Perform GO enrichment analysis
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP", # Biological Process
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

# View GO results
head(go_results)

# Create the dot plot
dotplot(go_results, showCategory = 10) +
  ggtitle("Top 10 GO Terms Enriched: H445Yifc vs Control Neutrophils")


# Function to plot GO enrichment results for a specific cluster
plot_go_enrichment <- function(go_results, cluster) {
  if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
    dotplot(go_results, showCategory = 10, x = "GeneRatio") +
      scale_color_gradient(low = "blue", high = "red", name = "p.adjust") +
      ggtitle(paste("GGO Enrichment Analysis for DE Genes: H445Yifc vs Control HSCs", cluster)) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  }
}






#### NOW compare HSCs

# Add metadata to indicate the condition
control_stem_cells$condition <- "control"
stem_cell_cluster$condition <- "Mutifc"

# Merge the control and HNifc neutrophils into one combined object
HSC_combined <- merge(control_stem_cells, y = stem_cell_cluster, add.cell.ids = c("control", "Mutifc"))

# Normalize the data
HSC_combined <- NormalizeData(HSC_combined)

# Find variable features
HSC_combined <- FindVariableFeatures(HSC_combined)

# Scale the data
HSC_combined <- ScaleData(HSC_combined)

# Perform PCA
HSC_combined <- RunPCA(HSC_combined)

# Perform UMAP for visualization
HSC_combined <- RunUMAP(HSC_combined, dims = 1:20)

# Find clusters
HSC_combined <- FindNeighbors(HSC_combined, dims = 1:20)
HSC_combined <- FindClusters(HSC_combined)

# Join data layers
HSC_combined <- JoinLayers(HSC_combined)

# UMAP plot colored by condition
DimPlot(HSC_combined, reduction = "umap", group.by = "condition")

# Perform DE analysis on neutrophils
de_genes_HSC <- FindMarkers(HSC_combined, ident.1 = "control", ident.2 = "Mutifc", group.by = "condition")

# View the top differentially expressed genes
head(de_genes_HSC)

## now for visualization 
# Create a volcano plot
EnhancedVolcano(de_genes_HSC,
                lab = rownames(de_genes_HSC),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Differential Expression of H445Yifc vs Ctrl HSC',
                pCutoff = 0.05,
                FCcutoff = 0.25,
                pointSize = 1.0,
                labSize = 3.5)


###Perform GOE analysis on HSCs

# Extract DE genes
de_genes <- rownames(de_genes_HSC[de_genes_HSC$p_val_adj < 0.05, ])

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(de_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Perform GO enrichment analysis
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP", # Biological Process
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

# View GO results
head(go_results)

# Dotplot for GO enrichment results
dotplot(go_results, showCategory = 8) +
  ggtitle("GO Enrichment Analysis for DE Genes: H445Yifc vs Control HSCs")
