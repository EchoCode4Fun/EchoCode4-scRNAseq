# Load required libraries
library(celldex)
library(Seurat)
library(Matrix)

# Load the ImmGen reference dataset
ref <- celldex::ImmGenData()
ref_genes <- rownames(ref)
if (length(ref_genes) == 0) {
  stop("Reference dataset does not contain gene names.")
}

# Define the path to your files
data_dir <- "scRNAseq2/"

# Read the matrix file
HNifc_matrix <- readMM(file.path(data_dir, "GSM5365327_HN878_matrix.mtx.gz"))

# Read the barcodes and features
HNifc_barcodes <- readLines(file.path(data_dir, "GSM5365327_HN878_barcodes.tsv.gz"))
HNifc_features <- readLines(file.path(data_dir, "GSM5365327_HN878_features.tsv"))

# Split the character vector by tabs
split_features <- strsplit(HNifc_features, "\t")

# Convert the list of split elements into a data frame
features_df <- do.call(rbind, split_features)
features_df <- as.data.frame(features_df, stringsAsFactors = FALSE)

# Assign column names for clarity (adjust the names as needed)
colnames(features_df) <- c("ID", "GeneSymbol", "Description")

# Check the structure of the resulting data frame
str(features_df)

# Extract gene symbols from the second column
HNifc_gene_symbols <- features_df$GeneSymbol

# Print the gene symbols to verify
print(HNifc_gene_symbols)

# Convert barcodes and gene symbols to proper format
colnames(HNifc_matrix) <- HNifc_barcodes
rownames(HNifc_matrix) <- HNifc_gene_symbols

# Make row names unique by appending a suffix
rownames(HNifc_matrix) <- make.unique(rownames(HNifc_matrix))

# Verify that all row names are now unique
print(any(duplicated(rownames(HNifc_matrix))))

# Create a Seurat object
HNifc_seurat_object <- CreateSeuratObject(counts = HNifc_matrix, project = "scRNAseq", min.cells = 3, min.features = 200)

##Quality check
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
HNifc_seurat_object[["percent.mt"]] <- PercentageFeatureSet(HNifc_seurat_object, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(HNifc_seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(HNifc_seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HNifc_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
HNifc_seurat_object <- subset(HNifc_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 5)

# Normalize the data
HNifc_seurat_object <- NormalizeData(HNifc_seurat_object)

# Find variable features (default nfeature = 2000)
HNifc_seurat_object <- FindVariableFeatures(HNifc_seurat_object, selection.method = "vst", nfeatures = 1000)
VariableFeaturePlot(HNifc_seurat_object)
#OR
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(HNifc_seurat_object), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(HNifc_seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale the data
HNifc_seurat_object <- ScaleData(HNifc_seurat_object, features = rownames(HNifc_seurat_object))


# Perform PCA
HNifc_seurat_object <- RunPCA(HNifc_seurat_object, features = VariableFeatures(object = HNifc_seurat_object))

#Visualiza and examine PCA
VizDimLoadings(HNifc_seurat_object, dims = 1:2, reduction = "pca")
DimPlot(HNifc_seurat_object, reduction = "pca") + NoLegend()
DimHeatmap(HNifc_seurat_object, dims = 1, cells = 500, balanced = TRUE) #can change the dim numebr
DimHeatmap(HNifc_seurat_object, dims = 1:15, cells = 500, balanced = TRUE) #this shows all dims


# **ELBOW PLOT FOR DIMENSION DETERMINATION**
# Calculate the percentage of variance explained by each PC
pca_variance <- HNifc_seurat_object@reductions$pca@stdev / sum(HNifc_seurat_object@reductions$pca@stdev) * 100

# Create the elbow plot
plot(pca_variance, type = "b", xlab = "Principal Component", ylab = "Percentage of Variance Explained", main = "Elbow Plot for PCA")

# (Optionally) Add a line for visual guidance
lines(pca_variance, type = "b")

ElbowPlot(HNifc_seurat_object)



# Find neighbors
HNifc_seurat_object <- FindNeighbors(HNifc_seurat_object, dims = 1:20)

# Cluster the cells
HNifc_seurat_object <- FindClusters(HNifc_seurat_object, resolution = 0.6)

# Run UMAP for visualization
library(ggplot2)
HNifc_seurat_object <- RunUMAP(HNifc_seurat_object, dims = 1:20)
DimPlot(HNifc_seurat_object, reduction = "umap", label = TRUE) + NoLegend() + ggtitle("HNifc Bone Marrow_WOAnnotation")

##   II: Assess the quality of your clustering --- run Silhoutte analysis
# 1. Extract Cluster Labels and UMAP Embeddings
cluster_labels <- Idents(HNifc_seurat_object)  # Get cluster assignments
umap_embeddings <- HNifc_seurat_object@reductions$umap@cell.embeddings  # Extract UMAP coordinates

# 2. Calculate Silhouette Scores
library(cluster)  # Load the 'cluster' package for silhouette() function
dist_matrix <- dist(umap_embeddings)  # Calculate distance matrix from UMAP embeddings
sil <- silhouette(as.numeric(cluster_labels), dist_matrix) # Perform silhouette analysis

# 3. Visualize and Interpret Silhouette Plot
plot(sil, main = "Silhouette Plot for Seurat Clustering") 

# 4. (Optional) Calculate Average Silhouette Score
avg_sil_score <- mean(sil[, "sil_width"])
print(paste("Average Silhouette Score:", avg_sil_score))

## If you just want the table, run these codes below
# Calculate average silhouette width for each cluster
cluster_silhouette <- aggregate(sil[, "sil_width"], by = list(Cluster = cluster_labels), FUN = mean)
# Rename the columns for clarity
colnames(cluster_silhouette) <- c("Cluster", "Average Silhouette Width")
# Print the table
print(cluster_silhouette)



# 1. Find all markers for each cluster
HNifc_markers <- FindAllMarkers(HNifc_seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Run ROC test to identity those subtle markers
# 2. Iterate through each cluster, get ROC results, and find overlaps 
for (cluster in levels(HNifc_seurat_object)) {
  
  # Get markers from FindAllMarkers for the current cluster
  cluster_markers <- HNifc_markers %>% 
    filter(cluster == !!cluster) 
  
  # Apply ROC test using FindMarkers
  roc_results <- FindMarkers(HNifc_seurat_object, 
                             ident.1 = cluster, 
                             features = cluster_markers$gene, 
                             test.use = "roc", 
                             only.pos = TRUE)
  
  # Find overlapping genes (genes significant in both tests)
  overlapping_genes <- intersect(cluster_markers$gene, 
                                 roc_results$gene[roc_results$p_val_adj < 0.05]) # Adjust p-value threshold as needed
  
  # Add "Overlap" column to HNifc_markers (TRUE if gene overlaps, FALSE otherwise)
  HNifc_markers <- HNifc_markers %>%
    mutate(Overlap = ifelse(gene %in% overlapping_genes & cluster == !!cluster, TRUE, FALSE))
}

##Visulization of the markers 
#Vinplot: shows expression probability 
VlnPlot(HNifc_seurat_object, features = c("Cd34", "Flt3","Ly6a","Ly6e","Slamf1","Slamf2","Kit","Flk2")) #to identify HSCs
VlnPlot(HNifc_seurat_object, features = c("Cd34", "Flt3","Cd48","Ly6a/e","Kit","Slamf1")) #assess LT-HSC (lack of Lin,cd34, cd135/Flt3, cd48; enriched for Sca-1+ and c-Kit+ and cd150/slamf1)
#assess ST-HSC (lack of lin and Flt3)
#assess MPP (lack of lin, but enriched for Cd34, Flt3, and Sca-1 and c-kit)
VlnPlot(HNifc_seurat_object, features = c("Ly6a","Ly6e","Kit")) 

VlnPlot(HNifc_seurat_object, features = c("Ptprc","Prom1","Slamf1","Cd48","Kit")) #assess HSCs subtypes 
VlnPlot(HNifc_seurat_object, features = c("Cd34", "Cd105","Endoglin","Itga2b")) #to identify Megakaryocyte-Erythroid progenitors
VlnPlot(HNifc_seurat_object, features = c("Cd34", "Cd16/32","Cd32","Csf3r")) #to identify Granulocyte-Monocyte progenitors
VlnPlot(HNifc_seurat_object, features = c("Cd34", "Flt3","Fcgr3")) #common myeloid progenitors
VlnPlot(HNifc_seurat_object, features = c("Cd34", "Flt3","Flk2","Cd117")) #multipotent progenitors
VlnPlot(HNifc_seurat_object, features = c("Flt3","Il7r")) #Common Lymphoid progenitors
VlnPlot(HNifc_seurat_object, features = c("S100a8", "S100a9","Elane","Mpo","Ly6g","Cxcr2","Cd11b")) #to identify neutrophils
VlnPlot(HNifc_seurat_object, features = c("S100a8", "S100a9","Elane","Mpo","Ly6g","Adgre1","Cxcr2","Cd11b","c-fms","Ccr2","Cd115"))
VlnPlot(HNifc_seurat_object, features = c("Adgre1", "Csf1r","Cd68","Itgam","Mertk","H2-Ab1","H2-Ea")) #to identify monocytes macrophages
VlnPlot(HNifc_seurat_object, features = c("Cd19", "Ms4a1","Cd79a")) #identify B cells
VlnPlot(HNifc_seurat_object, features = c("Cd4", "Cd8","Cd3e")) #identify T cells
VlnPlot(HNifc_seurat_object, features = c("Klrb1c", "Ncr1")) #identify NK cells
VlnPlot(HNifc_seurat_object, features = c("Cd11c", "H2-Ab1")) #identify DC cells
VlnPlot(HNifc_seurat_object, features = c("Siglec-F", "Ccr3")) #identify Eosinophils cells
VlnPlot(HNifc_seurat_object, features = c("Fcer1a", "Cpa3","Fcer1b")) #identify Basophils cells
VlnPlot(HNifc_seurat_object, features = c("Kit", "Cpa3","FcεRI")) #identify Mast cells
VlnPlot(HNifc_seurat_object, features = c("Pf4", "Itga2b","Vwf")) #identify Megakaryocytes
VlnPlot(HNifc_seurat_object, features = c("Hba")) #identify Erythrocytes 
VlnPlot(HNifc_seurat_object, features = c("Nt5e","Thy1","Eng","Runx2","Alpl","Col1a1")) #identify stromal cells
VlnPlot(HNifc_seurat_object, features = c("Pecam1","Cdh5","Vwf","Lepr","Kitl")) #endothelial cells 
VlnPlot(HNifc_seurat_object, features = c("Lin","Ly6a","Ckit"))
VlnPlot(HNifc_seurat_object, features = c("Cd47","Thbs1"))

#Vinplot: plot raw counts 
VlnPlot(HNifc_seurat_object, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

#Featureplot: shows global distribution of a series of markers
FeaturePlot(HNifc_seurat_object, features = c("Ms4a1", "Gnly", "Cd3e", "Cd14", "Fcer1a", "Fcgr3a", "Lyz", "Ppbp", "Cd8a"))

#Heatmap: show the top 10 markers (a good way to annotate clusters)
HNifc_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(HNifc_seurat_object, features = top10$gene) + NoLegend() + theme(axis.text.y = element_text(size = 6))


# Extract the top 10 markers for each cluster and save to Excel
library(dplyr)
HNifc_top10_markers <- HNifc_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) %>% 
  arrange(cluster, dplyr::desc(avg_log2FC))
print(HNifc_top10_markers)

library(writexl)
write_xlsx(HNifc_top10_markers, "scRNAseq#2_HNifc_top10_markers3.xlsx")


##Manually add annotation
# Combine and annotate clusters
cluster_annotations <- c(
  "0" = "Hematopoietic Stem/Progenitor Cells",
  "1" = "Neutrophils",
  "2" = "Neutrophils",
  "3" = "Neutrophils",
  "4" = "Macrophages",
  "5" = "Hematopoietic Stem/Progenitor Cells",
  "6" = "Megakaryocytes",
  "7" = "Neutrophils",
  "8" = "pro-B cells",
  "9" = "Basophils",
  "10" = "Erythroblasts",
  "11" = "Neutrophils"
)

##Manually add annotation
# Combine and annotate clusters
cluster_annotations <- c(
  "0" = "Dendritic Cells",
  "1" = "Megakaryocytes",
  "2" = "Phase I Neutrophils",
  "3" = "M1 Macrophages",
  "4" = "Neutrophils",
  "5" = "Macrophages",
  "6" = "Hematopoietic Stem/Progenitor Cells(erythroblasts)",
  "7" = "Hematopoietic Stem/Progenitor Cells",
  "8" = "Mast Cells",
  "9" = "Hematopoietic Stem/Progenitor Cells(erythroblasts)"
)

##Manually add annotation (for CellChat)
# Combine and annotate clusters
cluster_annotations <- c(
  "0" = "Dendritic Cells",
  "1" = "Megakaryocytes",
  "2" = "Neutrophils",
  "3" = "Macrophages",
  "4" = "Neutrophils",
  "5" = "Macrophages",
  "6" = "Hematopoietic Stem/Progenitor Cells",
  "7" = "Hematopoietic Stem/Progenitor Cells",
  "8" = "Mast Cells",
  "9" = "Hematopoietic Stem/Progenitor Cells"
)

# Apply annotations to the Seurat object
HNifc_seurat_object <- RenameIdents(HNifc_seurat_object, cluster_annotations)

# Visualize the annotated clusters
DimPlot(HNifc_seurat_object, reduction = "umap", label = TRUE, label.size = 5) + # Adjust label.size here
  NoLegend() + 
  ggtitle("HNifc Whole Bone Marrow") + 
  theme(text = element_text(size = 10, face = "bold"))



###Automatic addtion of annotation
library(SingleR)
# Extract gene symbols from test genes
test_genes <- rownames(HNifc_seurat_object)
ref_genes <- toupper(ref_genes)

# Find common genes
common_genes <- intersect(test_genes, ref_genes)
if (length(common_genes) == 0) {
  stop("No common genes found between test and reference datasets.")
} else {
  print(paste(length(common_genes), "common genes found."))
}

# Subset the Seurat object and reference data to the common genes
counts <- GetAssayData(HNifc_seurat_object, slot = "counts")[rownames(HNifc_seurat_object) %in% common_genes, ]
ref <- ref[rownames(ref) %in% common_genes, , drop = FALSE]
common_genes <- intersect(rownames(counts), rownames(ref))
counts <- counts[common_genes, ]
ref <- ref[common_genes, , drop = FALSE]

# Extract cluster labels
clusters <- Idents(HNifc_seurat_object)

# Perform SingleR annotation
singleR_results <- SingleR(test = counts, ref = ref, labels = ref$label.main, clusters = clusters)

# Create a named vector of SingleR annotations
cluster_annotations <- singleR_results$labels
names(cluster_annotations) <- levels(Idents(HNifc_seurat_object))

# Rename the identities in the Seurat object using the named vector
HNifc_seurat_object <- RenameIdents(HNifc_seurat_object, cluster_annotations)

# Plot UMAP with annotations
DimPlot(HNifc_seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() + ggtitle("HNifc Bone Marrow_Annotated")

# Save the Seurat object
saveRDS(HNifc_seurat_object, file = "annotated_HNifc_seurat_object.rds")

# Export cluster annotations
write.csv(Idents(seurat_object), file = "cluster_annotations.csv")

### NOW STUDY THE NEUTROPHILS SUBCLUSTER

# Subset the Neutrophil cluster
HNifc_neutrophils <- subset(HNifc_seurat_object, idents = "Neutrophils")

# Normalize the subset data
HNifc_neutrophils <- NormalizeData(HNifc_neutrophils)

# Find variable features in the subset
HNifc_neutrophils <- FindVariableFeatures(HNifc_neutrophils, selection.method = "vst", nfeatures = 2000)

# Scale the subset data
HNifc_neutrophils <- ScaleData(HNifc_neutrophils, features = rownames(HNifc_neutrophils))

# Perform PCA on the subset
HNifc_neutrophils <- RunPCA(HNifc_neutrophils, features = VariableFeatures(object = HNifc_neutrophils))

# Find neighbors in the subset
HNifc_neutrophils <- FindNeighbors(HNifc_neutrophils, dims = 1:10)

# Cluster the subset
HNifc_neutrophils <- FindClusters(HNifc_neutrophils, resolution = 0.5)

# Run UMAP for the subset
HNifc_neutrophils <- RunUMAP(HNifc_neutrophils, dims = 1:10)

# Plot UMAP for neutrophil subclusters
DimPlot(HNifc_neutrophils, reduction = "umap", label = TRUE, label.size = 5) + 
  ggtitle("post-HNifc Neutrophil Subclusters")

# Find markers for each neutrophil subcluster
HNifc_neutrophil_markers <- FindAllMarkers(HNifc_neutrophils, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Extract the top 10 markers for each subcluster
top10_HNifc_neutrophil_markers <- HNifc_neutrophil_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


# Write the top 10 markers for neutrophil subclusters to a CSV file
write.csv(top10_HNifc_neutrophil_markers, "scRNAseq2/top10_post_HNifc_neutrophil_markers.csv", row.names = FALSE)

# Draw a heatmap of the top markers for neutrophil subclusters
DoHeatmap(HNifc_neutrophils, features = top10_HNifc_neutrophil_markers$gene) + 
  ggtitle("Heatmap of Top Markers for post-HNifc Neutrophil Subclusters")

## To do GO enrichment analysis 

# Load necessary libraries
library(clusterProfiler)
library(dplyr)
library(ggplot2)

# Merge with top10_control_neutrophil_markers to ensure we have cluster information
top10_HNifc_neutrophil_markers <- merge(top10_HNifc_neutrophil_markers, gene_entrez_ids, by.x = "gene", by.y = "SYMBOL")

# Perform GO enrichment analysis for each cluster
go_results_list <- list()

for (cluster in unique(top10_HNifc_neutrophil_markers$cluster)) {
  cluster_genes <- top10_HNifc_neutrophil_markers %>%
    filter(cluster == !!cluster) %>%
    pull(ENTREZID)
  
  go_results <- enrichGO(
    gene = cluster_genes,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",  # BP: Biological Process, CC: Cellular Component, MF: Molecular Function
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  go_results_list[[as.character(cluster)]] <- go_results
}

# Function to plot GO enrichment results for a specific cluster
plot_go_enrichment <- function(go_results, cluster) {
  if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
    dotplot(go_results, showCategory = 10, x = "GeneRatio") +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(low = "blue", high = "red", name = "p.adjust") +
      scale_size(range = c(3, 8), name = "Gene Count") +
      ggtitle(paste("GO Enrichment for Subcluster", cluster)) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  }
}

# Plot GO enrichment results for each cluster
for (cluster in names(go_results_list)) {
  go_results <- go_results_list[[cluster]]
  print(plot_go_enrichment(go_results, cluster))
}

##### Move on to subset HSC

# Subset Stem cells
HNifc_stem_cells <- subset(HNifc_seurat_object, idents = "Hematopoietic Stem/Progenitor Cells")

# Normalize the data
HNifc_stem_cells <- NormalizeData(HNifc_stem_cells)

# Identify the most variable features
HNifc_stem_cells <- FindVariableFeatures(HNifc_stem_cells)

# Scale the data
HNifc_stem_cells <- ScaleData(HNifc_stem_cells)

# Perform PCA
HNifc_stem_cells <- RunPCA(HNifc_stem_cells)

# Perform UMAP dimensional reduction
HNifc_stem_cells <- RunUMAP(HNifc_stem_cells, dims = 1:10)

# Identify subclusters within Stem cells
HNifc_stem_cells <- FindClusters(HNifc_stem_cells, resolution = 0.5)  # Adjust resolution as needed

# Find markers for each subcluster
stem_cell_markers <- FindAllMarkers(HNifc_stem_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Extract the top 10 markers for each subcluster
top10_stem_cell_markers <- stem_cell_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Draw UMAP plot
DimPlot(HNifc_stem_cells, reduction = "umap", label = TRUE, pt.size = 0.5) + ggtitle("UMAP of post-HNifc HSCs")

# Draw Heatmap
top10_genes <- unique(top10_stem_cell_markers$gene)
DoHeatmap(HNifc_stem_cells, features = top10_genes) + ggtitle("Heatmap of Top 10 Markers per HNifc HSC Subcluster")


# Load necessary libraries
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(org.Mm.eg.db)
library(enrichplot)

# Convert gene symbols to Entrez IDs
gene_symbols <- top10_stem_cell_markers$gene
gene_entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Merge with top10_stem_cell_markers to ensure we have cluster information
top10_stem_cell_markers <- merge(top10_stem_cell_markers, gene_entrez_ids, by.x = "gene", by.y = "SYMBOL", all.x = TRUE)

# Rename ENTREZID column to avoid confusion
colnames(top10_stem_cell_markers)[which(colnames(top10_stem_cell_markers) == "ENTREZID.y")] <- "ENTREZID"

# Remove rows with NA values in the ENTREZID column
top10_stem_cell_markers <- top10_stem_cell_markers %>% filter(!is.na(ENTREZID))

# Perform GO enrichment analysis for each subcluster
go_results_list <- list()

for (cluster in unique(top10_stem_cell_markers$cluster)) {
  cluster_genes <- top10_stem_cell_markers %>% filter(cluster == !!cluster) %>% pull(ENTREZID)
  
  go_results <- enrichGO(
    gene = cluster_genes,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",  # BP: Biological Process, CC: Cellular Component, MF: Molecular Function
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  go_results_list[[as.character(cluster)]] <- go_results
}

# Function to plot GO enrichment results for a specific cluster
plot_go_enrichment <- function(go_results, cluster) {
  if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
    dotplot(go_results, showCategory = 10, x = "GeneRatio") +
      scale_color_gradient(low = "blue", high = "red", name = "p.adjust") +
      ggtitle(paste("GO Enrichment for Subcluster", cluster)) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  }
}

# Visualize GO enrichment results for each subcluster
for (cluster in names(go_results_list)) {
  go_results <- go_results_list[[cluster]]
  print(plot_go_enrichment(go_results, cluster))
}



#####Run CellChat

###Analyze cell-cell interaction 

# Load necessary libraries
library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(CellChat)
library(dplyr)
library(celldex)  # For loading reference datasets

# Assuming you have a Seurat object loaded
HNifc_seurat_object <- readRDS("annotated_HNifc_seurat_object.rds")

# Convert Seurat object to SingleCellExperiment
sce <- as.SingleCellExperiment(HNifc_seurat_object)

# Load ImmGen reference data for SingleR
immgen_ref <- celldex::ImmGenData()

# Recompute SingleR annotations using ImmGen reference data
singleR_results <- SingleR(test = sce, ref = immgen_ref, labels = immgen_ref$label.main)

# Extract the SingleR annotations and ensure they have cell names
singleR_annotations <- singleR_results$labels
names(singleR_annotations) <- colnames(sce)

# Extract counts and metadata from Seurat object
counts <- GetAssayData(HNifc_seurat_object, slot = "counts")
meta <- HNifc_seurat_object@meta.data

# Align SingleR annotations with Seurat metadata
common_cells <- intersect(names(singleR_annotations), rownames(meta))
singleR_annotations <- singleR_annotations[common_cells]
meta <- meta[common_cells, ]

# Add SingleR annotations to the metadata
meta$singleR_annotations <- singleR_annotations

# Update the Seurat object with the new metadata
HNifc_seurat_object@meta.data <- meta

# Verify the updated metadata
head(meta)

# Create a CellChat object
HNifc_cellchat <- createCellChat(object = counts, meta = meta, group.by = "singleR_annotations")

# Set the database to use for CellChat (e.g., for mouse)
CellChatDB <- CellChatDB.mouse  # or CellChatDB.human for human data

# Use the selected database
HNifc_cellchat@DB <- CellChatDB

# Preprocess the data
HNifc_cellchat <- subsetData(HNifc_cellchat)  # Subset the expression data of signaling genes
HNifc_cellchat <- identifyOverExpressedGenes(HNifc_cellchat)
HNifc_cellchat <- identifyOverExpressedInteractions(HNifc_cellchat)

# Compute the communication probability
HNifc_cellchat <- computeCommunProb(HNifc_cellchat)

# Filter out the interactions with low probability
HNifc_cellchat <- filterCommunication(HNifc_cellchat, min.cells = 10)

# Infer cell-cell communication network
HNifc_cellchat <- computeCommunProbPathway(HNifc_cellchat)

# Aggregate the cell-cell communication network
HNifc_cellchat <- aggregateNet(HNifc_cellchat)

# Check the cell-cell communication network
print(HNifc_cellchat@net)

# List available pathways
pathways.available <- names(HNifc_cellchat@netP$pathway)
print(pathways.available)



#Work by using CD52
netVisual_aggregate(cellchat, signaling = "CD52", layout = "chord")
netVisual_heatmap(cellchat, signaling = "CD52", color.heatmap = "Reds")

###Signalling pathway networks (3 visualization methods)

##1. Circle Plot
pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

##2. Chord diagram
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

##3. Heatmap
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

######What about ligand-pair?
netAnalysis_contribution(cellchat, signaling = pathways.show)

#then, visualize it. Here we have Cxcl2-Cxcr2
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")


#######We can also visualize the aggregated cell-cell communication network. For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 0.5)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", vertex.label.cex = 0.5)

#Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Alternatively, I want to focus on Stem cells, Neutrophils, Monocytes, and Macrophages
# Define group size
groupSize <- as.numeric(table(cellchat@idents))

# Define the weight matrix
mat <- cellchat@net$weight

# Specify the cell groups of interest
cell_groups <- c("Neutrophils", "Stem cells", "Monocytes", "Macrophages")

# Loop through each specified cell group
for (cell_group in cell_groups) {
  # Find the row index for the current cell group
  i <- which(rownames(mat) == cell_group)
  
  if (length(i) > 0) { # Ensure the cell group exists in the matrix
    # Initialize a new matrix for the current cell group
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    
    # Save the plot to a PNG file with larger dimensions
    png(filename = paste0(cell_group, "_plot.png"), width = 1200, height = 900, res = 150) # Adjust dimensions and resolution as needed
    
    # Plot the signaling network for the current cell group
    netVisual_circle(
      mat2, 
      vertex.weight = groupSize, 
      weight.scale = TRUE, 
      edge.weight.max = max(mat), 
      title.name = rownames(mat)[i],
      vertex.label.cex = 0.5  # Decrease the font size
    )
    
    # Close the PNG device
    dev.off()
  } else {
    message(paste("Cell group", cell_group, "not found in the matrix."))
  }
}


#####Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
##1.Bubble Plot
table(cellchat@idents) 

## Between Neutrophils and Stem cells 
# Define the indices for neutrophils and stem cells
source_index <- 9  # Neutrophils
target_index <- 12  # Stem cells

# Visualize the interactions using netVisual_bubble
png(filename = "neutrophils_to_stem_cells_interactions.png", width = 1200, height = 900, res = 150) # Adjust dimensions and resolution as needed
netVisual_bubble(cellchat, sources.use = source_index, targets.use = target_index, remove.isolate = FALSE)
dev.off()

# Define the indices for neutrophils and stem cells
source_index <- 12  # Stem cells
target_index <- 9  # Neutrophils

# Visualize the interactions using netVisual_bubble
png(filename = "stem_cells_to_neutrophils_interactions.png", width = 1200, height = 900, res = 150) # Adjust dimensions and resolution as needed
netVisual_bubble(cellchat, sources.use = source_index, targets.use = target_index, remove.isolate = FALSE)
dev.off()

##Between Stem cells and other cells 
# Define the index for stem cells
source_index <- 12  # Stem cells

# Define the indices for all other cell groups (excluding stem cells)
# Assuming there are 16 groups in total as per your provided data
target_indices <- setdiff(1:14, source_index)

# Visualize the interactions using netVisual_bubble
png(filename = "stem_cells_to_other_cells_interactions.png", width = 1200, height = 900, res = 150) # Adjust dimensions and resolution as needed
netVisual_bubble(cellchat, sources.use = source_index, targets.use = target_indices, remove.isolate = FALSE)
dev.off()

###Show the reverse 
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
#> Comparing communications on a single object

##For violin plot 
plotGeneExpression(cellchat, signaling = "CXCL")
plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)

###Compute and visualize the network centrality scores
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

###Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Load the necessary library
library(gridExtra)

#### Visualize the dominant senders (sources) and receivers (targets) in a 2D space

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

###Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2



save.image(file = "scRNAseq_2_HNifc.RData", version = NULL, ascii = FALSE, compress = TRUE, safe = TRUE)






