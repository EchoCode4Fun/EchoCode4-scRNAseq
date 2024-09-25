library(Seurat)
library(SingleR)
library(celldex)
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
matrix <- readMM(file.path(data_dir, "GSM5365326_H445Y_matrix.mtx.gz"))

# Read the barcodes and features
barcodes <- readLines(file.path(data_dir, "GSM5365326_H445Y_barcodes.tsv.gz"))
features <- readLines(file.path(data_dir, "GSM5365326_H445Y_features.tsv"))

# Split the character vector by tabs
split_features <- strsplit(features, "\t")

# Convert the list of split elements into a data frame
features_df <- do.call(rbind, split_features)

# Assign column names for clarity (adjust the names as needed)
colnames(features_df) <- c("ID", "GeneSymbol", "Description")

# Convert the data frame to appropriate types (if necessary)
features_df <- as.data.frame(features_df, stringsAsFactors = FALSE)

# Check the structure of the resulting data frame
str(features_df)

# Extract gene symbols from the second column
gene_symbols <- features_df$GeneSymbol

# Print the gene symbols to verify
print(gene_symbols)

# Convert barcodes and gene symbols to proper format
colnames(matrix) <- barcodes
rownames(matrix) <- gene_symbols

# Make row names unique by appending a suffix
rownames(matrix) <- make.unique(rownames(matrix))

# Verify that all row names are now unique
print(any(duplicated(rownames(matrix))))

# Create a Seurat object
MUTifc_seurat_object <- CreateSeuratObject(counts = matrix, project = "scRNAseq", min.cells = 3, min.features = 200)

# Normalize the data
MUTifc_seurat_object <- NormalizeData(MUTifc_seurat_object)

# Scale the data
MUTifc_seurat_object <- ScaleData(MUTifc_seurat_object, features = rownames(MUTifc_seurat_object))

# Find variable features
MUTifc_seurat_object <- FindVariableFeatures(MUTifc_seurat_object, selection.method = "vst", nfeatures = 2000)
# Visualize the variable feature plot
VariableFeaturePlot(MUTifc_seurat_object)

# Perform PCA
MUTifc_seurat_object <- RunPCA(MUTifc_seurat_object, features = VariableFeatures(object = MUTifc_seurat_object))

# **ELBOW PLOT FOR DIMENSION DETERMINATION**
# Calculate the percentage of variance explained by each PC
pca_variance <- MUTifc_seurat_object@reductions$pca@stdev / sum(MUTifc_seurat_object@reductions$pca@stdev) * 100

# Create the elbow plot
plot(pca_variance, type = "b", xlab = "Principal Component", ylab = "Percentage of Variance Explained", main = "Elbow Plot for PCA")

# (Optionally) Add a line for visual guidance
lines(pca_variance, type = "b")


# Find neighbors
MUTifc_seurat_object <- FindNeighbors(MUTifc_seurat_object, dims = 1:8)

# Cluster the cells
MUTifc_seurat_object <- FindClusters(MUTifc_seurat_object, resolution = 0.5)

# Run UMAP for visualization
MUTifc_seurat_object <- RunUMAP(MUTifc_seurat_object, dims = 1:8)
DimPlot(MUTifc_seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

##   II: Assess the quality of your clustering --- run Silhoutte analysis
# 1. Extract Cluster Labels and UMAP Embeddings
cluster_labels <- Idents(MUTifc_seurat_object)  # Get cluster assignments
umap_embeddings <- MUTifc_seurat_object@reductions$umap@cell.embeddings  # Extract UMAP coordinates

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



# Find markers for each cluster
MUTifc_markers <- FindAllMarkers(MUTifc_seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

##Visulization of the markers 
#Vinplot: shows expression probability 
VlnPlot(MUTifc_seurat_object, features = c("Cd34", "Flt3","Ly6a","Ly6e","Slamf1","Slamf2","Kit")) #to identify HSCs
VlnPlot(MUTifc_seurat_object, features = c("Ptprc","Prom1","Slamf1","Cd48","Kit")) #assess HSCs subtypes 
VlnPlot(MUTifc_seurat_object, features = c("Cd34", "Cd105","Endoglin","Itga2b")) #to identify Megakaryocyte-Erythroid progenitors
VlnPlot(MUTifc_seurat_object, features = c("Cd34", "Cd16/32","Cd32","Csf3r")) #to identify Granulocyte-Monocyte progenitors
VlnPlot(MUTifc_seurat_object, features = c("Cd34", "Flt3","Fcgr3")) #common myeloid progenitors
VlnPlot(MUTifc_seurat_object, features = c("Cd34", "Flt3","Flk2","Cd117")) #multipotent progenitors
VlnPlot(MUTifc_seurat_object, features = c("Flt3","Il7r")) #Common Lymphoid progenitors
VlnPlot(MUTifc_seurat_object, features = c("S100a8", "S100a9","Elane","Mpo","Ly6g","Cxcr2","Cd11b")) #to identify neutrophils
VlnPlot(MUTifc_seurat_object, features = c("S100a8", "S100a9","Elane","Mpo","Ly6g","Adgre1","Cxcr2","Cd11b","c-fms","Ccr2","Cd115"))
VlnPlot(MUTifc_seurat_object, features = c("Adgre1", "Csf1r","Cd68","Itgam","Mertk","H2-Ab1","H2-Ea")) #to identify monocytes macrophages
VlnPlot(MUTifc_seurat_object, features = c("Cd19", "Ms4a1","Cd79a")) #identify B cells
VlnPlot(MUTifc_seurat_object, features = c("Cd4", "Cd8","Cd3e")) #identify T cells
VlnPlot(MUTifc_seurat_object, features = c("Klrb1c", "Ncr1")) #identify NK cells
VlnPlot(MUTifc_seurat_object, features = c("Cd11c", "H2-Ab1")) #identify DC cells
VlnPlot(MUTifc_seurat_object, features = c("Siglec-F", "Ccr3")) #identify Eosinophils cells
VlnPlot(MUTifc_seurat_object, features = c("Fcer1a", "Cpa3","Fcer1b")) #identify Basophils cells
VlnPlot(MUTifc_seurat_object, features = c("Cpa3","FcεRI")) #identify Mast cells
VlnPlot(MUTifc_seurat_object, features = c("Pf4", "Itga2b","Vwf")) #identify Megakaryocytes
VlnPlot(MUTifc_seurat_object, features = c("Hba")) #identify Erythrocytes 
VlnPlot(MUTifc_seurat_object, features = c("Nt5e","Thy1","Eng","Runx2","Alpl","Col1a1")) #identify stromal cells
VlnPlot(MUTifc_seurat_object, features = c("Pecam1","Cdh5","Vwf","Lepr","Kitl")) #endothelial cells 
VlnPlot(MUTifc_seurat_object, features = c("Lin","Ly6a","Ckit"))
VlnPlot(MUTifc_seurat_object, features = c("Cd47","Thbs1"))



# Extract the top 10 markers for each cluster and save to Excel
library(dplyr)
MUTifc_top10_markers <- MUTifc_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) %>% 
  arrange(cluster, dplyr::desc(avg_log2FC))
print(MUTifc_top10_markers)

library(writexl)
write_xlsx(MUTifc_top10_markers, "scRNAseq#2_MUTifc_top10_markers.xlsx")


##Manually add annotation
# Combine and annotate clusters
cluster_annotations <- c(
  "0" = "Hematopoietic Stem/Progenitor Cells",
  "1" = "Hematopoietic Stem/Progenitor Cells(MEP)",
  "2" = "Hematopoietic Stem/Progenitor Cells(erythroblasts)",
  "3" = "Dendritic Cells",
  "4" = "Hematopoietic Stem/Progenitor Cells",
  "5" = "Hematopoietic Stem/Progenitor Cells",
  "6" = "Macrophages",
  "7" = "Phase I Neutrophils",
  "8" = "Hematopoietic Stem/Progenitor Cells(erythroblasts)",
  "9" = "Hematopoietic Stem/Progenitor Cells",
  "10" = "Hematopoietic Stem/Progenitor Cells",
  "11" = "Neutrophils",
  "12" = "Mast Cells",
  "13" = "Natural Killer Cells"
)

##Manually add annotation (for CellChat)
# Combine and annotate clusters
cluster_annotations <- c(
  "0" = "Hematopoietic Stem/Progenitor Cells",
  "1" = "Hematopoietic Stem/Progenitor Cells",
  "2" = "Hematopoietic Stem/Progenitor Cells",
  "3" = "Dendritic Cells",
  "4" = "Hematopoietic Stem/Progenitor Cells",
  "5" = "Hematopoietic Stem/Progenitor Cells",
  "6" = "Macrophages",
  "7" = "Neutrophils",
  "8" = "Hematopoietic Stem/Progenitor Cells",
  "9" = "Hematopoietic Stem/Progenitor Cells",
  "10" = "Hematopoietic Stem/Progenitor Cells",
  "11" = "Neutrophils",
  "12" = "Mast Cells",
  "13" = "Natural Killer Cells"
)

# Apply annotations to the Seurat object
MUTifc_seurat_object <- RenameIdents(MUTifc_seurat_object, cluster_annotations)

# Visualize the annotated clusters
DimPlot(MUTifc_seurat_object, reduction = "umap", label = TRUE, label.size = 2.5) + # Adjust label.size here
  NoLegend() + 
  ggtitle("MUTifc Whole Bone Marrow") + 
  theme(text = element_text(size = 10, face = "bold"))


# Extract gene symbols from test genes
test_genes <- rownames(MUTifc_seurat_object)
ref_genes <- toupper(ref_genes)

# Find common genes
common_genes <- intersect(test_genes, ref_genes)
if (length(common_genes) == 0) {
  stop("No common genes found between test and reference datasets.")
} else {
  print(paste(length(common_genes), "common genes found."))
}

# Subset the Seurat object and reference data to the common genes
counts <- GetAssayData(MUTifc_seurat_object, slot = "counts")[rownames(MUTifc_seurat_object) %in% common_genes, ]
ref <- ref[rownames(ref) %in% common_genes, , drop = FALSE]
common_genes <- intersect(rownames(counts), rownames(ref))
counts <- counts[common_genes, ]
ref <- ref[common_genes, , drop = FALSE]

# Extract cluster labels
clusters <- Idents(MUTifc_seurat_object)

# Perform SingleR annotation
singleR_results <- SingleR(test = counts, ref = ref, labels = ref$label.main, clusters = clusters)

# Create a named vector of SingleR annotations
cluster_annotations <- singleR_results$labels
names(cluster_annotations) <- levels(Idents(MUTifc_seurat_object))

# Rename the identities in the Seurat object using the named vector
MUTifc_seurat_object <- RenameIdents(MUTifc_seurat_object, cluster_annotations)

# Plot UMAP with annotations
DimPlot(MUTifc_seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Save the Seurat object
saveRDS(MUTifc_MUTifc_seurat_object, file = "annotated_MUTifc_seurat_object.rds")

# Export cluster annotations
write.csv(Idents(MUTifc_seurat_object), file = "cluster_annotations.csv")


###Analyse the neutrophil clusters 

# Load required libraries
library(Seurat)
library(dplyr)

# Extract raw data from all layers
counts_matrix <- as.matrix(GetAssayData(neutrophil_subset, slot = "counts"))
data_matrix <- as.matrix(GetAssayData(neutrophil_subset, slot = "data"))
scale_data_matrix <- as.matrix(GetAssayData(neutrophil_subset, slot = "scale.data"))

# Ensure unique gene names
rownames(counts_matrix) <- make.unique(rownames(counts_matrix))
rownames(data_matrix) <- make.unique(rownames(data_matrix))
rownames(scale_data_matrix) <- make.unique(rownames(scale_data_matrix))

# Check dimensions
cat("Counts matrix dimensions: ", dim(counts_matrix), "\n")
cat("Data matrix dimensions: ", dim(data_matrix), "\n")
cat("Scale data matrix dimensions: ", dim(scale_data_matrix), "\n")

# Find common gene names
common_genes <- Reduce(intersect, list(rownames(counts_matrix), rownames(data_matrix), rownames(scale_data_matrix)))

# Subset matrices to the common genes
counts_matrix <- counts_matrix[common_genes,]
data_matrix <- data_matrix[common_genes,]
scale_data_matrix <- scale_data_matrix[common_genes,]

# Check dimensions again to confirm alignment
cat("Counts matrix dimensions after subsetting: ", dim(counts_matrix), "\n")
cat("Data matrix dimensions after subsetting: ", dim(data_matrix), "\n")
cat("Scale data matrix dimensions after subsetting: ", dim(scale_data_matrix), "\n")

# Extract metadata
raw_metadata <- neutrophil_subset@meta.data

# Create a new Seurat object
new_neutrophil_subset <- CreateSeuratObject(counts = counts_matrix, meta.data = raw_metadata)

# Add the other data layers back to the Seurat object
new_neutrophil_subset <- SetAssayData(new_neutrophil_subset, slot = "data", new.data = data_matrix)
new_neutrophil_subset <- SetAssayData(new_neutrophil_subset, slot = "scale.data", new.data = scale_data_matrix)

# Proceed with normalization and downstream analysis
new_neutrophil_subset <- NormalizeData(new_neutrophil_subset)
new_neutrophil_subset <- FindVariableFeatures(new_neutrophil_subset, selection.method = "vst", nfeatures = 2000)
new_neutrophil_subset <- ScaleData(new_neutrophil_subset, features = rownames(new_neutrophil_subset))
new_neutrophil_subset <- RunPCA(new_neutrophil_subset, features = VariableFeatures(object = new_neutrophil_subset))
new_neutrophil_subset <- FindNeighbors(new_neutrophil_subset, dims = 1:10)
new_neutrophil_subset <- FindClusters(new_neutrophil_subset, resolution = 0.5)
new_neutrophil_subset <- RunUMAP(new_neutrophil_subset, dims = 1:10)

# Plot UMAP to visualize the sub-clusters within the neutrophil cluster
DimPlot(new_neutrophil_subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Find markers for each sub-cluster
neutrophil_markers <- FindAllMarkers(new_neutrophil_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View the top markers for each sub-cluster
top_markers <- neutrophil_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
print(top_markers)

# Load required libraries
library(Seurat)
library(dplyr)

# Assuming you have already created and processed the Seurat object `new_neutrophil_subset`
# and found markers using `FindAllMarkers`
# `neutrophil_markers` should contain the results from `FindAllMarkers`

# Extract the top 10 markers for each sub-cluster
top_markers <- neutrophil_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Print the top markers to check
print(top_markers)

# Ensure the Seurat object is correctly named
MUTifc_seurat_object <- new_neutrophil_subset

# Extract the top genes for visualization in the heatmap
top_genes <- top_markers$gene

# Create a heatmap for the top markers
DoHeatmap(MUTifc_seurat_object, features = top_genes) + NoLegend()


# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Assuming the Seurat object `new_neutrophil_subset` has been processed and is ready for plotting

# 1. Feature Plot
FeaturePlot(new_neutrophil_subset, features = c("GeneA")) + ggtitle("GeneA Expression")

# 2. Violin Plot
VlnPlot(new_neutrophil_subset, features = c("GeneA", "GeneB")) + ggtitle("Gene Expression Distribution")

# 3. Ridge Plot
RidgePlot(new_neutrophil_subset, features = c("GeneA")) + ggtitle("GeneA Expression Density")

# 4. Dot Plot
# Assuming `top_genes` has been defined previously
DotPlot(new_neutrophil_subset, features = top_genes) + ggtitle("Top Markers Expression")

# 5. Cell Type Proportions
cell_proportions <- as.data.frame(table(Idents(new_neutrophil_subset)))
colnames(cell_proportions) <- c("Cluster", "Count")

ggplot(cell_proportions, aes(x = Cluster, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ggtitle("Cell Type Proportions")

##Move on to analyse HSCs 
library(Seurat)
library(ggplot2)
library(dplyr)

# Subset the Seurat object to include only "Stem cells"
stem_cell_cluster <- subset(MUTifc_seurat_object, idents = "Stem cells")

# Verify the subset
print(stem_cell_cluster)

# Normalize, find variable features, and scale the data
stem_cell_cluster <- NormalizeData(stem_cell_cluster)
stem_cell_cluster <- FindVariableFeatures(stem_cell_cluster, selection.method = "vst", nfeatures = 2000)
stem_cell_cluster <- ScaleData(stem_cell_cluster)

# Perform PCA
stem_cell_cluster <- RunPCA(stem_cell_cluster, features = VariableFeatures(object = stem_cell_cluster))

# Find neighbors and perform clustering
stem_cell_cluster <- FindNeighbors(stem_cell_cluster, dims = 1:10)
stem_cell_cluster <- FindClusters(stem_cell_cluster, resolution = 0.5)  # Adjust resolution as needed

# Run UMAP for visualization
stem_cell_cluster <- RunUMAP(stem_cell_cluster, dims = 1:10)

# Identify marker genes for each subcluster
stem_cell_subcluster_markers <- FindAllMarkers(stem_cell_cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View top markers for each subcluster
head(stem_cell_subcluster_markers)

# Select top 10 marker genes for each cluster
top10 <- stem_cell_subcluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Create the heatmap
DoHeatmap(stem_cell_cluster, features = top10$gene) + NoLegend()


save.image(file = "scRNAseq_2_MUTifc.RData", version = NULL, ascii = FALSE, compress = TRUE, safe = TRUE)
