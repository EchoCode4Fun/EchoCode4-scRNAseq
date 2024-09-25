library(Matrix)

# Define the path to your files
data_dir <- "scRNAseq3/"

# Read the matrix file
control_matrix <- readMM(file.path(data_dir, "GSM5916440_17422filteredmatrix.mtx.gz"))

# Read the barcodes and features
control_barcodes <- readLines(file.path(data_dir, "GSM5916440_17422filteredbarcodes.tsv.gz"))
control_features <- readLines(file.path(data_dir, "GSM5916440_17422filteredfeatures.tsv.gz"))

# Split the character vector by tabs
split_features <- strsplit(control_features, "\t")

# Convert the list of split elements into a data frame
features_df <- do.call(rbind, split_features)
features_df <- as.data.frame(features_df, stringsAsFactors = FALSE)

# Assign column names for clarity (adjust the names as needed)
colnames(features_df) <- c("ID", "GeneSymbol", "Description")

# Check the structure of the resulting data frame
str(features_df)

# Extract gene symbols from the second column
control_gene_symbols <- features_df$GeneSymbol

# Print the gene symbols to verify
print(control_gene_symbols)

# Convert barcodes and gene symbols to proper format
colnames(control_matrix) <- control_barcodes
rownames(control_matrix) <- control_gene_symbols

# Make row names unique by appending a suffix
rownames(control_matrix) <- make.unique(rownames(control_matrix))

# Verify that all row names are now unique
print(any(duplicated(rownames(control_matrix))))

library(Seurat)

# Create a Seurat object
control_seurat_object <- CreateSeuratObject(counts = control_matrix, project = "scRNAseq", min.cells = 3, min.features = 200)
HNifc_seurat_object[["percent.mt"]] <- PercentageFeatureSet(HNifc_seurat_object, pattern = "^MT-")

# Normalize the data (account for library size)
control_seurat_object <- NormalizeData(control_seurat_object)

# Find variable features (choose either disp or vst;usually nfeature=2000)
control_seurat_object <- FindVariableFeatures(control_seurat_object, selection.method = "vst", nfeatures = 2000) #450
VariableFeaturePlot(control_seurat_object)
selected_genes <- VariableFeatures(control_seurat_object)
print(selected_genes)

# Scale the data (account for sequencing depth)
control_seurat_object <- ScaleData(control_seurat_object, features = rownames(control_seurat_object))

# Perform PCA
control_seurat_object <- RunPCA(control_seurat_object, features = VariableFeatures(object = control_seurat_object))

# **ELBOW PLOT FOR DIMENSION DETERMINATION**
# Calculate the percentage of variance explained by each PC
pca_variance <- control_seurat_object@reductions$pca@stdev / sum(control_seurat_object@reductions$pca@stdev) * 100

# Create the elbow plot
plot(pca_variance, type = "b", xlab = "Principal Component", ylab = "Percentage of Variance Explained", main = "Elbow Plot for PCA")

# (Optionally) Add a line for visual guidance
lines(pca_variance, type = "b")

#The results return several options
#Option1: PC 4 or 5: Start by trying PC 4 or 5, as the first elbow suggests. Examine your clustering results carefully.
#Option2: PC 10: Next, try using PC 10 to see if including those additional components improves your clustering in a meaningful way (e.g., more refined clusters, better separation of cell types).
#Option3: Intermediate Values: You could even explore intermediate values like PC 7 or 8 to see if they strike a good balance.
#the first one tried is 7. 

# Find neighbors to determine cell similarity
control_seurat_object <- FindNeighbors(control_seurat_object, dims = 1:7)

# Cluster the cells (start low and go high)
control_seurat_object <- FindClusters(control_seurat_object, resolution = 0.4)

# Run UMAP for visualization
control_seurat_object <- RunUMAP(control_seurat_object, dims = 1:7)
library(ggplot2)
DimPlot(control_seurat_object, reduction = "umap", label = TRUE) + NoLegend() + ggtitle("Control Whole Bone Marrow_WOAnnotation")

##   II: Assess the quality of your clustering --- run Silhoutte analysis
# 1. Extract Cluster Labels and UMAP Embeddings
cluster_labels <- Idents(control_seurat_object)  # Get cluster assignments
umap_embeddings <- control_seurat_object@reductions$umap@cell.embeddings  # Extract UMAP coordinates

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
seq3_control_markers <- FindAllMarkers(control_seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Run ROC test to identity those subtle markers

table(control_seurat_object@active.ident)
##Visulization of the markers 
VlnPlot(control_seurat_object, features = c("Cd34", "Flt3","Ly6a","Ly6e","Slamf1","Slamf2","Kit","Flk2")) #to identify HSCs
VlnPlot(control_seurat_object, features = c("Kit","Slamf1","Ly6a","Ly6e","Cd34", "Flt3","Cd48")) #assess LT-HSC (lack of Lin,cd34, cd135/Flt3, cd48; enriched for Sca-1+ and c-Kit+ and cd150/slamf1)

VlnPlot(control_seurat_object, features = c("Ly6a","Ly6e","Kit","Slamf1","Cd48","Flt3")) #assess ST-HSC (lack of lin and Flt3)
#assess MPP (lack of lin, but enriched for Cd34, Flt3, and Sca-1 and c-kit)
VlnPlot(control_seurat_object, features = c("Ly6a","Ly6e","Kit")) 
VlnPlot(control_seurat_object, features = c("Ly6a","Ly6e","Kit","Itga2b")) 

VlnPlot(control_seurat_object, features = c("Cd34", "Cd105","Endoglin","Itga2b")) #to identify Megakaryocyte-Erythroid progenitors
VlnPlot(control_seurat_object, features = c("Cd34", "Fcgr2b", "Fcgr3","Csf3r", "Kit")) #to identify Granulocyte-Monocyte progenitors
VlnPlot(control_seurat_object, features = c("Cd34", "Flt3","Fcgr3")) #common myeloid progenitors
VlnPlot(control_seurat_object, features = c("Cd34", "Flt3","Flk2","Cd117")) #multipotent progenitors
VlnPlot(control_seurat_object, features = c("Flt3","Il7r")) #Common Lymphoid progenitors
#neutrophils
VlnPlot(control_seurat_object, features = c("S100a8", "S100a9","Elane","Mpo","Ly6g","Cxcr2","Cd11b")) #to identify neutrophils
VlnPlot(control_seurat_object, features = c("S100a8", "S100a9","Elane","Mpo","Ly6g","Adgre1","Cxcr2","Cd11b","c-fms","Ccr2","Cd115"))
#macrophages
VlnPlot(control_seurat_object, features = c("Adgre1", "Csf1r","Cd68","Itgam","Mertk","H2-Ab1","H2-Ea")) #to identify monocytes macrophages
#B cells
#B
VlnPlot(control_seurat_object, features = c("Ptprc","Cd19")) #pan-B cell markers
VlnPlot(control_seurat_object, features = c("Il7r","Cd24a","Spn")) #pro-B cell marker
VlnPlot(control_seurat_object, features = c("Cd79a","Cd79b","Vpreb1", "Igll1")) #pre-B cell marker
VlnPlot(control_seurat_object, features = c("Ighd","Il2ra")) #mature B cells


VlnPlot(control_seurat_object, features = c("Cd4", "Cd8","Cd3e")) #identify T cells
VlnPlot(control_seurat_object, features = c("Klrb1c", "Ncr1")) #identify NK cells
VlnPlot(control_seurat_object, features = c("Cd11c", "H2-Ab1")) #identify DC cells
VlnPlot(control_seurat_object, features = c("Pf4", "Itga2b","Vwf")) #identify Megakaryocytes
VlnPlot(control_seurat_object, features = c("Hba-a1", "Hba-a2","Hbb-bs","Hbb-bt")) #identify Erythrocytes 
VlnPlot(control_seurat_object, features = c("Nt5e","Thy1","Eng","Runx2","Alpl","Col1a1")) #identify stromal cells
VlnPlot(control_seurat_object, features = c("Pecam1","Cdh5","Vwf","Lepr","Kitl")) #endothelial cells 
VlnPlot(control_seurat_object, features = c("Lin","Ly6a","Ckit"))
VlnPlot(control_seurat_object, features = c("Cd47","Thbs1"))

#Vinplot: plot raw counts 
VlnPlot(control_seurat_object, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

#Featureplot: shows global distribution of a series of markers
FeaturePlot(control_seurat_object, features = c("Ms4a1", "Gnly", "Cd3e", "Cd14", "Fcer1a", "Fcgr3a", "Lyz", "Ppbp", "Cd8a"))

#Heatmap: show the top 10 markers (a good way to annotate clusters)
HNifc_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(control_seurat_object, features = top10$gene) + NoLegend() + theme(axis.text.y = element_text(size = 6))




# Extract the top 10 markers for each cluster and save to Excel
library(dplyr)
top10_markers <- control_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)%>% arrange(cluster, desc(avg_log2FC))
print(top10_markers)

library(writexl)
write_xlsx(top10_marker, "scRNAseq#3_control_top10_markers.xlsx")


##Manually add annotation
# Combine and annotate clusters
cluster_annotations <- c(
  "0" = "GMP",
  "1" = "CD41+ HSC",
  "2" = "B Cells",
  "3" = "Neutrophils",
  "4" = "Macrophages Monocytes",
  "5" = "GMP",
  "6" = "CD41+ HSC",
  "7" = "CMP",
  "8" = "Pre-B Cells",
  "9" = "Macrophages Monocytes",
  "10" = "Erythrocytes",
  "11" = "CD41+ HSC",
  "12" = "NK Cells",
  "13" = "Neutrophils",
  "14" = "Pro-B Cells",
  "15" = "Neutrophils"
)

#for CellChat 
cluster_annotations <- c(
  "0" = "HSPCs",
  "1" = "HSPCs",
  "2" = "B Cells",
  "3" = "Neutrophils",
  "4" = "Macrophages Monocytes",
  "5" = "HSPCs",
  "6" = "HSPCs",
  "7" = "HSPCs",
  "8" = "Pre-B Cells",
  "9" = "Macrophages Monocytes",
  "10" = "Erythrocytes",
  "11" = "HSPCs",
  "12" = "NK Cells",
  "13" = "Neutrophils",
  "14" = "Pro-B Cells",
  "15" = "Neutrophils"
)

# Apply annotations to the Seurat object
control_seurat_object <- RenameIdents(control_seurat_object, cluster_annotations)

# Visualize the annotated clusters
DimPlot(control_seurat_object, reduction = "umap", label = TRUE) + NoLegend() + ggtitle("Control Whole Bone Marrow") + theme(text = element_text(size = 10, face = "bold"))




##Using ROGUE to asses the heterogeneity of the cluster
##If ROGUE score < = 0.9, it means high heterogeneity and subcluster exists
#to look as HSPC
hspc_cells <- subset(control_seurat_object, idents = "Hematopoietic Stem/Progenitor Cells")
hspc_expression_matrix <- GetAssayData(hspc_cells, slot = "counts")

ent.res <- SE_fun(hspc_expression_matrix)
head(ent.res)

SEplot(ent.res)

rogue.value <- CalculateRogue(ent.res, platform = "UMI")
rogue.value

#to look at macrophages monocytes  
macmono_cells <- subset(control_seurat_object, idents = "Macrophages Monocytes")
macmono_expression_matrix <- GetAssayData(macmono_cells, slot = "counts")

ent.res <- SE_fun(macmono_expression_matrix)
head(ent.res)

SEplot(ent.res)

rogue.value <- CalculateRogue(ent.res, platform = "UMI")
rogue.value

#to look at neutrophis 
aneu_cells <- subset(control_seurat_object, idents = "Antimicrobial Neutrophils")
aneu_expression_matrix <- GetAssayData(aneu_cells, slot = "counts")

ent.res <- SE_fun(aneu_expression_matrix)
head(ent.res)

SEplot(ent.res)

rogue.value <- CalculateRogue(ent.res, platform = "UMI")
rogue.value


###Automatic addition of annotation
library(SingleR)
library(celldex)

ref_data <- celldex::MouseRNAseqData()

# Extract gene symbols from test genes
test_genes <- rownames(control_seurat_object)
ref_genes <- rownames(ref_data)
ref_genes <- toupper(ref_genes)
test_genes <- toupper(rownames(control_seurat_object))

# Find common genes
common_genes <- intersect(test_genes, ref_genes)
if (length(common_genes) == 0) {
  stop("No common genes found between test and reference datasets.")
} else {
  print(paste(length(common_genes), "common genes found."))
}

#subset the Seurat and reference to the common genes 
counts <- GetAssayData(control_seurat_object, layer = "counts")[rownames(control_seurat_object) %in% common_genes, ]
ref_data_subset <- ref_data[rownames(ref_data) %in% common_genes, , drop = FALSE]

# Perform SingleR annotation
singleR_results <- SingleR(test = counts, ref = ref, labels = ref$label.main, clusters = clusters)

# Create a named vector of SingleR annotations
cluster_annotations <- singleR_results$labels
names(cluster_annotations) <- levels(Idents(control_seurat_object))

# Rename the identities in the Seurat object using the named vector
control_seurat_object <- RenameIdents(control_seurat_object, cluster_annotations)

# Plot UMAP with annotations
DimPlot(control_seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() + ggtitle("Control Whole Bone Marrow_Annotated")





### NOW STUDY THE NEUTROPHILS SUBCLUSTER

# Subset the Neutrophil cluster
control_neutrophils <- subset(control_seurat_object, idents = "Neutrophils")

# Normalize the subset data
control_neutrophils <- NormalizeData(control_neutrophils)

# Find variable features in the subset
control_neutrophils <- FindVariableFeatures(control_neutrophils, selection.method = "vst", nfeatures = 2000)

# Scale the subset data
control_neutrophils <- ScaleData(control_neutrophils, features = rownames(control_neutrophils))

# Perform PCA on the subset
control_neutrophils <- RunPCA(control_neutrophils, features = VariableFeatures(object = control_neutrophils))

# Find neighbors in the subset
control_neutrophils <- FindNeighbors(control_neutrophils, dims = 1:10)

# Cluster the subset
control_neutrophils <- FindClusters(control_neutrophils, resolution = 0.5)

# Run UMAP for the subset
control_neutrophils <- RunUMAP(control_neutrophils, dims = 1:10)

# Plot UMAP for neutrophil subclusters
DimPlot(control_neutrophils, reduction = "umap", label = TRUE, label.size = 5) + 
  ggtitle("Control Neutrophil Subclusters")

# Find markers for each neutrophil subcluster
control_neutrophil_markers <- FindAllMarkers(control_neutrophils, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Extract the top 10 markers for each subcluster
top10_control_neutrophil_markers <- control_neutrophil_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


# Write the top 10 markers for neutrophil subclusters to a CSV file
write.csv(top10_control_neutrophil_markers, "scRNAseq2/top10_neutrophil_markers.csv", row.names = FALSE)

# Draw a heatmap of the top markers for neutrophil subclusters
DoHeatmap(control_neutrophils, features = top10_control_neutrophil_markers$gene) + 
  ggtitle("Heatmap of Top Markers for Control Neutrophil Subclusters")

# Convert gene symbols to Entrez IDs
gene_symbols <- top10_control_neutrophil_markers$gene
gene_entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Merge with top10_control_neutrophil_markers to ensure we have cluster information
top10_control_neutrophil_markers <- merge(top10_control_neutrophil_markers, gene_entrez_ids, by.x = "gene", by.y = "SYMBOL")

go_results_list <- list()

for (cluster in unique(top10_control_neutrophil_markers$cluster)) {
  cluster_genes <- top10_control_neutrophil_markers %>% filter(cluster == !!cluster) %>% pull(ENTREZID)
  
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

library(ggplot2)

for (cluster in names(go_results_list)) {
  go_results <- go_results_list[[cluster]]
  
  if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
    print(barplot(go_results, showCategory = 10) + ggtitle(paste("GO Enrichment for Subcluster", cluster)))
  }
}

##### Move on to subset HSC

# Subset Stem cells
control_stem_cells <- subset(control_seurat_object, idents = "Hematopoietic Stem/Progenitor Cells")

# Normalize the data
control_stem_cells <- NormalizeData(control_stem_cells)

# Identify the most variable features
control_stem_cells <- FindVariableFeatures(control_stem_cells)

# Scale the data
control_stem_cells <- ScaleData(control_stem_cells)

# Perform PCA
control_stem_cells <- RunPCA(control_stem_cells)

# Perform UMAP dimensional reduction
control_stem_cells <- RunUMAP(control_stem_cells, dims = 1:10)

# Identify subclusters within Stem cells
control_stem_cells <- FindClusters(control_stem_cells, resolution = 0.5)  # Adjust resolution as needed

# Find markers for each subcluster
stem_cell_markers <- FindAllMarkers(control_stem_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Extract the top 10 markers for each subcluster
top10_stem_cell_markers <- stem_cell_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Draw UMAP plot
DimPlot(control_stem_cells, reduction = "umap", label = TRUE, pt.size = 0.5) + ggtitle("UMAP of Control HSCs")

# Draw Heatmap
top10_genes <- unique(top10_stem_cell_markers$gene)
DoHeatmap(control_stem_cells, features = top10_genes) + ggtitle("Heatmap of Top 10 Markers per Control HSC Subcluster")


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

# Visualize GO enrichment results
for (cluster in names(go_results_list)) {
  go_results <- go_results_list[[cluster]]
  
  if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
    print(barplot(go_results, showCategory = 10) + ggtitle(paste("GO Enrichment for Subcluster", cluster)))
  }
}

######################Now for infected group 
# Define the path to your files
data_dir <- "scRNAseq3/"

# Read the matrix file
MAifc_matrix <- readMM(file.path(data_dir, "GSM5916441_17423filteredmatrix.mtx.gz"))

# Read the barcodes and features
MAifc_barcodes <- readLines(file.path(data_dir, "GSM5916441_17423filteredbarcodes.tsv.gz"))
MAifc_features <- readLines(file.path(data_dir, "GSM5916441_17423filteredfeatures.tsv.gz"))

# Split the character vector by tabs
split_features <- strsplit(MAifc_features, "\t")

# Convert the list of split elements into a data frame
features_df <- do.call(rbind, split_features)
features_df <- as.data.frame(features_df, stringsAsFactors = FALSE)

# Assign column names for clarity (adjust the names as needed)
colnames(features_df) <- c("ID", "GeneSymbol", "Description")

# Check the structure of the resulting data frame
str(features_df)

# Extract gene symbols from the second column
MAifc_gene_symbols <- features_df$GeneSymbol

# Print the gene symbols to verify
print(MAifc_gene_symbols)

# Convert barcodes and gene symbols to proper format
colnames(MAifc_matrix) <- MAifc_barcodes
rownames(MAifc_matrix) <- MAifc_gene_symbols

# Make row names unique by appending a suffix
rownames(MAifc_matrix) <- make.unique(rownames(MAifc_matrix))

# Verify that all row names are now unique
print(any(duplicated(rownames(MAifc_matrix))))

# Create a Seurat object
MAifc_seurat_object <- CreateSeuratObject(counts = MAifc_matrix, project = "scRNAseq", min.cells = 3, min.features = 200)

##Quality check
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
MAifc_seurat_object[["percent.mt"]] <- PercentageFeatureSet(MAifc_seurat_object, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(MAifc_seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(MAifc_seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(MAifc_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
MAifc_seurat_object <- subset(MAifc_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

# Normalize the data
MAifc_seurat_object <- NormalizeData(MAifc_seurat_object)

# Find variable features
MAifc_seurat_object <- FindVariableFeatures(MAifc_seurat_object, selection.method = "vst", nfeatures = 2000)

# Scale the data
MAifc_seurat_object <- ScaleData(MAifc_seurat_object, features = rownames(MAifc_seurat_object))

# Perform PCA
MAifc_seurat_object <- RunPCA(MAifc_seurat_object, features = VariableFeatures(object = MAifc_seurat_object))

#Visualiza and examine PCA
VizDimLoadings(MAifc_seurat_object, dims = 1:2, reduction = "pca")
DimPlot(MAifc_seurat_object, reduction = "pca") + NoLegend()
DimHeatmap(MAifc_seurat_object, dims = 1, cells = 500, balanced = TRUE) #can change the dim numebr
DimHeatmap(MAifc_seurat_object, dims = 1:15, cells = 500, balanced = TRUE) #this shows all dims


ElbowPlot(MAifc_seurat_object)

# Find neighbors
MAifc_seurat_object <- FindNeighbors(MAifc_seurat_object, dims = 1:7) #1:7

# Cluster the cells
MAifc_seurat_object <- FindClusters(MAifc_seurat_object, resolution = 0.4) #0.4

# Run UMAP for visualization
MAifc_seurat_object <- RunUMAP(MAifc_seurat_object, dims = 1:7)
DimPlot(MAifc_seurat_object, reduction = "umap", label = TRUE) + NoLegend() + ggtitle("M.avium Infected Whole Bone Marrow_WOAnnotation") +
  theme(text = element_text(size = 10, face = "bold"))




##   II: Assess the quality of your clustering --- run Silhoutte analysis
# 1. Extract Cluster Labels and UMAP Embeddings
cluster_labels <- Idents(MAifc_seurat_object)  # Get cluster assignments
umap_embeddings <- MAifc_seurat_object@reductions$umap@cell.embeddings  # Extract UMAP coordinates

# 2. Calculate Silhouette Scores
library(cluster)  # Load the 'cluster' package for silhouette() function
dist_matrix <- dist(umap_embeddings)  # Calculate distance matrix from UMAP embeddings
sil <- silhouette(as.numeric(cluster_labels), dist_matrix) # Perform silhouette analysis

# 3. Visualize and Interpret Silhouette Plot
plot(sil, main = "Silhouette Plot for Seurat Clustering") 
text(0.5, 0.95, paste("Average Silhouette Score:", round(avg_sil_score, 3)), cex = 1.2)

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
MAifc_markers <- FindAllMarkers(MAifc_seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Run ROC test to identity those subtle markers

##Visulization of the markers 
#Vinplot: shows expression probability 
VlnPlot(MAifc_seurat_object, features = c("Cd34", "Flt3","Ly6a","Ly6e","Slamf1","Slamf2","Kit","Flk2")) #to identify HSCs
VlnPlot(MAifc_seurat_object, features = c("Ly6a","Ly6e","Cd34", "Flt3","Cd48","Ly6a/e","Kit","Slamf1")) #assess LT-HSC (lack of Lin,cd34, cd135/Flt3, cd48; enriched for Sca-1+ and c-Kit+ and cd150/slamf1)

VlnPlot(MAifc_seurat_object, features = c("Ly6a","Ly6e","Kit","Slamf1","Cd48","Flt3")) #assess ST-HSC (lack of lin and Flt3)
#assess MPP (lack of lin, but enriched for Cd34, Flt3, and Sca-1 and c-kit)
VlnPlot(MAifc_seurat_object, features = c("Ly6a","Ly6e","Kit")) 
VlnPlot(MAifc_seurat_object, features = c("Ly6a","Ly6e","Kit","Itga2b")) 

VlnPlot(MAifc_seurat_object, features = c("Ptprc","Prom1","Slamf1","Cd48","Kit")) #assess HSCs subtypes 
VlnPlot(MAifc_seurat_object, features = c("Cd34", "Cd105","Endoglin","Itga2b")) #to identify Megakaryocyte-Erythroid progenitors
VlnPlot(MAifc_seurat_object, features = c("Cd34", "Fcgr2b", "Fcgr3","Csf3r", "Kit")) #to identify Granulocyte-Monocyte progenitors
VlnPlot(MAifc_seurat_object, features = c("Cd34", "Flt3","Fcgr3")) #common myeloid progenitors
VlnPlot(MAifc_seurat_object, features = c("Cd34", "Flt3","Flk2","Cd117")) #multipotent progenitors
VlnPlot(MAifc_seurat_object, features = c("Flt3","Il7r")) #Common Lymphoid progenitors
#Neutrophils
VlnPlot(MAifc_seurat_object, features = c("S100a8", "S100a9","Elane","Mpo","Ly6g","Cxcr2","Cd11b")) #to identify neutrophils
VlnPlot(MAifc_seurat_object, features = c("S100a8", "S100a9","Elane","Mpo","Ly6g","Adgre1","Cxcr2","Cd11b","c-fms","Ccr2","Cd115"))
#macrophages
VlnPlot(MAifc_seurat_object, features = c("Adgre1", "Csf1r","Cd68","Itgam","Mertk","H2-Ab1","H2-Ea")) #to identify monocytes macrophages
VlnPlot(MAifc_seurat_object, features = c("Cd14", "Cd80","Cd86","Mrc1","Nos2","Arg1"))
#B
VlnPlot(MAifc_seurat_object, features = c("Ptprc","Cd19")) #pan-B cell markers
VlnPlot(MAifc_seurat_object, features = c("Il7r","Cd24a","Spn")) #pro-B cell marker
VlnPlot(MAifc_seurat_object, features = c("Cd79a","Cd79b","Vpreb1", "Igll1")) #pre-B cell marker
VlnPlot(MAifc_seurat_object, features = c("Ighd","Il2ra")) #mature B cells
#T
VlnPlot(MAifc_seurat_object, features = c("Cd4", "Cd8","Cd3e")) #identify T cells
#NK
VlnPlot(MAifc_seurat_object, features = c("Klrb1c", "Ncr1")) #identify NK cells
VlnPlot(MAifc_seurat_object, features = c("Cd11c", "H2-Ab1")) #identify DC cells
VlnPlot(MAifc_seurat_object, features = c("Siglec-F", "Ccr3")) #identify Eosinophils cells
VlnPlot(MAifc_seurat_object, features = c("Fcer1a", "Cpa3","Fcer1b")) #identify Basophils cells
VlnPlot(MAifc_seurat_object, features = c("Kit", "Cpa3","FcεRI")) #identify Mast cells
VlnPlot(MAifc_seurat_object, features = c("Pf4", "Itga2b","Vwf")) #identify Megakaryocytes
VlnPlot(MAifc_seurat_object, features = c("Hba-a1", "Hba-a2","Hbb-bs","Hbb-bt")) #identify Erythrocytes 
VlnPlot(MAifc_seurat_object, features = c("Nt5e","Thy1","Eng","Runx2","Alpl","Col1a1")) #identify stromal cells
VlnPlot(MAifc_seurat_object, features = c("Pecam1","Cdh5","Vwf","Lepr","Kitl")) #endothelial cells 
VlnPlot(MAifc_seurat_object, features = c("Lin","Ly6a","Ckit"))
VlnPlot(MAifc_seurat_object, features = c("Cd47","Thbs1"))

#Vinplot: plot raw counts 
VlnPlot(MAifc_seurat_object, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

#Featureplot: shows global distribution of a series of markers
FeaturePlot(MAifc_seurat_object, features = c("Ms4a1", "Gnly", "Cd3e", "Cd14", "Fcer1a", "Fcgr3a", "Lyz", "Ppbp", "Cd8a"))

#Heatmap: show the top 10 markers (a good way to annotate clusters)
MAifc_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(MAifc_seurat_object, features = top10$gene) + NoLegend() + theme(axis.text.y = element_text(size = 6))


# Extract the top 10 markers for each cluster
library(dplyr)
MAifc_top10_markers <- MAifc_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
print(MAifc_top10_markers)

#export the markers to an excel file
library(openxlsx)
write.xlsx(MAifc_top10_markers, "MAifc_top10_markers.xlsx")



#Manually annotate the cluster
# Combine and annotate clusters
# Define the cluster annotations
cluster_annotations <- c(
  "0" = "GMP",
  "1" = "Neutrophils",
  "2" = "LT-HSC",
  "3" = "GMP",
  "4" = "B Cells",
  "5" = "CD41+ HSC",
  "6" = "CD41+ HSC",
  "7" = "Macrophages Monocytes",
  "8" = "Macrophages Monocytes",
  "9" = "Erythrocytes",
  "10" = "Pro-B Cells",
  "11" = "Macrophages Monocytes",
  "12" = "ST-HSC",
  "13" = "NK Cells ",
  "14" = "Pre-B Cells"
)

#for CellChat 
cluster_annotations <- c(
  "0" = "HSPCs",
  "1" = "Neutrophils",
  "2" = "HSPCs",
  "3" = "HSPCs",
  "4" = "B Cells",
  "5" = "HSPCs",
  "6" = "HSPCs",
  "7" = "Macrophages Monocytes",
  "8" = "Macrophages Monocytes",
  "9" = "Erythrocytes",
  "10" = "Pro-B Cells",
  "11" = "Macrophages Monocytes",
  "12" = "HSPCs",
  "13" = "NK Cells",
  "14" = "Pre-B Cells"
)

# Apply annotations to the Seurat object
MAifc_seurat_object <- RenameIdents(MAifc_seurat_object, cluster_annotations)

# Visualize the annotated clusters
library(ggplot2)
library(ggrepel) # For non-overlapping labels

DimPlot(MAifc_seurat_object, 
        reduction = "umap", 
        label = TRUE) + 
  NoLegend() + 
  ggtitle("MAifc Whole Bone Marrow") + 
  theme(text = element_text(size = 10, face = "bold")) +
  geom_text_repel( # Use ggrepel for non-overlapping labels
    aes(label = ident), # Assuming 'cell.type' is your metadata column
    size = 2, # Adjust label size as needed
    box.padding = 0.5, # Adjust padding around labels
    point.padding = 0.2, # Adjust padding between labels and points 
    max.overlaps = Inf # Prevent labels from being hidden due to overlaps
  ) 


DimPlot(MAifc_seurat_object, 
        reduction = "umap", 
        label = TRUE, 
        pt.size = 0.3,
        label.size = 4) + # Adjust the label size here
  NoLegend() + 
  ggtitle("M.avium Infected Whole Bone Marrow_Annotated")



###Automatic addtion of annotation

# Extract gene symbols from test genes
test_genes <- rownames(MAifc_seurat_object)
ref_genes <- toupper(ref_genes)

# Find common genes
common_genes <- intersect(test_genes, ref_genes)
if (length(common_genes) == 0) {
  stop("No common genes found between test and reference datasets.")
} else {
  print(paste(length(common_genes), "common genes found."))
}

# Subset the Seurat object and reference data to the common genes
counts <- GetAssayData(MAifc_seurat_object, slot = "counts")[rownames(MAifc_seurat_object) %in% common_genes, ]
ref <- ref[rownames(ref) %in% common_genes, , drop = FALSE]
common_genes <- intersect(rownames(counts), rownames(ref))
counts <- counts[common_genes, ]
ref <- ref[common_genes, , drop = FALSE]

# Extract cluster labels
clusters <- Idents(MAifc_seurat_object)

# Perform SingleR annotation
singleR_results <- SingleR(test = counts, ref = ref, labels = ref$label.main, clusters = clusters)

# Create a named vector of SingleR annotations
cluster_annotations <- singleR_results$labels
names(cluster_annotations) <- levels(Idents(MAifc_seurat_object))

# Rename the identities in the Seurat object using the named vector
MAifc_seurat_object <- RenameIdents(MAifc_seurat_object, cluster_annotations)

# Plot UMAP with annotations
DimPlot(MAifc_seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() + ggtitle("M.avium Infected Whole Bone Marrow_Annotated")


##### Move on to subset HSC

# Subset Stem cells
MAifc_stem_cells <- subset(MAifc_seurat_object, idents = "Hematopoietic Stem/Progenitor Cells")

# Normalize the data
MAifc_stem_cells <- NormalizeData(MAifc_stem_cells)

# Identify the most variable features
MAifc_stem_cells <- FindVariableFeatures(MAifc_stem_cells)

# Scale the data
MAifc_stem_cells <- ScaleData(MAifc_stem_cells)

# Perform PCA
MAifc_stem_cells <- RunPCA(MAifc_stem_cells)

# Perform UMAP dimensional reduction
MAifc_stem_cells <- RunUMAP(MAifc_stem_cells, dims = 1:10)

# Identify subclusters within Stem cells
MAifc_stem_cells <- FindClusters(MAifc_stem_cells, resolution = 0.5)  # Adjust resolution as needed

# Find markers for each subcluster
MAifc_stem_cell_markers <- FindAllMarkers(MAifc_stem_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Extract the top 10 markers for each subcluster
top10_MAifc_stem_cell_markers <- MAifc_stem_cell_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Draw UMAP plot
DimPlot(MAifc_stem_cells, reduction = "umap", label = TRUE, pt.size = 0.5) + ggtitle("UMAP of MAifc HSCs")

# Draw Heatmap
top10_genes <- unique(top10_MAifc_stem_cell_markers$gene)
DoHeatmap(MAifc_stem_cells, features = top10_genes) + ggtitle("Heatmap of Top 10 Markers per MAifc HSC Subcluster") +
  theme(axis.text.y = element_text(size = 5))  # Adjust the size as needed

# Create the heatmap and customize the y-axis font size
heatmap_plot <- DoHeatmap(MAifc_stem_cells, features = top10_genes) + 
  ggtitle("Heatmap of Top 10 Markers per MAifc HSC Subcluster") +
  theme(axis.text.y = element_text(size = 8))  # Adjust the size as needed

# Print the heatmap
print(heatmap_plot)





# Convert gene symbols to Entrez IDs
gene_symbols <- top10_MAifc_stem_cell_markers$gene
gene_entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Merge with top10_stem_cell_markers to ensure we have cluster information
top10_MAifc_stem_cell_markers <- merge(top10_MAifc_stem_cell_markers, gene_entrez_ids, by.x = "gene", by.y = "SYMBOL", all.x = TRUE)

# Rename ENTREZID column to avoid confusion
colnames(top10_MAifc_stem_cell_markers)[which(colnames(top10_MAifc_stem_cell_markers) == "ENTREZID.y")] <- "ENTREZID"

# Remove rows with NA values in the ENTREZID column
top10_MAifc_stem_cell_markers <- top10_MAifc_stem_cell_markers %>% filter(!is.na(ENTREZID))

# Perform GO enrichment analysis for each subcluster
go_results_list <- list()

for (cluster in unique(top10_MAifc_stem_cell_markers$cluster)) {
  cluster_genes <- top10_MAifc_stem_cell_markers %>% filter(cluster == !!cluster) %>% pull(ENTREZID)
  
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

# Visualize GO enrichment results
# Function to plot GO enrichment results for a specific cluster
plot_go_enrichment <- function(go_results, cluster) {
  if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
    dotplot(go_results, showCategory = 10, x = "GeneRatio") +
      scale_color_gradient(low = "blue", high = "red", name = "p.adjust") +
      ggtitle(paste("GO BP Enrichment for Ma-infected HSC Subcluster", cluster)) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  }
}


###################Compare control HSC with Maifc HSCs

#### NOW compare HSCs

# Add metadata to indicate the condition
control_stem_cells$condition <- "control"
MAifc_stem_cells$condition <- "MAifc"

# Merge the control and MAifc neutrophils into one combined object
HSC_combined <- merge(control_stem_cells, y = MAifc_stem_cells, add.cell.ids = c("control", "MAifc"))

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
de_genes_HSC <- FindMarkers(HSC_combined, ident.1 = "control", ident.2 = "MAifc", group.by = "condition")

# View the top differentially expressed genes
head(de_genes_HSC)

library(EnhancedVolcano)
## now for visualization 
# Create a volcano plot
EnhancedVolcano(de_genes_HSC,
                lab = rownames(de_genes_HSC),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Differential Expression of MAifc vs Control HSC',
                pCutoff = 0.05,
                FCcutoff = 0.25,
                pointSize = 1.0,
                labSize = 3.5)

# Select top 20 DE genes based on adjusted p-value
top_genes <- head(rownames(de_genes_HSC[order(de_genes_HSC$p_val_adj), ]), 20)

# Extract normalized expression data for the top DE genes
expr_data <- FetchData(HSC_combined, vars = top_genes)

# Add condition metadata
expr_data$condition <- HSC_combined$condition

# Create a heatmap
pheatmap(t(expr_data[, -ncol(expr_data)]), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         annotation_col = data.frame(condition = expr_data$condition),
         use_raster = FALSE)


library(clusterProfiler)
library(org.Mm.eg.db) 

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
  ggtitle("GO BP Enrichment Analysis for DE Genes: MAifc vs Control HSCs")



### Move on to study neutrophils from the infected group 


### NOW STUDY THE NEUTROPHILS SUBCLUSTER

# Subset the Neutrophil cluster
MAifc_neutrophils <- subset(MAifc_seurat_object, idents = "Neutrophils")

# Normalize the subset data
MAifc_neutrophils <- NormalizeData(MAifc_neutrophils)

# Find variable features in the subset
MAifc_neutrophils <- FindVariableFeatures(MAifc_neutrophils, selection.method = "vst", nfeatures = 2000)

# Scale the subset data
MAifc_neutrophils <- ScaleData(MAifc_neutrophils, features = rownames(MAifc_neutrophils))

# Perform PCA on the subset
MAifc_neutrophils <- RunPCA(MAifc_neutrophils, features = VariableFeatures(object = MAifc_neutrophils))

# Find neighbors in the subset
MAifc_neutrophils <- FindNeighbors(MAifc_neutrophils, dims = 1:10)

# Cluster the subset
MAifc_neutrophils <- FindClusters(MAifc_neutrophils, resolution = 0.5)

# Run UMAP for the subset
MAifc_neutrophils <- RunUMAP(MAifc_neutrophils, dims = 1:10)

# Plot UMAP for neutrophil subclusters
DimPlot(MAifc_neutrophils, reduction = "umap", label = TRUE, label.size = 5) + 
  ggtitle("post-MAifc Neutrophil Subclusters")

# Find markers for each neutrophil subcluster
MAifc_neutrophil_markers <- FindAllMarkers(MAifc_neutrophils, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Extract the top 10 markers for each subcluster
top10_MAifc_neutrophil_markers <- MAifc_neutrophil_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Draw a heatmap of the top markers for neutrophil subclusters
DoHeatmap(MAifc_neutrophils, features = top10_MAifc_neutrophil_markers$gene) + 
  ggtitle("Heatmap of Top Markers for post-MAifc Neutrophil Subclusters")

## To do GO enrichment analysis 

# Load necessary libraries
library(clusterProfiler)
library(dplyr)
library(ggplot2)

# Merge with top10_control_neutrophil_markers to ensure we have cluster information
top10_MAifc_neutrophil_markers <- merge(top10_MAifc_neutrophil_markers, gene_entrez_ids, by.x = "gene", by.y = "SYMBOL")

# Perform GO enrichment analysis for each cluster
go_results_list <- list()

for (cluster in unique(top10_MAifc_neutrophil_markers$cluster)) {
  cluster_genes <- top10_MAifc_neutrophil_markers %>%
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
      ggtitle(paste("GO Enrichment for post-MAifc Neutrophils Subcluster", cluster)) +
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



###Checking cell-cell interaction 

# Load necessary libraries
library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(CellChat)
library(dplyr)
library(celldex)  # For loading reference datasets

# Assuming you have a Seurat object loaded
# MAifc_seurat_object <- readRDS("MAifc_seurat_object_with_SingleR_annotations.rds")

# Convert Seurat object to SingleCellExperiment
sce <- as.SingleCellExperiment(MAifc_seurat_object)

# Load ImmGen reference data for SingleR
immgen_ref <- celldex::ImmGenData()

# Recompute SingleR annotations using ImmGen reference data
singleR_results <- SingleR(test = sce, ref = immgen_ref, labels = immgen_ref$label.main)

# Extract the SingleR annotations and ensure they have cell names
singleR_annotations <- singleR_results$labels
names(singleR_annotations) <- colnames(sce)

# Extract counts and metadata from Seurat object
counts <- GetAssayData(MAifc_seurat_object, slot = "counts")
meta <- MAifc_seurat_object@meta.data

# Align SingleR annotations with Seurat metadata
common_cells <- intersect(names(singleR_annotations), rownames(meta))
singleR_annotations <- singleR_annotations[common_cells]
meta <- meta[common_cells, ]

# Add SingleR annotations to the metadata
meta$singleR_annotations <- singleR_annotations

# Update the Seurat object with the new metadata
MAifc_seurat_object@meta.data <- meta

# Verify the updated metadata
head(meta)

# Create a CellChat object
MAifc_cellchat <- createCellChat(object = counts, meta = meta, group.by = "singleR_annotations")

# Set the database to use for CellChat (e.g., for mouse)
CellChatDB <- CellChatDB.mouse  # or CellChatDB.human for human data

# Use the selected database
MAifc_cellchat@DB <- CellChatDB

# Preprocess the data
MAifc_cellchat <- subsetData(MAifc_cellchat)  # Subset the expression data of signaling genes
MAifc_cellchat <- identifyOverExpressedGenes(MAifc_cellchat)
MAifc_cellchat <- identifyOverExpressedInteractions(MAifc_cellchat)

# Compute the communication probability
MAifc_cellchat <- computeCommunProb(MAifc_cellchat)

# Filter out the interactions with low probability
MAifc_cellchat <- filterCommunication(MAifc_cellchat, min.cells = 10)

# Infer cell-cell communication network
MAifc_cellchat <- computeCommunProbPathway(MAifc_cellchat)

# Aggregate the cell-cell communication network
MAifc_cellchat <- aggregateNet(MAifc_cellchat)

# Check the cell-cell communication network
print(MAifc_cellchat@net)

# Visualize the communication network using a circle plot
netVisual_circle(MAifc_cellchat, weight.scale = TRUE, label.edge = FALSE)

# Optionally, visualize the contribution of each ligand-receptor pair
netVisual_bubble(MAifc_cellchat, sources.use = 1, targets.use = 2)

##################################(only use the below)
#Set up control and MAifc cellchat

library(CellChat)
library(Seurat)
library(patchwork)

# Convert Seurat to CellChat
cellchat_MAifc <- createCellChat(object = MAifc_seurat_object, group.by = "ident")
cellchat_Control <- createCellChat(object = control_seurat_object, group.by = "ident")

# Set the database for communication
CellChatDB <- CellChatDB.mouse # or CellChatDB.mouse if using mouse data
cellchat_MAifc@DB <- CellChatDB
cellchat_Control@DB <- CellChatDB

# Subset the expression data of signaling genes
cellchat_MAifc <- subsetData(cellchat_MAifc)
cellchat_Control <- subsetData(cellchat_Control)

# Identify overexpressed genes and interactions
cellchat_MAifc <- identifyOverExpressedGenes(cellchat_MAifc)
cellchat_Control <- identifyOverExpressedGenes(cellchat_Control)

cellchat_MAifc <- identifyOverExpressedInteractions(cellchat_MAifc)
cellchat_Control <- identifyOverExpressedInteractions(cellchat_Control)

#Compute commun prob
cellchat_MAifc <- computeCommunProb(cellchat_MAifc)
cellchat_Control <- computeCommunProb(cellchat_Control)

# Filter out low probability interactions
cellchat_MAifc <- filterCommunication(cellchat_MAifc, min.cells = 10)
cellchat_Control <- filterCommunication(cellchat_Control, min.cells = 10)

#Infer cellular communication networks 
cellchat_MAifc <- computeCommunProbPathway(cellchat_MAifc)
cellchat_Control <- computeCommunProbPathway(cellchat_Control)

cellchat_MAifc <- aggregateNet(cellchat_MAifc)
cellchat_Control <- aggregateNet(cellchat_Control)


# Compare communication networks between MAifc and Control CellChat

#########The code below is for two datasets with different cell compositions 

library(CellChat)
library(patchwork)

##The goal is to make two datasets (slightly different) contain the same cell composition 
##For vastly different cell groups, some functions do not work # Check and print cell types in MAifc_cellchat
##The idea is to lift up the ones with less compositions to match the one with more 

# Check and print cell types in MAifc_cellchat
MAifc_cell_types <- levels(cellchat_MAifc@idents)
print("Cell types in MAifc_cellchat:")
print(MAifc_cell_types)

# Check and print cell types in Control_cellchat
control_cell_types <- levels(cellchat_Control@idents)
print("Cell types in Control_cellchat:")
print(control_cell_types)

# Define the combined cell types from both objects
group.new <- union(levels(cellchat_MAifc@idents), levels(cellchat_Control@idents))

# Lift up the cellchat objects to match the combined cell types
cellchat_MAifc <- liftCellChat(cellchat_MAifc, group.new)
cellchat_Control <- liftCellChat(cellchat_Control, group.new)

# Merge the CellChat objects
object.list <- list(MAifc = cellchat_MAifc, Control = cellchat_Control)
cellchat <- mergeCellChat(object.list, add.names = c("MAifc", "Control"), cell.prefix = TRUE)

# Print the merged CellChat object
print(cellchat)

# Check and print cell types in MAifc_cellchat
MAifc_cell_types <- levels(cellchat_MAifc@idents)
print("Cell types in MAifc_cellchat:")
print(MAifc_cell_types)

# Check and print cell types in Control_cellchat
control_cell_types <- levels(cellchat_Control@idents)
print("Cell types in Control_cellchat:")
print(control_cell_types)

# Check if the cell types are the same
if (setequal(MAifc_cell_types, control_cell_types)) {
  print("Both CellChat objects have the same cell composition.")
} else {
  print("The CellChat objects have different cell compositions.")
}


#########The code below is for two datasets with the same cell composition

#Compare the total number of interactions and interaction strength
ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

##Compare the number of interactions and interaction strength among different cell populations
#Form 1: circle plot
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#Form2: heatmap
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#Form3: circle plot (for the number of interactions and interaction weight)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# Get the maximum weight for scaling the interaction weights
weight.max <- getMaxWeight(object.list, attribute = c("idents", "weight"))

# Set up the plotting area
par(mfrow = c(1, 2), xpd = TRUE)

# Loop through each CellChat object to visualize
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$weight,  # Use 'weight' instead of 'count'
    weight.scale = TRUE, 
    label.edge = FALSE, 
    edge.weight.max = weight.max[2], 
    edge.width.max = 12, 
    title.name = paste0("Interaction weights - ", names(object.list)[i])
  )
}


# Form 4: Circle plot (only shown the interactions between certain cell types)

#Form4: circle plot (only shown the interactions between certain cell types)
# Define the cell types of interest
group.cellType <- c(rep("HSPCs", 1), rep("Neutrophils", 1), rep("Macrophages Monocytes", 1))

# Set the factor levels to match the cell types you're interested in
group.cellType <- factor(group.cellType, levels = c("HSPCs", "Neutrophils", "Macrophages Monocytes"))

# Merge interactions for the specified cell types
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})

# Merge the CellChat objects
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Get the maximum weight for scaling the interaction weights in the merged data
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents", "count", "count.merged"))

# Set up the plotting area
par(mfrow = c(1, 2), xpd = TRUE)

# Loop through each CellChat object to visualize
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$count.merged, 
    weight.scale = TRUE, 
    label.edge = TRUE, 
    edge.weight.max = weight.max[3], 
    edge.width.max = 12, 
    title.name = paste0("Number of interactions - ", names(object.list)[i])
  )
}

# Calculate max interaction strength for EACH DATASET
weight.max.strength <- sapply(object.list, function(x) {
  max(x@net$weight.merged)
})

# Loop through each dataset in object.list
for (i in seq_along(object.list)) {
  # Plot interaction strength for the current dataset
  netVisual_circle(
    object.list[[i]]@net$weight.merged, 
    weight.scale = TRUE, 
    label.edge = TRUE, 
    edge.weight.max = weight.max.strength[i], # Use max strength for current dataset
    edge.width.max = 12, 
    title.name = paste0("Interaction strength - ", names(object.list)[i])
  )
}

##Compare the major sources and targets in a 2D space

#1. identify cell populations with significant changes in sending or receiving signals
#scatter plot
# Compute the network centrality scores for each object in the list
object.list <- lapply(object.list, netAnalysis_computeCentrality)

# Proceed with the scatter plot
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # Control the dot size in the different datasets
gg <- list()

for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

# Plot the results
gridExtra::grid.arrange(grobs = gg, ncol = 2)


#2. Identify the signalling changes of specific cell populations
##(the code below does not work...)
# Compute centrality scores for the CellChat object
cellchat <- netAnalysis_computeCentrality(cellchat)

# Check Neutrophils
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Neutrophils", signaling.exclude = "MIF")

# Check Stem cells
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "HSPCs", signaling.exclude = "MIF")

# Combine the plots
patchwork::wrap_plots(plots = list(gg1, gg2))


#3. Identify altered signaling with distinct interaction strength
##3.1Compare the overall information flow of each signaling pathway or ligand-receptor pair
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

##3.2 Compare outgoing (or incoming) signaling patterns associated with each cell population

library(ComplexHeatmap)
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#4.Identify the up-gulated and down-regulated signaling ligand-receptor pairs

##4.1.Identify dysfunctional signaling by comparing the communication probabities
#> Comparing communications on a merged object
#> > print(MAifc_cell_types)
#
#[1] "Mature Neutrophils"                    
#[2] "Activated Neutrophils"                 
#[3] "Hematopoietic Stem/Progenitor Cells"   
#[4] "Immature Neutrophils"                  
#[5] "pre B Cells"                           
#[6] "Megakaryocytes/Platelets"              
#[7] "Dendritic Cells"                       
#[8] "Fibroblasts or Mesenchymal Cells"      
#[9] "Erythroid Cells/Erythrocyte Precursors"
#[10] "Pro-B Cells"                           
#[11] "Macrophages Monocytes"                 
#[12] "Antimicrobial Neutrophils"             
#[13] "Undefined"                             
#[14] "Mast Cells/Progenitors"                
#[15] "Natural Killer Cells"                  
#[16] "Plasma Cells"       

# activated Neutrophils as senders 
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use =2, targets.use = c(1:8),  comparison = c(1, 2), angle.x = 45)



# Stem cells as senders 
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:8),  comparison = c(1, 2), angle.x = 45)

# Monocytes as senders 
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:8),  comparison = c(1, 2), angle.x = 45)



##More specifically for the upregulated and the downregulated 
gg1 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:8),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in MAifc group", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:8),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in MAifc group", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2


###Visually compare cell-cell communication with three different plots
##Compare the signalling gene expression distribution between different datasets 
#Violin plot
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("MAifc", "Control")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T, type = "violin")

#Circle plots
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#heatmap
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}

ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

#chord diagram 
# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}


save.image(file = "scRNAseq#3.RData")
