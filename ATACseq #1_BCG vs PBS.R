if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DiffBind", "ChIPseeker", "GenomicRanges", "rtracklayer", "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"))

library(DiffBind)
library(ChIPseeker)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

# Define your samples directly in R [always put control first in the c()!]# Save everything in your current workspace to a file named "my_workspace.RData"
samples <- data.frame(
  SampleID = c("PBS1", "PBS2", "PBS3", "BCG1", "BCG2", "BCG3"),
  Condition = c("Control", "Control", "Control","Treatment", "Treatment", "Treatment"),
  Replicate = c(1, 2, 3, 1, 2, 3),
  bamReads = c("scRNAseq/Samtools_sort_PBS1.bam", "scRNAseq/Samtools_sort_PBS2.bam", "scRNAseq/Samtools_sort_PBS3.bam", "scRNAseq/Samtools_sort_BCG1.bam", "scRNAseq/Samtools_sort_BCG2.bam", "scRNAseq/Samtools_sort_BCG3.bam"),
  Peaks = c("scRNAseq/MACS2_callpeak_PBS1.bed", "scRNAseq/MACS2_callpeak_PBS2.bed", "scRNAseq/MACS2_callpeak_PBS3.bed","scRNAseq/MACS2_callpeak_BCG1.bed", "scRNAseq/MACS2_callpeak_BCG2.bed", "scRNAseq/MACS2_callpeak_BCG3.bed"),
  PeakCaller = c("bed", "bed", "bed", "bed", "bed", "bed")
)

# Load data into DiffBind
dba <- dba(sampleSheet = samples)

#Count Reads
dba <- dba.count(dba, summits = 250)

#Differential analysis 
dba <- dba.contrast(dba, categories = DBA_CONDITION, minMembers = 2)
dba <- dba.analyze(dba)
dba.report <- dba.report(dba)

# View the first few rows
head(dba.report)

#Annotate peaks 
peakAnno <- annotatePeak(dba.report, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, 
                         annoDb = "org.Mm.eg.db")


# Assuming peakAnno is your object
peak_anno_df <- as.data.frame(peakAnno)

# View the first few rows
head(peak_anno_df)


#### FOR VISUALIZATION 

library(dplyr)
library(ggplot2)

# Assuming peak_anno_df has columns: Fold, FDR, annotation, SYMBOL
# Select top 10 most significant peaks based on FDR
top_peaks <- peak_anno_df %>%
  arrange(FDR) %>%
  head(10)

ggplot(peak_anno_df, aes(x = Fold, y = -log10(FDR))) +
  geom_point(aes(color = Fold > 0)) +
  scale_color_manual(values = c("TRUE" = "pink", "FALSE" = "lightblue")) +
  geom_text(aes(label = paste(annotation, SYMBOL, sep = " - ")), vjust = 1.5, check_overlap = TRUE,size=2) +
  theme_classic() +
  labs(x = "Fold Change", y = "-log10 FDR", color = "Direction")

ggplot(peak_anno_df, aes(x = Fold, y = -log10(FDR))) +
  geom_point(aes(color = Fold > 0)) +
  scale_color_manual(values = c("TRUE" = "pink", "FALSE" = "lightblue")) +
  geom_text(aes(label = paste(SYMBOL, sep = " - ")), vjust = 1.5, check_overlap = TRUE,size=2) +
  theme_classic() +
  labs(x = "Fold Change", y = "-log10 FDR", color = "Direction")


# Peak annotation plot (pie)
plotAnnoPie(peakAnno)

# Plot heatmap 
dba.plotHeatmap(dba)

##plot heatmap of the peak_anno_df (does not work)

library(pheatmap)

# Example: Create a count matrix (substitute with your actual data)
# Assuming peak_anno_df has columns: SYMBOL, and some numerical data for heatmap
# Select relevant columns and set SYMBOL as row names
data_matrix <- as.matrix(peak_anno_df[, c("Fold", "-log10(FDR)")])
rownames(data_matrix) <- peak_anno_df$SYMBOL

# Plot the heatmap
pheatmap(data_matrix, 
         main = "Heatmap of Fold Change and FDR",
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         scale = "row")

##plot GO 
# Install packages if needed
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Mm.eg.db")

library(clusterProfiler)
library(org.Mm.eg.db)

# Extract unique gene symbols
genes <- unique(peak_anno_df$SYMBOL)

# Run GO enrichment analysis
go_results <- enrichGO(gene         = genes,
                       OrgDb        = org.Mm.eg.db,
                       keyType      = "SYMBOL",
                       ont          = "BP", # Biological Process
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

# View results
summary(go_results)

# Plot the GO enrichment results
dotplot(go_results, showCategory = 10) + 
  ggtitle("GO Enrichment Analysis")



# Load necessary libraries
library(ggplot2)
library(pheatmap)

## Peak annotation plots

# 1. Bar plot of peak annotations
barplot(table(peakAnno$Annotation), main = "Peak Annotations", xlab = "Annotation", ylab = "Frequency")

# 2. Pie chart of peak annotations
pie(table(peakAnno$Annotation), main = "Peak Annotations", col = rainbow(length(unique(peakAnno$Annotation))))

# 3. Heatmap of peak annotations vs. genomic features
pheatmap(table(peakAnno$Annotation, peakAnno$GenomicFeature), main = "Peak Annotations vs. Genomic Features", cluster_rows = FALSE, cluster_cols = FALSE)

# DBA report plots

# Volcano plot (already discussed)
ggplot(dba.report, aes(x = Fold, y = -log10(FDR), color = Significant)) +
  geom_point() +
  scale_color_manual(values = c("No" = "gray", "Yes" = "red")) +
  theme_classic() +
  labs(x = "Fold Change", y = "-log10 FDR", color = "Significant")

# MA plot
plot(dba.report$Fold, dba.report$-log10(FDR), main = "MA Plot", xlab = "Fold Change", ylab = "-log10 FDR")

# Histogram of fold changes
hist(dba.report$Fold, breaks = 50, main = "Histogram of Fold Changes", xlab = "Fold Change", ylab = "Frequency")

# Histogram of FDR values
hist(dba.report$FDR, breaks = 50, main = "Histogram of FDR Values", xlab = "FDR", ylab = "Frequency")

# DBA plots

# Correlation heatmap of DBA samples
pheatmap(cor(dba$counts), main = "Correlation Heatmap of DBA Samples", cluster_rows = FALSE, cluster_cols = FALSE)

# Principal Component Analysis (PCA) plot of DBA samples
plot(prcomp(dba$counts)$x[, 1:2], main = "PCA Plot of DBA Samples", xlab = "PC1", ylab = "PC2")

# Hierarchical clustering of DBA samples
plot(hclust(dist(dba$counts)), main = "Hierarchical Clustering of DBA Samples", xlab = "", ylab = "")

# Boxplot of DBA sample counts
boxplot(dba$counts, main = "Boxplot of DBA Sample Counts", xlab = "Sample", ylab = "Count")

### Peak intensity plots
# Access the 'anno' slot and convert to data frame
peak_data <- as.data.frame(peakAnno@anno)

dev.new() 

# Boxplot of peak intensities
boxplot(peak_data$Conc, main = "Boxplot of Peak Intensities", xlab = "", ylab = "Intensity")

dev.new() 

# Histogram of peak intensities
hist(peak_data$Conc, breaks = 50, main = "Histogram of Peak Intensities", xlab = "Intensity", ylab = "Frequency")

dev.new() 

# Adjust plot margins
par(mar = c(5, 5, 4, 4))

# Scatterplot of peak intensities vs. fold changes (+colour)
plot(peak_data$Conc, peak_data$Fold, main = "Scatterplot of Peak Intensities vs. Fold Changes", xlab = "Intensity", ylab = "Fold Change")


# Genomic feature plots

# Bar plot of genomic feature frequencies
barplot(table(peakAnno$GenomicFeature), main = "Genomic Feature Frequencies", xlab = "Genomic Feature", ylab = "Frequency")

# Pie chart of genomic feature frequencies
pie(table(peakAnno$GenomicFeature), main = "Genomic Feature Frequencies", col = rainbow(length(unique(peakAnno$GenomicFeature))))

# Heatmap of genomic feature frequencies vs. peak annotations
pheatmap(table(peakAnno$GenomicFeature, peakAnno$Annotation), main = "Genomic Feature Frequencies vs. Peak Annotations", cluster_rows = FALSE, cluster_cols = FALSE)

# Heatmap of binding affinity
dba.plotHeatmap(dba)


# Save everything in your current workspace to a file named "my_workspace.RData"
save.image(file = "BCG vs PBS_BMDM_ATAC.RData")


