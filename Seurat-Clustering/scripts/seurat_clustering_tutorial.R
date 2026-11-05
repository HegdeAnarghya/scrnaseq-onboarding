# Seurat- Guided Clustering tutorial

# Load necessary packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Dataset: PBMC (Peripheral Blood Mono-nuclear Cells)(Immune cells with single-round nucleus;
# PBMC populations consists of T cells, B cells, Natural Killer (NK) cells, Monocytes and dendritic cells)

# Load the dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Examine few genes in the first 30 cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# Check for memory requirement of the PBMC gene expression data in dense and sparse format
# dense matrix = stores 0s and non-zero values; sparse matrix stores non-zero values and their positions
dense.size <- object.size(as.matrix(pbmc.data))
dense.size

sparse.size <- object.size(pbmc.data)
sparse.size

dense.size/sparse.size

# ==============================================================

# Pre-processing
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# View QC for first 5 cells
head(pbmc@meta.data,5)

# Visualizing QC metrics
# filtering unique feature counts over 2500 or less than 200;
# filtering cells with > 5% percent.mt

# Violin Plot
qc_plot <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave("PBMC_QC_Violin_Plot.png", width = 12, height = 6, dpi = 300)

# FeatureScatter
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

ggsave("PBMC_QC_Scatter_Plot.png", width = 12, height = 6, dpi = 300)

# Perform quality control filtering on the PBMC dataset
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# ===============================================================

# Normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


# Feature Selection (Identifying highly variable features)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

ggsave("PBMC_VariableFeaturesPlot.png", width = 12, height =6, dpi = 300)

# ===============================================================

# Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Remove unwanted sources of variation
# pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

# ===============================================================

# Linear Dimensional reduction

pbmc <-RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results in different ways

# 1
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# 2
p <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
ggsave("VizDimLoadings_PC1_PC2.png", plot = p, width = 12, height =6, dpi = 300)

# 3
d <- DimPlot(pbmc, reduction = "pca") + NoLegend()
ggsave("DimPlot_PC1_PC2.png", plot = d, width = 12, height =6, dpi = 300)

# 4 
b <- DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE, fast = FALSE)
ggsave("DimHeatmap_PC1.png", plot = b, width = 12, height =6, dpi = 300)

a <- DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)
ggsave("DimHeatmap_PC15.png", plot = a, width = 34, height =28, dpi = 300)

# =================================================================

# Determining the dimensionality of the dataset

k <- ElbowPlot(pbmc)
ggsave("ElbowPlot.png", plot = k, width = 12, height =6, dpi = 300)

# =================================================================

# Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# view the cluster IDs of first 5 cells
head(Idents(pbmc), 5)

# =================================================================

# Run non-linear dimensional reduction (UMAP/t=SNE)
pbmc <- RunUMAP(pbmc, dims = 1:10)

P <- DimPlot(pbmc, reduction = "umap")
ggsave("UMAP.png", plot = P, width = 12, height =6, dpi = 300)

# save the work at this point
saveRDS(pbmc, file = "../pbmc_tutorial.rds")

# =================================================================

# Determining the differentially expressed features (cluster biomarkers)

# Find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n=5)

# Finding all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n=5)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# Visualization

v <- VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
ggsave("VlnPlot_clusters.png", plot = v, width = 12, height =6, dpi = 300)

# Plotting raw counts
V <- VlnPlot(pbmc, features = c("NKG7", "PF4"), layer = "counts", log = TRUE)
ggsave("VlnPlot_counts.png", plot = V, width = 12, height =6, dpi = 300)

F <- FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))
ggsave("FeaturePlot.png", plot = F, width = 12, height =6, dpi = 300)

# Generating an expression heatmap for given cells and features
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC >1) %>%
  slice_head(n=10) %>%
  ungroup() -> top10
D <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()
ggsave("DoHeatmap.png", plot = D, width = 12, height =6, dpi = 300)

# ===================================================================

# Assigning cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
Pr <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("Cell_type_assigned-umap.png", plot = Pr, width = 12, height =6, dpi = 300)

# plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
#  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
# ggsave("Cell_type_assigned-umap1.png", plot = plot, width = 12, height =6, dpi = 300)

saveRDS(pbmc, file = "../pbmc3k_final.rds")
