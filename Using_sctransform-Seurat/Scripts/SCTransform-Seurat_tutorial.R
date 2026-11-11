# Using sctransform - Seurat tutorial

# Load required packages
library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)

# packageVersion("Seurat") 

# install.packages('BiocManager')
# BiocManager::install('glmGamPoi')

# Load data and create Seurat object
pbmc_data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc_data)

# Apply sctransform normalization

# First, remove the confounding sources of variation: mitochondrial mapping percentage
# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "MT-", col.name = 'percent.mt')

# run sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)

# Dimensionality reduction (PCA) and UMAP embedding

pbmc <- RunPCA(pbmc, verbose = FALSE)
# notice: more PCs chosen
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)

p <- DimPlot(pbmc, label = TRUE)
ggsave("Dimplot.png", plot = p, width = 12, height = 8)

# Visualize canonical marker genes as violin plots
d <- VlnPlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"), pt.size = 0.2, ncol = 4)
ggsave("VlnPlots.png", plot = d, width = 12, height = 8)

# Visualize canonical marker genes on the sctransform embedding
f <- FeaturePlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7"), pt.size = 0.2, ncol = 3)
ggsave("FeaturePlot.png", plot = f, width = 12, height = 8)
g <- FeaturePlot(pbmc, features = c("CD3D", "ISG15", "TCL1A", "FCER2", "XCL1", "FCGR3A"), pt.size = 0.2, ncol = 3)
ggsave("FeaturePlot2.png", plot = g, width = 12, height = 8)

