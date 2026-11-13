# Using Seurat with multimodal data

# Load Required Libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)

# Loading the RNA UMI matric
# Considering that this dataset also contains ~5% of mouse cells, which can be used as negative controls for the protein measurements
# For this reason, the gene expression matrix has HUMAN_ or MOUSE_ appended to the beginning of each gene
cbmc.rna <- as.sparse(read.csv(file = "Data/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", header = TRUE, row.names = 1))
# Discard all but the top 100 most highly expressed mouse genes, and remove the "HUMAN_" from the CITE-seq prefix
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)

# Loading the ADT UMI matrix
cbmc.adt <- as.sparse(read.csv(file = "Data/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", header = TRUE, row.names = 1))

# check if the two matrices have identical column names, since measurements were made in the same cells
all.equal(colnames(cbmc.rna), colnames(cbmc.adt))

# Set up a Seurat object with the RNA and ADT data (Seurat object with RNA data and add ADT as a second assay)
cbmc <- CreateSeuratObject(counts = cbmc.rna)
Assays(cbmc)

# Creating new assay to tore the ADT information
adt_assay <- CreateAssay5Object(counts = cbmc.adt)

# Add this assay to previously created seurat object
cbmc[["ADT"]] <- adt_assay

# Validate that the object now contains multiple assays
Assays(cbmc)

# Extract a list of features measured in the ADT assay
rownames(cbmc[["ADT"]])

# List the current default assay
DefaultAssay(cbmc)

# switch default assay to ADT
DefaultAssay(cbmc) <- 'ADT'
DefaultAssay(cbmc)

# setting it back to RNA
DefaultAssay(cbmc) <- 'RNA'

# ==========================================================

# Cluster cells on the basis of their scRNA-seq profiles

# Performing visualization and clustering steps
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = FALSE)
cbmc <- FindNeighbors(cbmc, dims = 1:30)
cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:30)

d <- DimPlot(cbmc, label = TRUE)
ggsave("scRNA-seq-cluster-Dimplot.png", plot = d, width = 12, height = 8)

# ==========================================================

# Visualizing multiple modalities side by side

# Normalize ADT data
DefaultAssay(cbmc) <- "ADT"
cbmc <- NormalizeData(cbmc, normalization.method = 'CLR', margin = 2)
DefaultAssay(cbmc) <- 'RNA'

# Alternative
cbmc <- NormalizeData(cbmc, normalization.method = 'CLR', margin = 2, assay = 'ADT')

# Visualize CD19 levels for RNA and protein
DefaultAssay(cbmc) <- 'ADT'
p1 <- FeaturePlot(cbmc, "CD19", cols = c("lightgrey","darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(cbmc) <- 'RNA'
p2 <- FeaturePlot(cbmc, "CD19") + ggtitle("CD19 RNA")

# Place plots side-by-side
p <- p1 | p2
ggsave("CD19-protein-RNA-UMAP.png", plot = p, width=12, height=8)

# Alternatively, using specific assay keys to specify a specific modality
# IDentify the key for the RNA and protein assays
Key(cbmc[["RNA"]])
Key(cbmc[["ADT"]])

# Now, we can include the key in the feature name, which overrides the default assay
p1 <- FeaturePlot(cbmc, "adt_CD19",cols = c("lightgrey","darkgreen")) + ggtitle("CD19 protein")
p2 <- FeaturePlot(cbmc, "rna_CD19") + ggtitle("CD19 RNA")
p1 | p2

# ==========================================================

# Identify cell surface markers for scRNA-seq clusters

# As we know that CD19 is a B cell marker, we can identify cluster 6 as expressing CD19 on the surface
pr <- VlnPlot(cbmc,"adt_CD19")
ggsave("adt_CD19-VlnPlot.png", plot = pr, width = 12, height = 8)

# We can also identify alternative protein and RNA markers for this cluster through differential expression
adt_markers <- FindMarkers(cbmc, ident.1 = 6, assay = "ADT")
rna_markers <- FindMarkers(cbmc, ident.1 = 6, assay = "RNA")

head(adt_markers)
head(rna_markers)

# Additional visualizations of multimodal data

# ADT scatter plots
f <- FeatureScatter(cbmc, feature1 = 'adt_CD19', feature2 = 'adt_CD3')
ggsave("featureScatterPlot-adtCD3-adtCD19.png", plot = f, width = 12, height = 8)

# view relationship between protein and RNA
fr <- FeatureScatter(cbmc, feature1 = 'adt_CD3', feature2 = 'rna_CD3E')
ggsave("featureScatterPlot-adtCD3-rnaCD3.png", plot = fr, width = 12, height = 8)

d <- FeatureScatter(cbmc, feature1 = 'adt_CD4', feature2 = 'adt_CD8')
ggsave("featureScatterPlot-adtCD4-adtCD8.png", plot = d, width = 12, height = 8)


#looking at the raw (non-normalized) ADT counts. The values are quite high, particularly in comparison to RNA values
# This is due to the significantly higher protein copy number in cells, which significantly reduces 'drop-out' in ADT data  
c <- FeatureScatter(cbmc, feature1 = 'adt_CD4', feature2 = 'adt_CD8', slot = 'counts')
ggsave("featureScatterPlot-adtCD4-adtCD8-rawADTCounts.png", plot = c, width = 12, height = 8)


## ========================================================

# Loading data from 10X multi-modal experiments

# Recreating the above plots using the 10x Genomics Flex 2, 16k Human PBMCs, TotalSeq-C 4-plex dataset
pbmc10k.data <- Read10X(data.dir = 'Data/sample_filtered_feature_bc_matrix/')

# Clean the names of the proteins (features) in the Antibody Capture matrix by removing unwanted text from their names
rownames(x = pbmc10k.data[['Antibody Capture']]) <- gsub(
  pattern = '_[control_]*TotalSeqB',
  replacement = '',
  x = rownames(x = pbmc10k.data[['Antibody Capture']])
)

names(pbmc10k.data)
head(rownames(pbmc10k.data[['Antibody Capture']]))

# Create Seurat object
pbmc10K <- CreateSeuratObject(counts = pbmc10k.data[['Gene Expression']], min.cells = 3, min.features = 200)
pbmc10k <- NormalizeData(pbmc10K)
pbmc10k[['ADT']] <- CreateAssay5Object(pbmc10k.data[['Antibody Capture']][, colnames(x=pbmc10k)])
pbmc10k <- NormalizeData(pbmc10k, assay = 'ADT', normalization.method = 'CLR')

# check (Different ADT features' names from those in the tutorial)
head(rownames(pbmc10k[["ADT"]]),30)
Key(pbmc10k[["ADT"]])

grep("CD19|CD3|CD4|CD8",
     rownames(pbmc10k[["ADT"]]),
     value = TRUE)

# visualization
plot1 <- FeatureScatter(pbmc10k, feature1 = 'adt_Hu.CD19', feature2 = 'adt_Hu.CD3-UCHT1', pt.size = 1)
plot2 <- FeatureScatter(pbmc10k, feature1 = 'adt_Hu.CD4-RPA.T4', feature2 = 'adt_Hu.CD8', pt.size = 1)
plot3 <- FeatureScatter(pbmc10k, feature1 = 'adt_Hu.CD3-UCHT1', feature2 = 'CD3E', pt.size = 1)
p <- (plot1 + plot2 + plot3) & NoLegend()
ggsave('FeatureScatterPlots-pbmc10k-1-2-3.png', plot = p, width = 12, height = 8)


pr <- FeatureScatter(
  pbmc10k,
  feature1 = "adt_Hu.CD19",
  feature2 = "adt_Hu.CD3-UCHT1"
) +
  NoLegend() + 
  theme (
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )

print(pr)

# Check for clusters with higher CD19 or CD3 expression

# Cluster cells in the pbmc10k dataset using RNA data first
DefaultAssay(pbmc10k) <- "RNA"

pbmc10k <- NormalizeData(pbmc10k)

pbmc10k <- FindVariableFeatures(pbmc10k)

pbmc10k <- ScaleData(pbmc10k)

pbmc10k <- RunPCA(pbmc10k, verbose = FALSE)

pbmc10k <- FindNeighbors(pbmc10k, dims = 1:30)

pbmc10k <- FindClusters(
  pbmc10k,
  resolution = 0.8,
  verbose = FALSE
)

pbmc10k <- RunUMAP(pbmc10k, dims = 1:30)

# Visualize the clusters
d <- DimPlot(pbmc10k, reduction = "umap", label = TRUE)
ggsave("pbmc10k-Dimplot.png", plot = d, width = 12, height = 8)  

v <- VlnPlot(
  pbmc10k,
  features = c("Hu.CD19", "Hu.CD3-UCHT1"),
  assay = "ADT",
  group.by = "seurat_clusters",
  pt.size = 0.1
)
ggsave("pbmc10k-Vlnplot.png", plot = v, width = 12, height = 8)  
