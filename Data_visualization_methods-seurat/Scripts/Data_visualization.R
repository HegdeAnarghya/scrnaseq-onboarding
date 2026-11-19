# Data visualization methods in Seurat

# Note: Using the data - 2700 PBMC from guided clustering tutorial.

# Installing the required packages and dataset:
install.packages("remotes")
remotes::install_github("satijalab/seurat-data")

library(SeuratData)
SeuratData::InstallData("pbmc3k")

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)

pbmc3k.final <- LoadData("pbmc3k", type = "pbmc3k.final")
pbmc3k.final[["RNA"]] <- as(pbmc3k.final[["RNA"]], Class = "Assay5")
pbmc3k.final <- NormalizeData(pbmc3k.final)
pbmc3k.final <- FindVariableFeatures(pbmc3k.final)
pbmc3k.final <- ScaleData(pbmc3k.final)
pbmc3k.final$groups <- sample(c("group1", "group2"), size = ncol(pbmc3k.final), replace = TRUE)
features <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")
pbmc3k.final

# Visualizations of marker feature expression

# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
r <- RidgePlot(pbmc3k.final, features = features, ncol = 2)
ggsave("RidgePlot.png", plot = r, width = 12, height = 8)

# Violin plot - Visualize single cell expression distributions in each cluster
v <- VlnPlot(pbmc3k.final, features = features)
ggsave("VlnPlot.png", plot = v, width = 12, height = 8)

# Feature plot - visualize feature expression in low-dimensional space
f <- FeaturePlot(pbmc3k.final, features = features)
ggsave("FeaturePlot.png", plot = f, width = 12, height = 8)

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the feature in each cluster.
# The color represents the average expression level

d <- DotPlot(pbmc3k.final, features = features) + RotatedAxis()
ggsave("Dotplot.png", plot = d, width = 12, height = 8)

# Single cell heatmap of feature expression
dh <- DoHeatmap(subset(pbmc3k.final, downsample = 100), features = features, size = 3)
ggsave("DoHeatmap.png", plot = dh, width = 12, height = 8)

# Customizing FeaturePlot
# Plot a legend to map colors to expression levels
fe <- FeaturePlot(pbmc3k.final, features = 'MS4A1')
ggsave("Customized-featurePlot.png", plot = fe, width = 12, height = 8)

# Adjust the contrast in the plot
fe1 <- FeaturePlot(pbmc3k.final, features = 'MS4A1', min.cutoff = 1, max.cutoff = 3)
ggsave("Customized-featurePlot1.png", plot = fe1, width = 12, height = 8)

# Calculate feature-specific contrast levels based on quantiles of non-zero expression
fe2 <- FeaturePlot(pbmc3k.final, features = c('MS4A1', "PTPRCAP"), min.cutoff = "q10", max.cutoff = "q90")
ggsave("Feature-specific-contrast.png", plot = fe2, width = 12, height = 8)

# Visualize co-expression of two features simultaneously
fe3 <- FeaturePlot(pbmc3k.final, features = c('MS4A1', 'CD79A'), blend = TRUE)
ggsave("Features-co-expression.png", plot = fe3, width = 12, height = 8)

# Split visualization to view expression by groups (replaces FeatureHeatmap)
fe4 <- FeaturePlot(pbmc3k.final, features = c('MS4A1', 'CD79A'), split.by = 'groups')
ggsave("groups-expressions.png", plot = fe4, width = 12, height = 8)

# Other visualization functionality
# Violin plots can also be split on some variable. Simply add the splitting variable to object metadata and pass it to the 'split.by' argument
v1 <- VlnPlot(pbmc3k.final, features = 'percent.mt', split.by = 'groups')
ggsave("vlnplot-split.png", plot = v1, width = 12, height = 8)

# SplitDotPlotGG is replaced with the `split.by` parameter for DotPlot
d1 <- DotPlot(pbmc3k.final, features = features, split.by = 'groups') + RotatedAxis()
ggsave("DotPlot-split.png", plot = d1, width = 12, height = 8)

# DimPlot replaces TSNEPlot, PCAPlot, etc., and plots either "umap", "tsne", or "pca" by default, in that order
D <- DimPlot(pbmc3k.final)
ggsave("Dimplot.png", plot = D, width = 12, height = 8)

pbmc3k.final.no.umap <- pbmc3k.final
pbmc3k.final.no.umap[['umap']] <- NULL
D1 <- DimPlot(pbmc3k.final.no.umap) + RotatedAxis()
ggsave("Dimplot-no-umap.png", plot = D1, width = 12, height = 8)

# DoHeatmap now shows a grouping bar, splitting the heatmap into groups or clusters. This can be changed with the `group.by` parameter
dh <- DoHeatmap(pbmc3k.final, features = VariableFeatures(pbmc3k.final)[1:100], cells = 1:500, size = 4, angle = 90) + NoLegend()
ggsave("DoHeatmap-group.png", plot = dh, width = 11, height = 14)

# Applying themes to plots:
baseplot <- DimPlot(pbmc3k.final, reduction = 'umap')
# Add custom labels and titles
baseplot + labs(title = 'Clustering of 2,700 PBMCs')

bp <- baseplot + DarkTheme()
ggsave("baseplot.png", plot = bp, width = 11, height = 14)

bp1 <- baseplot + FontSize(x.title = 20, y.title = 20) + NoLegend()
ggsave("baseplot1.png", plot = bp1, width = 11, height = 14)

plot1 <- baseplot + DarkTheme() + 
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size = 10)))

ggsave("baseplot2.png", plot = plot1, width = 12, height = 8)

# Interactive plotting:
# Include additional data to display alongside cell names by passing in a data frame of information.
# Works well when using FetchData
plot2 <- FeaturePlot(pbmc3k.final, features = 'MS4A1')
pt2 <- HoverLocator(plot = plot2, information = FetchData(pbmc3k.final, vars = c('ident', 'PC_1', 'nFeature_RNA')))

# saving as an interactive html file
htmlwidgets::saveWidget(
  widget = pl2,
  file = "interactive-featureplot.html"
)

# Manually select cells for further investigation
# assuming that DCs had merged with monocytes in the clustering, but we wanted to see what was unique about them based on their position in the tSNE plot
pbmc3k.final <- RenameIdents(pbmc3k.final, "DC" = "CD14+ Mono")
plot <- DimPlot(pbmc3k.final, reduction = 'umap')
select.cells <- CellSelector(plot = plot)

# Change the identity of these cells to turn them into their own mini-cluster
head(select.cells)

Idents(pbmc3k.final, cells = select.cells) <- 'NewCells'

# find markers that are specific to the new cells, and find clear DC markers
newcells.markers <- FindMarkers(pbmc3k.final, ident.1 = "NewCells", ident.2 = "CD14+ Mono", min.diff.pct = 0.3, only.pos = TRUE)
head(newcells.markers)

# Using cellSelector to automatically assign cell identities
pbmc3k.final <- CellSelector(plot = plot, object = pbmc3k.final, ident = 'selected')
levels(pbmc3k.final)

# Plotting accessories (Functions to manipulate and combine plots)
# LabelClusters and LabelPoints will label clusters (a coloring variable) or individual points on a ggplot2-based scatter plot
plot3 <- DimPlot(pbmc3k.final, reduction = 'pca') + NoLegend()
pt3 <- LabelClusters(plot = plot3, id = 'ident')
ggsave("Labelled-Clusters.png", plot = pt3, width = 9, height = 5)

# Both functions support 'repel', which will intelligently stagger labels and draw connecting lines from the labels to the points or clusters
pt3a <- LabelPoints(plot = plot, points = TopCells(object = pbmc3k.final[['pca']]), repel = TRUE)
ggsave("Labelled-points-cluster.png", plot = pt3a, width = 12, height = 8)

# Combining two plots
plot4 <- DimPlot(pbmc3k.final)
# Create scatter plot with the Pearson correlation value as the title
plot5 <- FeatureScatter(pbmc3k.final, feature1 = 'LYZ', feature2 = 'CCL5')
# combine two plots
pt4 <- plot4 + plot5
ggsave("CombinePlots.png", plot = pt4, width = 12, height = 8)
# removing the legend from all plots
pt4a <- (plot4 + plot5) & NoLegend()
ggsave("CombinePlots-NoLegend.png", plot = pt4a, width = 12, height = 8)
