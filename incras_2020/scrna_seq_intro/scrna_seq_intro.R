## installing and loading packages

if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("MAST", quietly = TRUE)) BiocManager::install("MAST")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")

library(Seurat)
library(MAST)
library(ggplot2)
library(dplyr)
library(Matrix)


## loading count matrix

data <- Read10X("filtered_feature_bc_matrix/")
dim(data)



## visualizing UMI distribution

plotData <- data.frame(
  umis <- colSums(data)
)
ggplot(data=plotData, aes(x=umis)) +
  geom_histogram() + theme_bw()



## Creating Seurat object 
## only keeping genes that are detected in at least 10 cells
## and keeping cells that have at least 10 genes

seurat <- CreateSeuratObject(data, min.cells = 10, min.features = 10)
dim(seurat)



## Showing UMI/gene relationship
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
FeatureScatter(seurat, "nCount_RNA", "nFeature_RNA") + scale_x_log10() + scale_y_log10()


## Showing relationship between UMIs/genes and percentage of mitochondrial RNA
FeatureScatter(seurat, "nCount_RNA", "percent.mt") + scale_x_log10()
FeatureScatter(seurat, "nFeature_RNA", "percent.mt") + scale_x_log10()


## Filtering cells with at least 300 umis and less than 25% mirochondrial reads
seurat <- subset(seurat, subset = nFeature_RNA > 300 & percent.mt < 25)
dim(seurat)


## Old way to normalize and scale the data

## seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
## seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
## seurat <- ScaleData(seurat)



## Better way to normalize and scale the data
seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", verbose = FALSE)



## Running PCA / showing standard deviation for each PCA component
seurat <- RunPCA(seurat, verbose = FALSE)
ElbowPlot(seurat, ndims = 50)


## Running TSNE dimensionality reduction / showing TSNE plot 
seurat <- RunTSNE(seurat, dims=1:20)
DimPlot(seurat, reduction = "tsne") + NoLegend()


## Running UMAP dimensionality reduction / showing UMAP plot
seurat <- RunUMAP(seurat, dims=1:20)
DimPlot(seurat, reduction = "umap") + NoLegend()


## Comparing two plots
DimPlot(seurat, reduction = "tsne") + NoLegend()
DimPlot(seurat, reduction = "umap") + NoLegend()


## Finding neighbots and performing clustering
## Also showing clusters on top of TSNE/UMAP plots
seurat <- FindNeighbors(seurat, dims = 1:20, verbose = FALSE)
seurat <- FindClusters(seurat, resolution=0.6, verbose = FALSE)
DimPlot(seurat, reduction = "tsne", label = TRUE) + NoLegend()
DimPlot(seurat, reduction = "umap", label = TRUE) + NoLegend()


## Showing expression of canonical makreks on top of UMAP plot
FeaturePlot(seurat, c("CD14", "CD79A", "CD3D"), cols=c("grey", "red"), reduction="umap", ncol=3)


## Performing DE for each cluster vs all other clusters
## max cells per ident is only set to speed up the whole thing
allMarkers <- FindAllMarkers(seurat, max.cells.per.ident = 100, test.use = "MAST", only.pos = T)
head(allMarkers)

## choosing the best marker for each cluster by log fold change
goodMarkers <- allMarkers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC) %>% pull(gene)
goodMarkers


## showing expression of these matkers on top of TSNE
FeaturePlot(seurat, goodMarkers[1:3], cols=c("grey", "red"), reduction="umap", ncol=3)


## shwoing violint plot
VlnPlot(seurat, goodMarkers[1:3], pt.size = 0.1)

