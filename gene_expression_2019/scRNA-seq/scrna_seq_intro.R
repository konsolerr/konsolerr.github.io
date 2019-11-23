## ----message=F, warning=F------------------------------------------------
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("MAST", quietly = TRUE)) BiocManager::install("MAST")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")

library(Seurat)
library(Matrix)
library(MAST)
library(ggplot2)
library(dplyr)


## ------------------------------------------------------------------------

data <- Read10X("filtered_feature_bc_matrix/")
dim(data)



## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F-----------------

plotData <- data.frame(
  umis <- colSums(data)
)
ggplot(data=plotData, aes(x=umis)) +
  geom_histogram() + theme_bw()



## ------------------------------------------------------------------------

seurat <- CreateSeuratObject(data, min.cells = 10, min.features = 10)
dim(seurat)



## ----fig.height=3, fig.width=4.5, dev='png', message=F, fig.show='hold', dpi=100----
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
FeatureScatter(seurat, "nCount_RNA", "nFeature_RNA") + scale_x_log10() + scale_y_log10()


## ----fig.height=3, fig.width=4.5, dev='png', message=F, fig.show='hold', dpi=100----
FeatureScatter(seurat, "nCount_RNA", "percent.mt") + scale_x_log10()
FeatureScatter(seurat, "nFeature_RNA", "percent.mt") + scale_x_log10()


## ----message=F, warning=F------------------------------------------------
seurat <- subset(seurat, subset = nFeature_RNA > 300 & percent.mt < 25)
dim(seurat)


## ------------------------------------------------------------------------

## seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
## seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
## seurat <- ScaleData(seurat)



## ----message=F, warning=F, cache=T---------------------------------------
seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", verbose = FALSE)



## ----fig.height=3, fig.width=6, dev='svg', message=F, fig.show='hold'----
seurat <- RunPCA(seurat, verbose = FALSE)
ElbowPlot(seurat, ndims = 50)


## ----fig.height=3, fig.width=4.5, dev='png', message=F, warning=F, fig.show='hold', dpi=100----
seurat <- RunTSNE(seurat, dims=1:20)
DimPlot(seurat, reduction = "tsne") + NoLegend()


## ----fig.height=3, fig.width=4, dev='png', message=F, warning=F, fig.show='hold', dpi=120----
seurat <- RunUMAP(seurat, dims=1:20)
DimPlot(seurat, reduction = "umap") + NoLegend()


## ----fig.height=3, fig.width=3, dev='png', message=F, fig.show='hold', dpi=120----
DimPlot(seurat, reduction = "tsne") + NoLegend()
DimPlot(seurat, reduction = "umap") + NoLegend()


## ----fig.height=2.6, fig.width=2.6, dev='png', message=F, warning=F, fig.show='hold', dpi=120----
seurat <- FindNeighbors(seurat, dims = 1:20, verbose = FALSE)
seurat <- FindClusters(seurat, resolution=0.6, verbose = FALSE)
DimPlot(seurat, reduction = "tsne", label = TRUE) + NoLegend()
DimPlot(seurat, reduction = "umap", label = TRUE) + NoLegend()


## ----fig.height=2.4, fig.width=8, dev='png', message=F, warning=F, fig.show='hold', dpi=120----
FeaturePlot(seurat, c("CD14", "CD79A", "CD3D"), cols=c("grey", "red"), reduction="umap", ncol=3)


## ----fig.height=2.4, fig.width=6, dev='png', message=F, warning=F, fig.show='hold', dpi=120, cache=T----
# max cells per ident is only seed to speed up the whole thing
allMarkers <- FindAllMarkers(seurat, max.cells.per.ident = 100, test.use = "MAST", only.pos = T)
goodMarkers <- allMarkers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC) %>% pull(gene)
goodMarkers


## ----fig.height=2.4, fig.width=8, dev='png', message=F, warning=F, fig.show='hold', dpi=120----
FeaturePlot(seurat, goodMarkers[1:3], cols=c("grey", "red"), reduction="umap", ncol=3)


## ----fig.height=2.4, fig.width=8, dev='png', message=F, warning=F, fig.show='hold', dpi=120----
VlnPlot(seurat, goodMarkers[1:3], pt.size = 0.1)

