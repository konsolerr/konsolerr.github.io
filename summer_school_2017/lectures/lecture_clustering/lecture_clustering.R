# Lecture "Clustering and dimensionality reduction". R script

library(ggplot2)
library(Rtsne)
library(dbscan)
library(pheatmap)
library(amap)


data <- read.csv("GSE89225_illumina_counts_preprocessed.csv", row.names=1)
conditions <- read.csv("conditions.csv", row.names=1)

## data <- data[, !colnames(data) %in% c("treg_NBP_patient3")]
## conditions <- conditions[!rownames(conditions) %in% c("treg_NBP_patient3"), ]

head(data)
head(conditions)


# Hclust euclidean distance
dists <- dist(t(data))
plot(hclust(dists, method="average"))
# plot(hclust(dists, method="ward.D2"))

# Hclust correlation
cors <- cor(data)
dists <- 1 - cors
dists <- as.dist(dists)
plot(hclust(dists, method="average"))

# Kmeans correlation
clustering <- Kmeans(data, 8, iter.max=20000, method="correlation")
toVisualize <- data[order(clustering$cluster), order(conditions[, 2])]
rowAnnot <- data.frame(cluster=as.factor(clustering$cluster), row.names=names(clustering$cluster))
png("heatmap_large.png", width=8, height=12, units="in", res=300)
pheatmap(toVisualize, show_rownames = F, show_colnames = F, cluster_rows = F, cluster_cols = T,
         annotation_col = conditions, scale = "row", annotation_row = rowAnnot)
dev.off()

# Performing PCA on samples

dataForPCA <- t(data)
dataForPCA <- scale(dataForPCA)

pcaData <- prcomp(dataForPCA)
percents <- pcaData$sdev^2 / sum(pcaData$sdev^2)
toPlot <- dataForPCA %*% pcaData$rotation

gdata <- data.frame(
  x=toPlot[, 1],
  y=toPlot[, 2],
  tissue=conditions[, 1],
  cells=conditions[, 2],
  name=rownames(conditions)
)

ggplot(data=gdata, aes(x=x, y=y, color=cells, shape=tissue, text=name)) +
  geom_point(size=3) + theme_bw()  +
  xlab(paste0("PC", 1, ": ", formatC(100 * percents[1], digits=4), "%")) +
  ylab(paste0("PC", 2, ": ", formatC(100 * percents[2], digits=4), "%"))

# Loading large single cell data

dataSc <- read.csv(gzfile("back_tmp/counts_2000.tsv.gz"), sep="\t", row.names=1)
dataSc <- dataSc[complete.cases(dataSc), ]
dataSc[1:5, 1:5]

# Running TSNE
tsneRes <- Rtsne(t(dataSc), check_duplicates=F, initial_dims = 10)
toPlot <- data.frame(TSNE1=tsneRes$Y[, 1], TSNE2=tsneRes$Y[, 2])
ggplot(toPlot, aes(x=TSNE1, y=TSNE2)) + geom_point() +
  theme_classic()

# Running and visualising HDBSCAN
clustering <- hdbscan(tsneRes$Y, 20)
clusters <- clustering$cluster
clusters[clusters == 0] <- NA
toPlot$cluster <- as.factor(clusters)
ggplot(toPlot, aes(x=TSNE1, y=TSNE2, color=cluster)) + geom_point() +
  theme_classic()