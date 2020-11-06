## ----message=F, warning=F------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("apeglm", quietly = TRUE)) BiocManager::install("apeglm")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("vsn", quietly = TRUE)) BiocManager::install("vsn")


## ----message=F, warning=F------------------------------------------------
library(DESeq2)
library(ggplot2)
library(vsn)
library(pheatmap)
library(ggrepel)
library(dplyr)
set.seed(1)

files <- list.files("WT_countdata", full.names = T)
sampleNames <- gsub(".*/(WT_rep\\d+_MID\\d+).*", "\\1", files)

counts <- lapply(files, function(file){
  read.table(file, sep="\t", row.names = 1, header = 0)
})


## ----message=F, warning=F------------------------------------------------
counts <- do.call(cbind, counts)
colnames(counts) <- sampleNames
counts[c(1:6, (nrow(counts) - 6):nrow(counts)), 1:3]


## ------------------------------------------------------------------------
counts <- counts[1:(nrow(counts) - 5), ]
dim(counts)


## ----message=F, warning=F, fig.width=6, fig.height=3, fig.fullwidth=T, dev='svg'----

ggplot(data=data.frame(librarySize=colSums(counts)), aes(x=librarySize)) +
  geom_histogram() + theme_bw() + scale_x_log10()


## ------------------------------------------------------------------------
sampleByCoverage <- colnames(counts)[order(colSums(counts))]
lowcov <- sampleByCoverage[1]
highcov <- sampleByCoverage[length(sampleByCoverage)]


## ----message=F, warning=F, fig.width=6, fig.height=3, fig.fullwidth=T, dev='svg'----
ggplot(data=data.frame(highcov=counts[, highcov], 
                       lowcov=counts[, lowcov]), 
       aes(x=lowcov, y=highcov)) +
  geom_point() + theme_bw() + scale_x_log10() + scale_y_log10()


## ----message=F, warning=F, fig.width=6, fig.height=3, fig.fullwidth=T, dev='svg'----
ggplot(data=data.frame(highcov= 1000000 * counts[, highcov] / sum(counts[, highcov]), 
                       lowcov=  1000000 * counts[, lowcov] / sum(counts[, lowcov])), 
       aes(x=lowcov, y=highcov)) +
  geom_point() + theme_bw() + scale_x_log10() + scale_y_log10()


## ----message=F, warning=F, fig.width=8, fig.height=3, fig.fullwidth=T, dev='svg'----
dds <- DESeqDataSetFromMatrix(counts, 
                              data.frame(strain=rep("WT", ncol(counts)),
                                         libsize=colSums(counts),
                                         sampleName=sampleNames,
                                         row.names = sampleNames), design=~1)
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds)


## ----message=F, warning=F, fig.width=4, fig.height=3, fig.fullwidth=T, dev='svg', fig.show="hold"----
hist(log10(colSums(counts(dds))), main = "Library size")
hist(log10(colSums(counts(dds, normalized=T))), main = "Adjusted library size")


## ---- message=F, warning=F, fig.width=4, fig.height=3, fig.fullwidth=T, dev='svg', fig.show="hold"----

ntd <- normTransform(dds)
meanSdPlot(assay(ntd), ranks = F)
meanSdPlot(assay(ntd)) + theme_bw()


## ---- message=F, warning=F, fig.width=4, fig.height=3, fig.fullwidth=T, dev='svg', fig.show="hold", cache=T----

rld <- rlogTransformation(dds)
meanSdPlot(assay(rld), ranks = F)
meanSdPlot(assay(rld)) + theme_bw()


## ---- message=F, warning=F, fig.width=4, fig.height=3, fig.fullwidth=T, dev='svg', fig.show="hold", cache=T----

vsd <- varianceStabilizingTransformation(dds)
meanSdPlot(assay(vsd), ranks = F) 
meanSdPlot(assay(vsd)) + theme_bw()


## ---- message=F, warning=F, fig.width=3, fig.height=3, fig.fullwidth=T, dev='svg', fig.show="hold"----

plotPCA(ntd, "libsize") + theme_bw()
plotPCA(rld, "libsize") + theme_bw()
plotPCA(vsd, "libsize") + theme_bw()


## ---- message=F, warning=F, fig.width=5, fig.height=4, dev='svg', fig.show="hold"----

plotPCA(vsd, "libsize") + geom_text_repel(aes(label=name), size=2) + theme_bw()


## ------------------------------------------------------------------------
filesSnf2 <- list.files("Snf2_countdata", full.names = T)
sampleNamesSnf2 <- gsub(".*/(Snf2_rep\\d+_MID\\d+).*", "\\1", filesSnf2)

countsSnf2 <- lapply(filesSnf2, function(file){
  read.table(file, sep="\t", row.names = 1, header = 0)
})

countsSnf2 <- do.call(cbind, countsSnf2)
colnames(countsSnf2) <- sampleNamesSnf2
countsSnf2 <- countsSnf2[1:(nrow(countsSnf2) - 5), ]
dim(countsSnf2)


## ----message=F, warning=F, cache=T---------------------------------------
ddsSnf2 <- DESeqDataSetFromMatrix(countsSnf2, 
                              data.frame(strain=rep("snf2_ko", ncol(countsSnf2)),
                                         libsize=colSums(countsSnf2),
                                         sampleName=sampleNamesSnf2,
                                         row.names = sampleNamesSnf2), design=~1)
ddsSnf2 <- ddsSnf2[rowSums(counts(ddsSnf2)) > 0, ]
ddsSnf2 <- DESeq(ddsSnf2)
ntdSnf2 <- normTransform(ddsSnf2)
rldSnf2 <- rlogTransformation(ddsSnf2)
vsdSnf2 <- varianceStabilizingTransformation(ddsSnf2)


## ---- message=F, warning=F, fig.width=3, fig.height=3, fig.fullwidth=T, dev='svg', fig.show="hold"----

plotPCA(ntdSnf2, "libsize") + theme_bw()
plotPCA(rldSnf2, "libsize") + theme_bw()
plotPCA(vsdSnf2, "libsize") + theme_bw()


## ---- message=F, warning=F, fig.width=5, fig.height=5, dev='svg'---------
plotPCA(vsdSnf2, "libsize") + geom_text_repel(aes(label=name), size=2) + theme_bw()


## ----message=F, warning=F------------------------------------------------
counts <- counts[, -c(21, 22, 25, 28)]
counts <- counts[, sample(ncol(counts), 5)]
countsSnf2 <- countsSnf2[, -c(6, 13, 35)]
countsSnf2 <- countsSnf2[, sample(ncol(countsSnf2), 5)]

ddsMerged <- DESeqDataSetFromMatrix(
  cbind(counts, countsSnf2), 
  data.frame(strain=c(rep("wt", ncol(counts)), rep("snf2_ko", ncol(countsSnf2))),
            libsize=c(colSums(counts), colSums(countsSnf2)),
            sampleName=c(colnames(counts), colnames(countsSnf2))), 
  design=~strain)
ddsMerged <- ddsMerged[rowSums(counts(ddsMerged)) > 0, ]
ddsMerged <- DESeq(ddsMerged)
ntdMerged <- normTransform(ddsMerged)
rldMerged <- rlogTransformation(ddsMerged)
vsdMerged <- varianceStabilizingTransformation(ddsMerged)


## ---- message=F, warning=F, fig.width=4, fig.height=4, fig.fullwidth=T, dev='svg', fig.show="hold"----
plotPCA(ntdMerged, "strain") + theme_bw()
plotPCA(vsdMerged, "strain") + theme_bw()


## ----message=F, warning=F, fig.width=8, fig.height=3, fig.fullwidth=T, dev='svg'----
plot(dnbinom(1:40, size=5, prob=1/3))


## ----message=F, warning=F, fig.width=8, fig.height=3, fig.fullwidth=T, dev='svg'----
plot(dnbinom(1:40, size=10, prob=1 - 1/3))


## ----message=F, warning=F, fig.width=8, fig.height=4, fig.fullwidth=T, dev='svg'----
hist(rnbinom(100, size=1e+7, prob=1 - 1e-6))


## ----message=F, warning=F, fig.width=8, fig.height=4, fig.fullwidth=T, dev='svg'----
hist(rnbinom(100, size=1e+6, prob=1 - 1e-6))


## ----message=F, warning=F, fig.width=8, fig.height=4, fig.fullwidth=T, dev='svg'----
hist(rnbinom(100, size=1e+5, prob=1 - 1e-6))


## ----message=FALSE, warning=F, fig.show = 'hold', fig.width=6, fig.height=3, fig.fullwidth=T, dev='svg'----
res <- results(ddsMerged)
res <- as.data.frame(res)
res$gene <- rownames(res)
ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=padj < 0.05)) +
  geom_point() + theme_bw() + scale_color_manual(values=c("black", "red"))
  # geom_text_repel(data=res %>% dplyr::filter(padj < 1e-40), aes(label=gene, color=NULL))



## ----message=FALSE, warning=F, fig.show = 'hold', fig.width=6, fig.height=3, fig.fullwidth=T, dev='svg'----
res <- lfcShrink(ddsMerged, coef="strain_wt_vs_snf2_ko", type="apeglm")
res <- as.data.frame(res)
res$gene <- rownames(res)
ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=padj < 0.05)) +
  geom_point() + theme_bw() + scale_color_manual(values=c("black", "red"))
  # geom_text_repel(data=res %>% dplyr::filter(padj < 1e-40), aes(label=gene, color=NULL))

