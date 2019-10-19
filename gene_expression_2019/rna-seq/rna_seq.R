## ----message=F-----------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("apeglm", quietly = TRUE)) BiocManager::install("apeglm")
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db")




## ----message=FALSE, warning=FALSE----------------------------------------
library(ggplot2)
library(DESeq2)
library(apeglm)
library(ggrepel)
library(dplyr)
library(org.Mm.eg.db)


countFiles <- list.files("GSE116239_RAW", full.names = T)
countFiles


## ----message=FALSE-------------------------------------------------------

counts <- lapply(countFiles, function(countsFile) {
  read.table(countsFile, sep="\t", header=1, row.names = 1, stringsAsFactors = F, comment.char = "")
})

head(counts[[1]])



## ----message=FALSE-------------------------------------------------------

counts <- lapply(counts, function(countsTable) countsTable[, "Count", drop=F])
counts <- do.call(cbind, counts)
colnames(counts) <- gsub(".*(GSM\\d+).*", "\\1", countFiles)
head(counts)



## ----message=FALSE-------------------------------------------------------

coldata <- data.frame(
  gsm=gsub(".*(GSM\\d+).*", "\\1", countFiles),
  foam=gsub(".*(GSM\\d+)_(foam|non_foam).*", "\\2", countFiles),
  row.names =gsub(".*(GSM\\d+).*", "\\1", countFiles)
)

coldata



## ----message=FALSE-------------------------------------------------------

# only keeping genes that have at least 10 reads
dds <- DESeqDataSetFromMatrix(countData = counts[rowSums(counts) > 10, ],
                              colData = coldata,
                              design= ~ foam)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients



## ----fig.show = 'hold', fig.width=3, fig.height=3, fig.fullwidth=T, dev='svg'----

vst <- varianceStabilizingTransformation(dds)
plotPCA(vst, intgroup=c("foam"))




## ----message=FALSE, warning=FALSE----------------------------------------

res <- lfcShrink(dds, coef="foam_non_foam_vs_foam", type="apeglm")
head(res)


## ----message=FALSE, warning=FALSE----------------------------------------

keytypes(org.Mm.eg.db)
res$Gene.symbol <- mapIds(org.Mm.eg.db, rownames(res), column="SYMBOL", "ENSEMBL")



## ----message=FALSE, warning=F, fig.show = 'hold', fig.width=6, fig.height=3, fig.fullwidth=T, dev='svg'----
res <- as.data.frame(res)
ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=padj < 0.05)) +
  geom_point() + theme_bw() + scale_color_manual(values=c("black", "red")) +
  geom_text_repel(data=res %>% dplyr::filter(padj < 1e-20), aes(label=Gene.symbol, color=NULL))



