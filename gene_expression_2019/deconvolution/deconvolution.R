## ----message=F, warning=F------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("Biobase", quietly = TRUE)) BiocManager::install("Biobase")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("linseed", quietly = TRUE)) devtools::install_github("https://github.com/ctlab/LinSeed")

library(GEOquery)
library(Biobase)
library(limma)
library(dplyr)
library(pheatmap)
library(linseed)
library(ggplot2)
library(reshape2)



## ----warning=FALSE, message=FALSE, cache=TRUE----------------------------


gse19830 <- getGEO("GSE19830", AnnotGPL = TRUE)[[1]]
head(pData(gse19830)[, c("geo_accession", "characteristics_ch1")])



## ----warning=FALSE, message=FALSE----------------------------------------

pdata <- pData(gse19830)
proportionsCharacter <- pdata$characteristics_ch1
proportions <- data.frame(
  liver=as.numeric(gsub(".* (\\d+) % Liver.*", "\\1", proportionsCharacter)),
  brain=as.numeric(gsub(".* (\\d+) % Brain.*", "\\1", proportionsCharacter)),
  lung=as.numeric(gsub(".* (\\d+) % Lung.*", "\\1", proportionsCharacter)),
  row.names=rownames(pdata)
)





## ----warning=FALSE, message=FALSE----------------------------------------

head(proportions)



## ----warning=FALSE, message=FALSE----------------------------------------

head(exprs(gse19830)[, 1:5])



## ------------------------------------------------------------------------

fData(gse19830) <- fData(gse19830)[, c("ID", "Gene symbol", "Gene ID")]
head(fData(gse19830))



## ------------------------------------------------------------------------

exprs(gse19830) <- normalizeBetweenArrays(exprs(gse19830), method="quantile")



## ------------------------------------------------------------------------
gse19830 <- gse19830[!grepl("///", fData(gse19830)$`Gene symbol`), ]
gse19830 <- gse19830[fData(gse19830)$`Gene symbol` != "", ]

fData(gse19830)$mean_expression <- apply(exprs(gse19830), 1, mean)
gse19830 <- gse19830[order(fData(gse19830)$mean_expression, decreasing = TRUE), ]
gse19830 <- gse19830[!duplicated(fData(gse19830)$`Gene ID`), ]
gse19830 <- gse19830[seq_len(12000), ]
dim(gse19830)


## ------------------------------------------------------------------------
pure <- gse19830[, 1:9]
pureProportions <- proportions[1:9, ]
pData(pure)$Liver <- c("NonLiver", "Liver")[as.numeric(pureProportions$liver == 100) + 1]
pData(pure)$Brain <- c("NonBrain", "Brain")[as.numeric(pureProportions$brain == 100) + 1]
pData(pure)$Lung <- c("NonLung", "Lung")[as.numeric(pureProportions$lung == 100) + 1]


## ------------------------------------------------------------------------
liver.model <- model.matrix(~0 + Liver, data=pData(pure))
fit <- lmFit(pure, liver.model)

fit2 <- contrasts.fit(fit, makeContrasts(LiverLiver - LiverNonLiver, levels=liver.model))
fit2 <- eBayes(fit2, trend = T)

de <- topTable(fit2, adjust.method="BH", number=Inf, sort.by = "P")


## ------------------------------------------------------------------------
head(de)


## ----warning=FALSE, message=FALSE----------------------------------------


markers.liver <- de %>% filter(adj.P.Val < 0.01) %>% top_n(50, logFC)
head(markers.liver)



## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F-----------------
markers.liver <- markers.liver %>% pull(ID)
pheatmap(pure[markers.liver, ], scale = "row", cluster_rows = F, cluster_cols = F)


## ------------------------------------------------------------------------
brain.model <- model.matrix(~0 + Brain, data=pData(pure))
fit <- lmFit(pure, brain.model)

fit2 <- contrasts.fit(fit, makeContrasts(BrainBrain - BrainNonBrain, levels=brain.model))
fit2 <- eBayes(fit2, trend = T)

de <- topTable(fit2, adjust.method="BH", number=Inf, sort.by = "P")
markers.brain <- de %>% filter(adj.P.Val < 0.01) %>% top_n(50, logFC)
markers.brain <- markers.brain %>% pull(ID)


## ------------------------------------------------------------------------
lung.model <- model.matrix(~0 + Lung, data=pData(pure))
fit <- lmFit(pure, lung.model)

fit2 <- contrasts.fit(fit, makeContrasts(LungLung - LungNonLung, levels=lung.model))
fit2 <- eBayes(fit2, trend = T)

de <- topTable(fit2, adjust.method="BH", number=Inf, sort.by = "P")
markers.lung <- de %>% filter(adj.P.Val < 0.01) %>% top_n(50, logFC)
markers.lung <- markers.lung %>% pull(ID)


## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F-----------------

pheatmap(pure[c(markers.liver, markers.brain, markers.lung), ], scale = "row", cluster_rows = F, cluster_cols = F)


## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F-----------------
pheatmap(t(proportions), cluster_rows = F, cluster_cols = F)


## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F-----------------
pheatmap(2^exprs(gse19830)[c(markers.liver[1:3], markers.brain[1:3], markers.lung[1:3]), ], scale="row", cluster_rows = F, cluster_cols = F)


## ----fig.height=2.5, fig.fullwidth=T, dev='svg', message=F---------------
mix <- exprs(gse19830)[, 10:42]
mix <- 2^mix # going to linear scale
markers <- list(
  liver=markers.liver,
  brain=markers.brain,
  lung=markers.lung
)
hist(mix[, 1])



## ------------------------------------------------------------------------

dsaResults <- linseed:::fastDSA(mix, markers)
head(dsaResults$W)


## ------------------------------------------------------------------------

head(dsaResults$H[, 1:6])


## ------------------------------------------------------------------------

actual <- as.matrix(proportions[10:42, ]) / 100
predicted <- t(dsaResults$H)
props <- data.frame(
  Actual=as.numeric(actual),
  Predicted=as.numeric(predicted),
  Tissue=c(rep("liver", ncol(mix)), rep("brain", ncol(mix)), rep("lung", ncol(mix)))
)
head(props)



## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F-----------------
ggplot(data=props, aes(x=Actual, y=Predicted, color=Tissue)) +
  geom_point() + geom_abline(slope=1, intercept = 0, lty=2) +
  facet_grid(~Tissue) + theme_bw() + theme(aspect.ratio = 1)



## ------------------------------------------------------------------------
# fancy plots: basis
actualBasis <- cbind(
  2^rowMeans(exprs(pure)[, 1:3]),
  2^rowMeans(exprs(pure)[, 4:6]),
  2^rowMeans(exprs(pure)[, 7:9])
)
colnames(actualBasis) <- c("Liver", "Brain", "Lung")
predictedBasis <- dsaResults$W

plotData <- data.frame(
  Actual=as.numeric(actualBasis),
  Predicted=as.numeric(predictedBasis),
  Tissue=c(rep("Liver", nrow(exprs(pure))), rep("Brain", nrow(exprs(pure))), rep("Lung", nrow(exprs(pure))))
)


## ----fig.height=3, fig.fullwidth=T, dev='svg', message=FALSE, warning=FALSE----
ggplot(data=plotData, aes(x=Actual, y=Predicted)) + 
  geom_hex(bins=100) + theme_bw() + facet_grid(~Tissue) +
  geom_abline(slope=1, intercept = 0, lty=2) + scale_x_log10() + scale_y_log10() +
  theme(aspect.ratio = 1)

