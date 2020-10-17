## ----message=F------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("Biobase", quietly = TRUE)) BiocManager::install("Biobase")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")

library(GEOquery)
library(Biobase)
library(ggplot2)
library(reshape2)
library(limma)
library(MASS)


## ----cache=T, message=F---------------------------------------------------------------

GSE53986 <- getGEO("GSE53986", AnnotGPL = TRUE)[[1]]



## -------------------------------------------------------------------------------------

colnames(pData(GSE53986))



## ----message=F------------------------------------------------------------------------
pData(GSE53986)$rep <- gsub(".*, (\\d)$", "rep \\1", pData(GSE53986)$title)
pData(GSE53986) <- pData(GSE53986)[, c("cell type:ch1", "treatment:ch1", "rep")]

colnames(pData(GSE53986)) <- c("Cell", "Treatment", "Replicate")
head(pData(GSE53986))



## -------------------------------------------------------------------------------------
pData(GSE53986)$LPS <- as.factor(c("no", "yes")[grepl("LPS", pData(GSE53986)$Treatment) + 1])
pData(GSE53986)$IFNg <- as.factor(c("no", "yes")[grepl("IFNg", pData(GSE53986)$Treatment) + 1])
head(pData(GSE53986))


## ----cache=T, message=F---------------------------------------------------------------
colnames(fData(GSE53986))



## ----message=F------------------------------------------------------------------------
fData(GSE53986) <- fData(GSE53986)[, c("ID", "Gene symbol", "Gene ID")]
head(fData(GSE53986))



## ----fig.height=3, fig.fullwidth=T, dev='svg'-----------------------------------------

ggplot(data=data.frame(expression=exprs(GSE53986)[, 1]),
       aes(x=expression)) +
  geom_histogram()



## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F------------------------------

ggplot(data=data.frame(expression_log2=log2(exprs(GSE53986)[, 1])),
       aes(x=expression_log2)) +
  geom_histogram()



## -------------------------------------------------------------------------------------
min(exprs(GSE53986))


## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F------------------------------

twoSamples <- melt(exprs(GSE53986[, 1:2]))
twoSamples$value <- log2(twoSamples$value)

ggplot(data=twoSamples, aes(x=value)) +
  facet_grid(~Var2) + geom_histogram()



## -------------------------------------------------------------------------------------
colSums(exprs(GSE53986))


## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F------------------------------

exprs(GSE53986) <- normalizeBetweenArrays(log2(exprs(GSE53986)+1), method="quantile")
twoSamples <- melt(exprs(GSE53986[, 1:2]))

ggplot(data=twoSamples, aes(x=value)) +
  facet_grid(~Var2) + geom_histogram()



## ----eval=F---------------------------------------------------------------------------
## head(fData(GSE53986), 1000)


## -------------------------------------------------------------------------------------
GSE53986 <- GSE53986[!grepl("///", fData(GSE53986)$`Gene symbol`), ]
GSE53986 <- GSE53986[fData(GSE53986)$`Gene symbol` != "", ]

fData(GSE53986)$mean_expression <- apply(exprs(GSE53986), 1, mean)
GSE53986 <- GSE53986[order(fData(GSE53986)$mean_expression, decreasing = TRUE), ]
GSE53986 <- GSE53986[!duplicated(fData(GSE53986)$`Gene ID`), ]
GSE53986 <- GSE53986[seq_len(12000), ]
dim(GSE53986)


## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F------------------------------
pcas <- prcomp(t(exprs(GSE53986)), scale. = T)
plotData <- cbind(pcas$x[, 1:2], pData(GSE53986))
ggplot(plotData, aes(x=PC1, y=PC2, color=Treatment)) +
  geom_point() + theme_bw() + theme(aspect.ratio = 1)



## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F------------------------------
ggplot(plotData, aes(x=PC1, y=PC2, color=Treatment)) +
  geom_point() +
  geom_text(aes(label=rownames(plotData))) + theme_bw() + theme(aspect.ratio = 1)



## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F------------------------------
variance <- pcas$sdev^2
ggplot(data=data.frame(component=1:16, variance=variance),
       aes(x=component, y=variance)) +
  geom_point() + geom_line() + theme_bw()



## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F------------------------------
variance <- variance / sum(variance)
ggplot(data=data.frame(component=1:16, percent=variance * 100),
       aes(x=component, y=percent)) +
  geom_point() + geom_line() + theme_bw()



## -------------------------------------------------------------------------------------

GSE53986.design <- model.matrix(~0+LPS+IFNg, data=pData(GSE53986))
colnames(GSE53986.design) <- c("LPSno", "LPSyes", "IFNgyes")

fit <- lmFit(GSE53986, GSE53986.design)

fit2 <- contrasts.fit(fit, makeContrasts(LPSyes - LPSno, levels=GSE53986.design))
fit2 <- eBayes(fit2, trend = T)

de <- topTable(fit2, adjust.method="BH", number=Inf, sort.by = "P")



## -------------------------------------------------------------------------------------
head(de)


## -------------------------------------------------------------------------------------

GSE53986.design <- model.matrix(~0+IFNg+LPS, data=pData(GSE53986))
colnames(GSE53986.design) <- c("IFNgno", "IFNgyes", "LPSyes")

fit <- lmFit(GSE53986, GSE53986.design)

fit2 <- contrasts.fit(fit, makeContrasts(IFNgyes-IFNgno, levels=GSE53986.design))
fit2 <- eBayes(fit2, trend = T)

de <- topTable(fit2, adjust.method="BH", number=Inf, sort.by = "P")



## -------------------------------------------------------------------------------------
head(de)


## BELOW ARE ILLUSTRATIONS ON HOW LINEAR MODELS WORK


## -------------------------------------------------------------------------------------
lps_model <- model.matrix(~0 + LPS, data = pData(GSE53986))
colnames(lps_model) <- c("LPS_no", "LPS_yes")
lps_model


## -------------------------------------------------------------------------------------
## il1rn
exprs(GSE53986)["1451798_at", ]
linear_fit <- lm.fit(lps_model, exprs(GSE53986)["1451798_at", ])
linear_fit$coef


## -------------------------------------------------------------------------------------
ifng_model <- model.matrix(~0 + IFNg, data = pData(GSE53986))
colnames(ifng_model) <- c("ifng_no", "ifng_yes")
ifng_model


## -------------------------------------------------------------------------------------
## il1rn
exprs(GSE53986)["1451798_at", ]
linear_fit <- lm.fit(ifng_model, exprs(GSE53986)["1451798_at", ])
linear_fit$coef


## -------------------------------------------------------------------------------------
## Serpina3g
exprs(GSE53986)["1424923_at", ]
linear_fit <- lm.fit(ifng_model, exprs(GSE53986)["1424923_at", ])
linear_fit$coef


## -------------------------------------------------------------------------------------
treatment_model <- model.matrix(~0 + LPS + IFNg, data = pData(GSE53986))
colnames(treatment_model) <- c("LPS_no", "LPS_yes", "IFNg_yes")
treatment_model


## -------------------------------------------------------------------------------------
## il1rn
exprs(GSE53986)["1451798_at", ]
linear_fit <- lm.fit(treatment_model, exprs(GSE53986)["1451798_at", ])
linear_fit$coef


## -------------------------------------------------------------------------------------
# serpina3g
exprs(GSE53986)["1424923_at", ]
linear_fit <- lm.fit(treatment_model, exprs(GSE53986)["1424923_at", ])
linear_fit$coef


## -------------------------------------------------------------------------------------
linear_fit <- lm.fit(treatment_model, exprs(GSE53986)["1451798_at", ])
linear_fit$coef

linear_fit <- lm.fit(treatment_model, exprs(GSE53986)["1424923_at", ])
linear_fit$coef


## -------------------------------------------------------------------------------------
full_model <- model.matrix(~1 + LPS + IFNg, data = pData(GSE53986))
colnames(full_model) <- c("Intercept", "LPS_yes", "IFNg_yes")
full_model


## -------------------------------------------------------------------------------------
linear_fit <- lm.fit(full_model, exprs(GSE53986)["1451798_at", ])
linear_fit$coef

linear_fit <- lm.fit(full_model, exprs(GSE53986)["1424923_at", ])
linear_fit$coef

