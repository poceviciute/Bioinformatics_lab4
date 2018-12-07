## ------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GEOquery", version = "3.8")

library(GEOquery)

## ---- cache=TRUE---------------------------------------------------------
x = getGEOSuppFiles("GSE20986")
x
setwd("C:/Users/Milda/Documents/Documents/Master degree/Bioinformatics/Bioinformatics_lab4")
## ------------------------------------------------------------------------
untar("GSE20986/GSE20986_RAW.tar", exdir = "data")

## ------------------------------------------------------------------------
cels = list.files("data/", pattern = "[gz]")
sapply(paste("data", cels, sep = "/"), gunzip)

## ------------------------------------------------------------------------
phenodata = matrix(rep(list.files("data"), 2), ncol =2)
class(phenodata)
phenodata <- as.data.frame(phenodata)
colnames(phenodata) <- c("Name", "FileName")
phenodata$Targets <- c("iris", 
                       "retina", 
                       "retina", 
                       "iris", 
                       "retina", 
                       "iris", 
                       "choroid", 
                       "choroid", 
                       "choroid", 
                       "huvec", 
                       "huvec", 
                       "huvec")
write.table(phenodata, "data/phenodata.txt", quote = F, sep = "\t", row.names = F)



## ------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("simpleaffy", version = "3.8")
library(simpleaffy)

celfiles <- read.affy(covdesc = "phenodata.txt", path = "data")
boxplot(celfiles)
install.packages("RColorBrewer")
library(RColorBrewer)
cols = brewer.pal(8, "Set1")
eset <- exprs(celfiles)
samples <- celfiles$Targets
colnames(eset)
colnames(eset) <- samples

boxplot(celfiles, col = cols, las = 2)
distance <- dist(t(eset), method = "maximum")
clusters <- hclust(distance)
plot(clusters)



## ------------------------------------------------------------------------

library(simpleaffy)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("affyPLM", version = "3.8")
library(affyPLM)

celfiles.gcrma = gcrma(celfiles)
par(mfrow=c(1,2))
boxplot(celfiles.gcrma, col = cols, las = 2, main = "Post-Normalization");
boxplot(celfiles, col = cols, las = 2, main = "Pre-Normalization")
dev.off()

distance2 <- dist(t(exprs(celfiles.gcrma)), method = "maximum")
clusters2 <- hclust(distance2)
plot(clusters2)

## ------------------------------------------------------------------------
library(limma)

phenodata
samples <- as.factor(samples)
design <- model.matrix(~0+samples)
colnames(design)
colnames(design) <- c("choroid", "huvec", "iris", "retina")
design
contrast.matrix = makeContrasts(
    huvec_choroid = huvec - choroid, 
    huvec_retina = huvec - retina, 
    huvec_iris <- huvec - iris, 
    levels = design)

fit = lmFit(celfiles.gcrma, design)
huvec_fit <- contrasts.fit(fit, contrast.matrix)
huvec_ebay <- eBayes(huvec_fit)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("hgu133plus2.db", version = "3.8")
library(hgu133plus2.db)
library(annotate)
probenames.list <- rownames(topTable(huvec_ebay, number = 100000))
getsymbols <- getSYMBOL(probenames.list, "hgu133plus2")
results <- topTable(huvec_ebay, number = 100000, coef = "huvec_choroid")
results <- cbind(results, getsymbols)

## ------------------------------------------------------------------------
summary(results)

results$threshold <- "1"
a <- subset(results, adj.P.Val < 0.05 & logFC > 5)
results[rownames(a), "threshold"] <- "2"
b <- subset(results, adj.P.Val < 0.05 & logFC < -5)
results[rownames(b), "threshold"] <- "3"
table(results$threshold)

## ------------------------------------------------------------------------
library(ggplot2)
volcano <- ggplot(data = results, 
                  aes(x = logFC, y = -1*log10(adj.P.Val), 
                      colour = threshold, 
                      label = getsymbols))

volcano <- volcano + 
    geom_point() + 
    scale_color_manual(values = c("black", "red", "green"), 
                       labels = c("Not Significant", "Upregulated", "Downregulated"), 
                       name = "Key/Legend")

volcano + 
    geom_text(data = subset(results, logFC > 5 & -1*log10(adj.P.Val) > 5), aes(x = logFC, y = -1*log10(adj.P.Val), colour = threshold, label = getsymbols)  )

