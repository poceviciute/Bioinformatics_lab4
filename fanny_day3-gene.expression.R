## ------------------------------------------------------------------------
library(GEOquery)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GEOquery", version = "3.8")
## ---- cache=TRUE---------------------------------------------------------
x = getGEOSuppFiles("GSE20986")
x

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
library(simpleaffy)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("simpleaffy", version = "3.8")

celfiles <- read.affy(covdesc = "phenodata.txt", path = "data")
#dim(celfiles)
boxplot(celfiles)
library(RColorBrewer)
cols = brewer.pal(8, "Set1")
eset <- exprs(celfiles)
#dim(eset)
samples <- celfiles$Targets
colnames(eset)
colnames(eset) <- samples

boxplot(celfiles, col = cols, las = 2)
distance <- dist(t(eset), method = "maximum") #finds largest distance between columns
#length(distance)
clusters <- hclust(distance) 
#This function performs a hierarchical cluster analysis using a set of dissimilarities 
#for the n objects being clustered. 
#Initially, each object is assigned to its own cluster and then the algorithm 
#proceedSs iteratively, at each stage joining the two most similar clusters, 
#continuing until there is just a single cluster. 
#At each stage distances between clusters are recomputed by the Lance-Williams 
#dissimilarity update formula according to the particular clustering method being used.
#The complete linkage method finds similar clusters.
plot(clusters)



## ------------------------------------------------------------------------
require(simpleaffy)
require(affyPLM)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("affyPLM", version = "3.8")
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
              huvec_iris = huvec - iris, 
              levels = design)

fit = lmFit(celfiles.gcrma, design)
#Fit linear model for each gene given a series of arrays
huvec_fit <- contrasts.fit(fit, contrast.matrix)
#Given a linear model fit to microarray data, compute estimated 
#coefficients and standard errors for a given set of contrasts.
huvec_ebay <- eBayes(huvec_fit)
#Given a microarray linear model fit, compute moderated t-statistics, 
#moderated F-statistic, and log-odds of differential expression by 
#empirical Bayes moderation of the standard errors towards a common value.

library(hgu133plus2.db)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("hgu133plus2.db", version = "3.8")
library(annotate)
probenames.list <- rownames(topTable(huvec_ebay, number = 100000))
#Extract a table of the top-ranked genes from a linear model fit.
getsymbols <- getSYMBOL(probenames.list, "hgu133plus2")
#The functions documented here are intended to make it easier to map 
#from a set of manufacturers identifiers (such as you will get from the chips etc) 
#to other identifiers.
results <- topTable(huvec_ebay, number = 100000, coef = "huvec_choroid")
results <- cbind(results, getsymbols)
#head(results)

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


########## Q2

library(preprocessCore)
#1=iris, 2=retina
#7=choroid, huvec=10
ir <- which(colnames(eset)=="iris")[1]
ret <- which(colnames(eset)=="retina")[1]
cho <- which(colnames(eset)=="choroid")[1]
huv <- which(colnames(eset)=="huvec")[1]
log2eset <- log2(eset)
# Normalization
eset_norm<- exprs(celfiles.gcrma)
colnames(eset_norm) <- samples
log2eset_norm <- log2(eset_norm)

par(mfrow=c(2,2))
#ir, ret
ma.plot(rowMeans(log2eset[,c(ir,ret)]), log2eset[,ir]-log2eset[,ret], main="Pre-norm (ir,ret)", cex=1)
ma.plot(rowMeans(log2eset_norm[,c(ir, ret)]), log2eset_norm[,ir]-log2eset_norm[,ret], main="Post-norm (ir,ret)",cex=1)

#ir, cho
ma.plot(rowMeans(log2eset[,c(ir,cho)]), log2eset[,ir]-log2eset[,cho], main="Pre-norm (ir,cho)", cex=1)
ma.plot(rowMeans(log2eset_norm[,c(ir, cho)]), log2eset_norm[,ir]-log2eset_norm[,cho], main="Post-norm (ir,cho)",cex=1)

#ir, huv
ma.plot(rowMeans(log2eset[,c(ir,huv)]), log2eset[,ir]-log2eset[,huv], main="Pre-norm (ir,huv)", cex=1)
ma.plot(rowMeans(log2eset_norm[,c(ir, huv)]), log2eset_norm[,ir]-log2eset_norm[,huv], main="Post-norm (ir,huv)",cex=1)

#ret, cho
ma.plot(rowMeans(log2eset[,c(ret,cho)]), log2eset[,ret]-log2eset[,cho], main="Pre-norm (ret,cho)", cex=1)
ma.plot(rowMeans(log2eset_norm[,c(ret, cho)]), log2eset_norm[,ret]-log2eset_norm[,cho], main="Post-norm (ret,cho)",cex=1)

#ret, huv
ma.plot(rowMeans(log2eset[,c(ret,huv)]), log2eset[,ret]-log2eset[,huv], main="Pre-norm (ret,huv)", cex=1)
ma.plot(rowMeans(log2eset_norm[,c(ret, huv)]), log2eset_norm[,ret]-log2eset_norm[,huv], main="Post-norm (ret,huv)",cex=1)

#cho, huv
ma.plot(rowMeans(log2eset[,c(cho,huv)]), log2eset[,cho]-log2eset[,huv], main="Pre-norm (cho, huv)", cex=1)
ma.plot(rowMeans(log2eset_norm[,c(cho, huv)]), log2eset_norm[,cho]-log2eset_norm[,huv], main="Post-norm (cho, huv)",cex=1)
