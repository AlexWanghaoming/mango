load("/Users/alexwang/r_env/r_data/DEG.RData")
# load("/Users/alexwang/r_env/.RData")
library(edgeRun)
library(GenomicFeatures)
rawReadsCounts <- read.delim("/Users/alexwang/data/diffExpression/readcount.xls",header = T, stringsAsFactors = F)
rownames(rawReadsCounts) <- rawReadsCounts[,1]   
# raw_data2 <- read.delim("/home/wanghm/wys/wanghmDE_seq/id_fpkm1.xls",header = T)
# rawReadsCounts <- raw_data1[raw_data2$geneID,]
# write.table(rawReadsCounts, "rawReadsCounts.txt",quote = F)
rawReadsCounts[,1] <- NULL  ## 33124 genes
xx <- sapply(strsplit(colnames(rawReadsCounts), "_"), function(x) x[1]) 
# a <- DGEList(counts = rawReadsCounts, group = xx)   # 实验设计需要几组数据，rowReadsCount就提取需要的几组数据
grp <- xx[c(3,4,7,8,15,16,19,20)]
a <- DGEList(counts = rawReadsCounts[,c(3,4,7,8,15,16,19,20)], group = grp)
order(grp)

# colnames(a);rm(y)
keep <- rowSums(cpm(a)>1) >=2
y <- a[keep,]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y,method = "TMM")
y <- estimateCommonDisp(y, verbose = TRUE)
y <- estimateTagwiseDisp(y) 

# QC plot
out_dir <- getwd()
pdf(paste(out_dir, "/res/QC.pdf", sep=""))
pairs(a$counts, pch=".", upper.panel = NULL) ## 给出八种样品见的关系一览概况
plotBCV(y)
plotMDS(y, method = "bcv")
dev.off()

# get DEGS
exactT
aaa_uc <- UCexactTest(y,pair = c("r1P24h","WTP24h"))
re <- decideTestsDGE(aaa_uc, adjust.method = "BH", p.value = 0.05, lfc)
summary(re)
detags <- rownames(y)[as.logical(re)]
uptags<-rownames(y)[as.logical(re>0)]
downtags<-rownames(y)[as.logical(re<0)]
## 准备火山图数据
GENE <- rownames(aaa_uc)
Pvalue <- aaa_uc$table$PValue
FC <- exp(aaa_uc$table$logFC)
FDR <- p.adjust(aaa_uc$table$PValue,method = "BH")
absFC <- abs(FC)
logFC <- aaa_uc$table$logFC
# # plot result
# out_dir <- getwd()
# pdf(paste(out_dir, "/res/res.pdf", sep=""))
# plotSmear(aaa_uc, de.tags=detags)
# abline(h=c(-1, 1), col="blue")
# EnhancedVolcano(volcanodata,x = "logFC",y = "FDR",lab = rownames(volcanodata),ylim = c(0,6))
# dev.off()

# Import data
volcanodata <- data.frame(GENE,Pvalue,FC,FDR,absFC,logFC)
colnames(volcanodata)
# Set up edgeR compatible column names
colnames(volcanodata) <- c("Gene","pvalue","FC","FDR","absFC","logFC")
# Set transcript names as row names
row.names(volcanodata)<-make.names(volcanodata[,1], unique=TRUE)
EnhancedVolcano(volcanodata,x = "logFC",y = "FDR",lab = rownames(volcanodata),ylim = c(0,6), FCcutoff = 1)

### multigroup output txt
# all_combination <- combn(1:10,2)  # free combination
# ss <- function(x){
#   aaa <- exactTest(y, pair = c(x[1], x[2]))
#   bb <- topTags(aaa, p.value = 0.01)
#   return(bb)
# }
# xx <- apply(all_combination,2, ss)
# sink("wys_diffGenes")
# for (i in xx) {
#   print(i)
# }
# sink()


#######*********############
# Analysis with DESeq2 ----- ource: github:stephenturner's DESeq2 pipline
library(DESeq2)
xx <- sapply(strsplit(colnames(rawReadsCounts), "_"), function(x) x[1])
xx[12] <- "rtp1M5h"
yy <- c("r1P24h","r1P24h", "WTP24h", "WTP24h")
condition <- factor(yy)
countData <- rawReadsCounts[,c(7,8,19,20)]
colData <- data.frame(row.names = colnames(countData), condition) 
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
# 对原始dds进行标准化
dds2 <- DESeq(dds)
summary(dds2)
# get DEGs
sampleA <- "r1P24h"
sampleB <- "WTP24h"
contrastV <- c("condition", sampleA, sampleB)
res <- results(dds2,contrast = contrastV,lfcThreshold = 0)
res <- results(dds2,lfcThreshold = 0)
# 校正后p-value为NA的复制为1
res$padj[is.na(res$padj)] <- 1
# 按pvalue排序, 把差异大的基因放前面
res <- res[order(res$pvalue),]
subset(res, res$log2FoldChange >= 1)
summary(res)

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds2)
head(assay(rld))
hist(assay(rld),plot = T)

# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))]

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Principal components analysis
## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup="condition")
## I like mine better:
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  textx
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()


# Get differential expression results
res <- results(dds2)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
summary(res)
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results.csv")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## Examine independent filtering
attr(res, "filterThreshold")
plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

# install.packages("calibrate")
## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()



