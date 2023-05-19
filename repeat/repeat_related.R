## repeat region associate accessory chromosome and 2-speed regions
strain <- "yn55"
cc <- c(rgb(8,186,255,maxColorValue = 255), rgb(255,128,2, maxColorValue = 255))
library(regioneR)
library(GenomicFeatures)
repeat.gr <- toGRanges(read.table(paste0("/Users/alexwang/0data/0mango/genome/repeat/",strain,".repeat.bed"), header = T))
df.speed <- read.table(paste0("/Users/alexwang/0data/0mango/popGenome/2speed/", strain,"_state.tsv"))
# df.speed <- read.table(paste0("/Users/alexwang/0data/0mango/genome/2speed/", strain,"_EMstate.tsv"))
df.speed$V1 <-  gsub(pattern = paste0(strain,"_"), replacement = "", df.speed$V1)
speed.gr <- toGRanges(df.speed)
fast.gr <- speed.gr[speed.gr$V4 == 2,]
slow.gr <- speed.gr[speed.gr$V4 == 1,]

# calculate overlapped base of repeat && fast/slow regions
ovBase.fast <- overlapRegions(repeat.gr, fast.gr, get.bases = T)$ov.bases
# (sum(ovBase.fast)*10000)/sum(width(fast.gr))
m.fast <- mean(ovBase.fast)
ovBase.slow <- overlapRegions(repeat.gr, slow.gr, get.bases = T)$ov.bases
# (sum(ovBase.slow)*10000)/sum(width(slow.gr))
m.slow <- mean(ovBase.slow)

# chisq.test
p <- c(1/2, 1/2)
x <- c(m.fast, m.slow)
chisq.test(x = x, p = p)
## plot bar
cc <- c(rgb(8,186,255,maxColorValue = 255), rgb(255,128,2, maxColorValue = 255))
par(mar=c(1.5,2,5,2))
barplot(c(m.fast/1000, m.slow/1000), horiz = T, space = 0.5,col = cc, las=1, xaxt="n")
legend(x = 0.4, y = 0.5, legend = c("Fast", "Slow"), 
       col = cc, xpd = T, pch=15, bty = "n", ncol = 2)
axis(side = 3, cex=0.5)
text(x = 0.4, y=1.7, "*** X-squared p<2.2e-16")
title("Number of repeats / 10 kb")



### calculate per 10k repeat base on accessory and core chromosomes
strain = "gd10"
strains = c("fj11", "yn55", "yn56", "hn47", "qz", "gz15", "gd10")
cc <- c(rgb(8,186,255,maxColorValue = 255), rgb(255,128,2, maxColorValue = 255))

chrom.core.fj11 <- 1:10
chrom.mini.fj11 <- 11:16
chrom.core.yn55 <- 1:10
chrom.mini.yn55 <- 11:13
chrom.core.yn56 <- 1:10
chrom.mini.yn56 <- 11:17
chrom.core.hn47 <- 1:10
chrom.mini.hn47 <- 11:12
chrom.core.qz <- 1:10
chrom.mini.qz <- 11:13
chrom.core.gz15 <- 1:10
chrom.mini.gz15 <- 11:16
chrom.core.gd10 <- 1:10
chrom.mini.gd10 <- 11:14

library(regioneR)
library(GenomicFeatures)
repeat.gr <- toGRanges(read.table(paste0("/Users/alexwang/0data/0mango/genome/repeat/",strain,".repeat.bed"), header = T))
genome <- read.table(paste0("/Users/alexwang/0data/0mango/genome/mummer/genome/",strain,".genome.txt"), header = T)
genome$start = 1
genome <- genome[,c(1,3,2)]
genome.gr <- toGRanges(genome)
accessory.gr <- genome.gr[seqnames(genome.gr)  %in% get0(paste0("chrom.mini.", strain))]
core.gr <- genome.gr[seqnames(genome.gr)  %in% get0(paste0("chrom.core.", strain))]
ovBase.accessory <- overlapRegions(repeat.gr, accessory.gr, get.bases = T)$ov.bases
ovBase.core <- overlapRegions(repeat.gr, core.gr, get.bases = T)$ov.bases
a1 <- sum(ovBase.accessory)*10000/sum(width(accessory.gr))
a2 <- sum(ovBase.core)*10000/sum(width(core.gr))
chisq.test(x=c(a1,a2), p = c(0.5,0.5))

pdf(file = "Accessory_repeat.pdf", width = 4, height = 6)
par(mar=c(1.5, 2, 3, 2))
barplot(c(a1/1000, a2/1000), horiz = F, space = 0.5,col = cc,xaxt="n", ylim = c(0, 1))
legend(x = 2, y = 0.5, legend = c("Accessory", "Core"),
       col = cc, xpd = T, pch=15, bty = "n", ncol = 2)
axis(side = 2, cex=0.5)
lines(c(1, 2.5), c(0.15,0.15), lwd = 2)
lines(c(1,1), c(0.17, 0.13), lwd = 2)
lines(c(2.5, 2.5), c(0.17, 0.13), lwd = 2)
title("Number of repeats / 10 kb")
dev.off()

### permutation test of the location correlation secretome to closest repeat
library(GenometriCorr)
library(regioneR)
strain <- "yn55"
repeat.df <- readr::read_table(paste0("/Users/alexwang/0data/0mango/LTR/LTR_repeatmasker/",strain,".fasta.out"),
                                         skip = 1)[,c(5,6,7)]
repeat.gr <- GRanges(seqnames = repeat.df$sequence, ranges = IRanges(start = repeat.df$begin, end = repeat.df$end),strand = "*")

txdb <- makeTxDbFromGFF(paste0("/Users/alexwang/0data/0mango/GTF/",strain,"_EVM.gtf"))
gene <- genes(txdb)
# effector.id <- read.table(paste0("/Users/alexwang/0data/0mango/inf_associated_gene/secretome/",strain,"_effector.geneID"))[,1]
effector.id <- read.table(paste0("/Users/alexwang/0data/0mango/inf_associated_gene/CSEP/",strain,".CSEP.geneID"))[,1]

effector.gr <- gene[gene$gene_id %in% effector.id]

sm.geneID <- read.table(paste0("/Users/alexwang/0data/0mango/inf_associated_gene/SM/",strain,"_SM.geneID"),fill = T)[,1]
sm.gr <- gene[gene$gene_id %in% sm.geneID]
gs <- read.table(paste0("~/0data/0mango/genome/mummer/genome/", strain, ".genome.txt"), header = T)
gs.len <- gs$size
names(gs.len) <- paste0(strain,"_",gs$chr)

### 统计检验

geCor <- GenometriCorr::GenometriCorrelation(query = effector.gr, reference = repeat.gr,
                                             permut.number = 100, keep.distributions = TRUE,
                                             chromosomes.to.proceed = names(gs.len), chromosomes.length = gs.len)
graphical.report(geCor, pdffile = paste0(strain,"_repeatPermutationTest2.SM.pdf"))

########################################################################## plot LTR percent
cc <- colorspace::divergex_hcl(2)
data = cbind(fj11=c(9.49,14.63-9.49), 
             yn55=c(8.08,14.14-8.08), 
             gz15=c(3.89, 6.18-3.89),
             hn47=c(0.45,4.14-0.45),
             yn56=c(2.42, 5.15-2.42),
             gd10=c(1.7,3.54-1.7))
mid = c(9.49/2, 8.08/2, 3.89/2, 0.45/2, 2.42/2, 1.7/2)
interval <- c(0.5,0.1,0.5,0.1,0.5,0.1)
barplot(data,space = interval,yaxt = "n",col = cc)
axis(side = 2,at=c(2,4,6,8,10,12,14),
     labels = c("%2","%4","%6","%8","%10","%12","%14"))
# s <- cumsum(interval[2:length(interval)])
text(c(0.5+0.5,2+0.1, 3+0.6, 4+0.7, 5+1.2, 6+1.3),mid, c(2049,1292,1029,179,734,489))
legend("topright", legend = c("LTR Repeat","Other Repeat"),
       border = F,pch = 15, col = cc, box.lwd = NA)



