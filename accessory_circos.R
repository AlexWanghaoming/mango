library(circlize)
library(GenomicFeatures)
library(RColorBrewer)
mini.fai <- read.table("/Users/alexwang/0data/0mango/synteny/accessory/accessory.genome.all.fasta.fai")
mini.gs <- mini.fai$V2
txdb <- makeTxDbFromGFF("/Users/alexwang/0data/0mango/synteny/accessory/accessory.all.gtf", format = "gtf")
ge <- genes(txdb)

chro <- mini.fai$V1
starts = rep(0, length(chro))
ends = mini.gs
genoCir <- data.frame(chr=chro, start=starts, end=ends)
genoCir$chr <- as.vector(genoCir[,1])

## start plot ##########################
circos.clear()
circos.par(start.degree = 87, track.height = 0.02, cell.padding = c(0,0,0,0), gap.degree = c(rep(1,nrow(genoCir)-1), 22))

## set plotType = "labels" if you would like to show names of contigs
circos.genomicInitialize(data = genoCir,
                         sector.names = chro,
                         labels.cex = 0.5, track.height = 0.05, plotType = NULL)
### a: ideagram of minichromosome
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  idx = which(chro == sector.index)
  circos.axis(labels.cex = 0.2,direction = "outside", labels.niceFacing = T, labels = "", minor.ticks = 5, lwd = 1, 
              major.at = c(0, mini.gs[idx]), major.tick.length = 0.8)
  }, track.height = 0.01, bg.border = NA)

####
# circos.trackPlotRegion(ylim = c(0, 0.5), ,panel.fun = function(x, y) {
#   sector.index = get.cell.meta.data("sector.index")
#   idx = which(chro == sector.index)
#   circos.genomicTrackPlotRegion(bg.col = cc1, track.height = 0.02)
# }, track.height = 0.01, bg.border = NA)

# extract color from RColorBrewer
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n_times.rep <- rle(gsub("_.*","",genoCir$chr))$lengths
s = length(n_times.rep)
cc1 = rep(sample(col_vector,size = s, replace = F), times=n_times.rep)
circos.genomicTrackPlotRegion(ylim = c(0,0.5),bg.col = cc1, bg.border = 1, track.height = 0.02)

# b: gene ortholog
###### extract cds
#python -m jcvi.formats.gff bed --type=gene --key=gene_id accessory.all.gtf -o accessory.bed
#gffread accessory.all.gtf -g accessory.genome.all.fasta -x accessory.cds.fasta
#awk '/>/{split($0,a,"=");print ">"a[2];next}{print}' accessory.cds.fasta > accessory.cds
#rm accessory.cds.fasta
###### get anchors
# python -m jcvi.compara.catalog ortholog accessory accessory --cscore=.90
# python -m jcvi.compara.synteny screen --minspan=10 --simple accessory.accessory.anchors accessory.accessory.anchors.new

ge.df <- data.frame(chr=seqnames(ge), start = start(ge), end = end(ge), id = ge$gene_id)

#### Using >10 gene block anchors
last.df <- read.table("/Users/alexwang/0data/0mango/synteny/accessory/accessory.accessory.anchors.simple")
last.df <- last.df[which(substr(last.df$V1, start = 1,stop = 5) != substr(last.df$V3, start = 1,stop = 5)),]

r1.start <- ge.df[match(last.df$V1, ge.df$id),]$start
r1.end <- ge.df[match(last.df$V2, ge.df$id), ]$end
r1.chr <- ge.df[match(last.df$V1, ge.df$id),]$chr
rr1 <- data.frame(chr = r1.chr, start = r1.start, r1.end = r1.end)
r2.start <- ge.df[match(last.df$V3, ge.df$id),]$start
r2.end <- ge.df[match(last.df$V4, ge.df$id),]$end
r2.chr <- ge.df[match(last.df$V3, ge.df$id),]$chr
rr2 <- data.frame(chr = r2.chr, start = r2.start, r1.end = r2.end)
#### Or using 1-to-1 gene anchors
# last.df <- read.table("/Users/alexwang/0data/0mango/synteny/accessory/accessory.accessory.anchors")
# g1 <- last.df$V1
# g2 <- last.df$V2
# r1 <- ge.df[match(g1, ge.df$id), c(1,2,3)]
# r2 <- ge.df[match(g2, ge.df$id), c(1,2,3)]
# bind.df <- cbind(r1,r2)
# rr1 <- bind.df[,c(1,2,3)]
# rr2 <- bind.df[,c(4,5,6)]

# assign color to ribbons
aa <- rr1$chr
color.map <- data.frame(chr=genoCir$chr , col = cc1)
ribbons.col <- sapply(aa, function(x) color.map[which(x == color.map$chr),][,2])
circos.genomicLink(rr1, rr2, border = NA, col =  paste0(ribbons.col, "50"), lwd = 2)

#### add legend
par(xpd=T)
legend(x = 1,y = 0.5, legend = unique(gsub("_.*","",color.map$chr)),
       fill = unique(cc1), border = NA,bty = "n", y.intersp = 0.55, x.intersp = 0.2)
## Or using ComplexHeatmap which is compatible with Circlize
# library(ComplexHeatmap)
# lgd <- Legend(labels = unique(gsub("_.*","",color.map$chr)), legend_gp = gpar(fill=unique(cc1), row=1:19))
# draw(lgd, x = unit(28, "cm"), y = unit(15, "cm"))

library(ggtree)
library(ggplot2)
library(ape)
library(magrittr)
tree <- read.tree(file = "/Users/alexwang/0data/0mango/SpeciesTreeAlignment.trimAl.phy_phyml_tree.txt")
#tree$tip.label
cladeA = c("yn33_prote", "Csp_protei", "hn42_prote", "Ctr_protei", "yn31_prote")
tree <- drop.tip(tree, cladeA)
tree <- phytools::midpoint.root(tree)

r <- ggtree(tree, branch.length = "none", layout = "circular") + scale_x_reverse(limits=c(100, 0)) + 
  geom_text(aes(label=node),color="red", size=3)
# ratate circular tree
r <- rotate_tree(r, 85)
# exchange two node
r <- flip(r, node1 = 41, node2 = 39)

scaleClade(r, node = 23, scale = 2.3) %>% 
  scaleClade(node = 36, scale = 1.8) %>% scaleClade( node = 33, scale = 1.2) %>% scaleClade(node = 32, scale = 2.0)%>% 
  scaleClade(node=35, scale=1.2)  %>% scaleClade(node = 33, scale=0.6) %>% scaleClade(node = 27, scale = 0.9) %>%
  scaleClade(node = 30, scale = 0.7) %>% scaleClade(node = 31, scale = 1.0) %>% scaleClade(node = 28, scale = 1.1) %>% 
  scaleClade(node = 29, scale = 1.2) %>% scaleClade(node = 26, scale = 1.4) %>% scaleClade(node = 37, scale = 1.4) %>% 
  scaleClade(node = 31, scale = 0.75) %>% scaleClade(node = 23, scale = 0.6) %>% scaleClade(node = 24, scale = 0.8) %>% 
  scaleClade(node = 41, scale=0.6) %>% scaleClade(node = 40, scale = 0.6) %>% scaleClade(node = 38, scale = 0.3) %>% 
  scaleClade(node = 39, scale = 2) %>% scaleClade(node = 25, scale = 1.1)

## calculate circos' scale factors
na <- unlist(lapply(strsplit(mini.fai$V1,split = "_"), function(x) x[1]))
na.uniq <- unique(na)
ggs <- aggregate(x = mini.gs, by = list(na), sum)
ggs.res <- ggs[order(match(ggs$Group.1,na)),]
scale.factors <- ggs.res$x / max(ggs.res$x)







