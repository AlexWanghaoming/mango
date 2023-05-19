plot.new()
## normalizing data to [m,M]
normlize <- function(data,m,M){
  k <- (M-m)/(max(data)-min(data))
  re <- sapply(data, function(x) m+k*(x-min(data)))
  return(re)
}

# global settings
max_refseq_num <- 10
height <- 0.12
space <- 0.01
intvl <- height + space
gap=0.4
# cols <- c("wheat1",RColorBrewer::brewer.pal(6,"Reds"))
#-------------------------------------------------------------------------
# Add track Function
#-------------------------------------------------------------------------
addtrack <- function(track, galign_info_file, tag, min_ali_length=10000, ...) {
  # galign <- read.table("/Users/alexwang/0data/0mango/genome/mummer/gz15_gz15.bundle_telo2.txt", header = T)
  galign <- read.table(galign_info_file, header = T)
  
  galign$reorder <- sapply(galign$ref, FUN = function(x) which(dd$chr == x))
  galign <- galign[abs(galign$s1-galign$e1)>min_ali_length & abs(galign$s2-galign$e2)>min_ali_length,]
  # add alignment
  galign_length = dim(galign)[1]
  aa <- unique(galign$ali)
  colnum <- length(aa)
  co <- rainbow(colnum)
  for (i in 1:galign_length) {
    a <- galign[i,]
    reorder_ali_len <- abs(a$e1 -a$s1)
    ali_ali_len <- abs(a$e2 -a$s2)
    if (reorder_ali_len > min_ali_length & ali_ali_len > min_ali_length) {
      idx <- a$reorder
      if (idx <= max_refseq_num) {
        rect(a$s1, idx-height+track*intvl+idx*gap, a$e1, idx-space+track*intvl+idx*gap, border = NA, lwd = 0.4, col = co[which(aa==a$ali)])
        if(a$bundle_size == 1){
          arrows(x0 = a$s1+reorder_ali_len/2.5, y0 = idx+track*intvl-(height+space)/2+idx*gap, x1 = a$e1-reorder_ali_len/2.5, idx+track*intvl-(height+space)/2+idx*gap,length = 0.05,lwd = 1.2,col = "#2C434B")
        }
        else{
          arrows(x0 = a$e1-reorder_ali_len/2.5, y0 = idx+track*intvl-(height+space)/2+idx*gap, x1 = a$s1+reorder_ali_len/2.5, idx+track*intvl-(height+space)/2+idx*gap,length = 0.05,lwd = 1.2,col = "#2C434B")
        }
      }
    }
  }
  ##  add telo
  for (i in 1:galign_length) {
    a <- galign[i,]
    reorder_ali_len <- abs(a$e1 -a$s1)
    ali_ali_len <- abs(a$e2 -a$s2)
    if (reorder_ali_len > min_ali_length & ali_ali_len > min_ali_length) {
      idx <- a$reorder
      if (idx <= max_refseq_num) {
        points(a$telo1, idx+track * intvl-(height+space)/2+idx*gap, pch=16, col = 2, cex = 0.4)
        points(a$telo2, idx+track * intvl-(height+space)/2+idx*gap, pch=16, col = 2, cex = 0.4)
      }
    }
  }
  legend(5000000+400000*track, 14.5, legend = sort(aa), col = co[order(aa)], box.lwd = NA, title = tag, pch=15, cex=0.6)
}


#-------------------------------------------------------------------------
# Plot Tracks
#-------------------------------------------------------------------------

## genome track

genome_info_file <- "/Users/alexwang/0data/0mango/genome/mummer/genome/gz15.genome.txt"
repeat_file <- "/Users/alexwang/0data/0mango/genome/repeat/gz15.repeat.bed"

dd <- read.table(genome_info_file, header = T)
dd <- dd[order(dd$size,decreasing = T),]
refseq <- dd[1:max_refseq_num,]
repeat_loci <- read.table(repeat_file, header = T)
# repeat.LTR <- repeat_loci[repeat_loci$type == "LTR/Gypsy" | repeat_loci$type == "LTR/Copia",]
# repeat.oth <- repeat_loci[repeat_loci$type != "LTR/Gypsy" & repeat_loci$type != "LTR/Copia",]
pdf("gz15_col.pdf", width = 15, height = 8)
# prepare
par(mar = c(3,1.5,0,0), mgp=c(1.5,0.5,0), xpd = T)
max_x <- max(refseq$size) + 1000
max_y <- max_refseq_num + 5

plot(1,1, pch=NA, xlim = c(1, max_x), ylim = c(1, max_y), axes = F, xlab = "genomic position (kb)", ylab = "", cex.lab = 0.8)
axis(1, at = c(0, 1500000, 3000000, 4500000, 6000000, 7500000), 
     labels = c("0", "1500", "3000", "4500", "6000", "7500"), cex.axis = 0.6)

# genome layout
cov <- read.table("/Users/alexwang/0data/0mango/genome/coverage/gz15_bin.txt")
colnames(cov) <- c("chr","start","end","me","max","min")
cov$chr <- as.numeric(sub("gz15_","",cov$chr))
for (idx in seq(1, max_refseq_num)) {
  rect(1, idx-height+gap*idx, refseq$size[idx],  idx-space+gap*idx, col = cols[1], border = "gray", lwd = 0.4)    # 画参考基因组
  chr.num <- refseq$chr[idx]
  cov.tmp <- cov[cov$chr==chr.num,]
  rect(cov.tmp$start, idx-height+gap*idx-0.45, cov.tmp$end,idx-height+gap*idx-0.45 + normlize(cov.tmp$mi,0,0.45), col="#BEBADA",border = NA)
  if (sum(repeat_loci$ref == idx) > 0) {
    # rect(repeat.oth$s[repeat.oth$ref == refseq$chr[idx]], idx - height+gap*idx, repeat.oth$e[repeat.oth$ref == refseq$chr[idx]], idx - space+gap*idx, border = NA, col = 2)
    # rect(repeat.LTR$s[repeat.LTR$ref == refseq$chr[idx]], idx - height+gap*idx, repeat.LTR$e[repeat.LTR$ref == refseq$chr[idx]], idx - space+gap*idx, border = NA, col = 3)
    rect(repeat_loci$s[repeat_loci$ref == refseq$chr[idx]], idx - height+gap*idx, repeat_loci$e[repeat_loci$ref == refseq$chr[idx]], idx - space+gap*idx, border = NA, col = 3)
  }
  text(-200000, idx + -0.5 * intvl+gap*idx,  paste("gz15", refseq$chr[idx], sep = "_"), adj = 0, cex = 0.5)
  text(-200000, idx +  0.5 * intvl+gap*idx,  "yn55", adj = 0, cex = 0.5)
  text(-200000, idx +  1.5 * intvl+gap*idx,  "fj11", adj = 0, cex = 0.5)
  text(-200000, idx +  2.5 * intvl+gap*idx,  "qz", adj = 0, cex = 0.5)
  text(-200000, idx +  3.5 * intvl+gap*idx,  "yn56", adj = 0, cex = 0.5)
  text(-200000, idx +  4.5 * intvl+gap*idx,  "gd10", adj = 0, cex = 0.5)
  text(-200000, idx +  5.5 * intvl+gap*idx,  "hn47", adj = 0, cex = 0.5)
}
addtrack(1, "/Users/alexwang/0data/0mango/genome/mummer/gz15_yn55.bundle_telo2.txt","yn55")
addtrack(2, "/Users/alexwang/0data/0mango/genome/mummer/gz15_fj11.bundle_telo2.txt", "fj11")
addtrack(3, "/Users/alexwang/0data/0mango/genome/mummer/gz15_qz.bundle_telo2.txt", "qz")
addtrack(4, "/Users/alexwang/0data/0mango/genome/mummer/gz15_yn56.bundle_telo2.txt", "yn56")
addtrack(5, "/Users/alexwang/0data/0mango/genome/mummer/gz15_gd10.bundle_telo2.txt", "gd10")
addtrack(6, "/Users/alexwang/0data/0mango/genome/mummer/gz15_hn47.bundle_telo2.txt", "hn47")
dev.off()
