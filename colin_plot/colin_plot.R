# global settings
max_refseq_num <- 10
height <- 0.12
space <- 0.01
intvl <- height + space
cols <- c("wheat1", "#D4AF37", "#40BCD8", "#39A9DB", "#FBCAEF", "#DC98D9", "#F865B0")

#-------------------------------------------------------------------------
# Add track Function
#-------------------------------------------------------------------------
addtrack <- function(track, galign_info_file, track_col,  min_ali_length=10000, ...) {
  
  galign <- read.table(galign_info_file, header = T)
  galign$reorder <- sapply(galign$ref, FUN = function(x) which(dd$chr == x))
  
  # add alignment
  galign_length = dim(galign)[1]
  for (i in 1:galign_length) {
    a <- galign[i,]
    reorder_ali_len <- abs(a$e1 -a$s1)
    ali_ali_len <- abs(a$e2 -a$s2)
    if (reorder_ali_len > min_ali_length & ali_ali_len > min_ali_length) {
      if (a$reorder <= max_refseq_num) {
        rect(a$s1, a$reorder - height + track * intvl, a$e1, a$reorder - space + track * intvl, border = "gray", lwd = 0.4, col = track_col)
        if(a$bundle_size == 1){
          arrows(x0 = a$s1+reorder_ali_len/2.5, y0 = a$reorder+track*intvl-(height+space)/2, x1 = a$e1-reorder_ali_len/2.5, a$reorder+track*intvl-(height+space)/2,length = 0.05,lwd = 1.2,col = "#2C434B")
          }
        else{
            arrows(x0 = a$e1-reorder_ali_len/2.5, y0 = a$reorder+track*intvl-(height+space)/2, x1 = a$s1+reorder_ali_len/2.5, a$reorder+track*intvl-(height+space)/2,length = 0.05,lwd = 1.2,col = "#2C434B")
          }
        text((a$s1+a$e1)/2, a$reorder - space + (track - 0.5) * intvl, a$ali, cex = 0.5)
         }
    }
  }
  ##  add telo
  for (i in 1:galign_length) {
    a <- galign[i,]
    reorder_ali_len <- abs(a$e1 -a$s1)
    ali_ali_len <- abs(a$e2 -a$s2)
    if (reorder_ali_len > min_ali_length & ali_ali_len > min_ali_length) {
      if (a$reorder <= max_refseq_num) {
        points(a$telo1, a$reorder + track * intvl-(height+space)/2, pch=16, col = 2, cex = 0.4)
        points(a$telo2, a$reorder + track * intvl-(height+space)/2, pch=16, col = 2, cex = 0.4)
      }
    }
  }
  ##  add scaffold ends
  # for (i in 1:galign_length) {
  #   a <- galign[i,]
  #   ref_ali_len <- abs(a$e1 -a$s1)
  #   ali_ali_len <- abs(a$e2 -a$s2)
  #   if (ref_ali_len > min_ali_length & ali_ali_len > min_ali_length) {
  #     if (a$ref <= max_refseq_num) {
  #       points(a$s3, a$ref + (track - 0.5) * intvl, pch=1, col = "gray", cex = 0.4)
  #       points(a$e3, a$ref + (track - 0.5) * intvl, pch=1, col = "gray", cex = 0.4)
  #     }
  #   }
  # }
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
repeat.LTR <- repeat_loci[repeat_loci$type == "LTR/Gypsy" | repeat_loci$type == "LTR/Copia",]
repeat.oth <- repeat_loci[repeat_loci$type != "LTR/Gypsy" & repeat_loci$type != "LTR/Copia",]
# pdf("gz15.pdf", width = 7, height = 8)
# prepare
par(mar = c(3,1.5,0,0), mgp=c(1.5,0.5,0), xpd = T)
max_x <- max(refseq$size) + 1000
max_y <- max_refseq_num + 0.5

plot(1,1, pch=NA, xlim = c(1, max_x), ylim = c(1, max_y), axes = F, xlab = "genomic position (kb)", ylab = "", cex.lab = 0.8)
axis(1, at = c(0, 1500000, 3000000, 4500000, 6000000, 7500000), 
     labels = c("0", "1500", "3000", "4500", "6000", "7500"), cex.axis = 0.6)

# genome layout
for (idx in seq(1, max_refseq_num)) {
  rect(1, idx - height , refseq$size[idx],  idx - space, col = cols[1], border = "gray", lwd = 0.4)
  if (sum(repeat_loci$ref == idx) > 0) {
    rect(repeat.oth$s[repeat.oth$ref == refseq$chr[idx]], idx - height, repeat.oth$e[repeat.oth$ref == refseq$chr[idx]], idx - space, border = NA, col = 2)
    rect(repeat.LTR$s[repeat.LTR$ref == refseq$chr[idx]], idx - height, repeat.LTR$e[repeat.LTR$ref == refseq$chr[idx]], idx - space, border = NA, col = 3)
    }
  # track names
  text(-200000, idx + -0.5 * intvl,  paste("gz15", refseq$chr[idx], sep = "_"), adj = 0, cex = 0.6)
  text(-200000, idx +  0.5 * intvl,  "hn47", adj = 0, cex = 0.6)
  text(-200000, idx +  1.5 * intvl,  "gd10", adj = 0, cex = 0.6)
  text(-200000, idx +  2.5 * intvl,  "yn56", adj = 0, cex = 0.6)
  text(-200000, idx +  3.5 * intvl,  "yn55", adj = 0, cex = 0.6)
  text(-200000, idx +  4.5 * intvl,  "qz", adj = 0, cex = 0.6)
  text(-200000, idx +  5.5 * intvl,  "fj11", adj = 0, cex = 0.6)
}
addtrack(1, "/Users/alexwang/0data/0mango/genome/mummer/gz15_hn47.bundle_telo.txt", cols[2])
addtrack(2, "/Users/alexwang/0data/0mango/genome/mummer/gz15_gd10.bundle_telo.txt", cols[3])
addtrack(3, "/Users/alexwang/0data/0mango/genome/mummer/gz15_yn56.bundle_telo.txt", cols[4])
addtrack(4, "/Users/alexwang/0data/0mango/genome/mummer/gz15_yn55.bundle_telo.txt", cols[5])
addtrack(5, "/Users/alexwang/0data/0mango/genome/mummer/gz15_qz.bundle_telo.txt", cols[6])
addtrack(6, "/Users/alexwang/0data/0mango/genome/mummer/gz15_fj11.bundle_telo.txt", cols[6])

legend(max_x * 0.9, max_y, c("Telomere", "Scaffold End"), pch = c(16, 1), col = c(2, "gray"), cex = 0.4, box.lwd = NA)

# dev.off()

