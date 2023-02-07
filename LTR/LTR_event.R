###################################### LTR from LTR_harvest deprecated ！
library(ape)
library(Biostrings)
library(GenomicFeatures)
library(Rsamtools)
library(systemPipeR)
yn56.genome <- "/Users/alexwang/0data/0mango/genome/yn56.fasta"
yn55.genome <- "/Users/alexwang/0data/0mango/genome/yn55.fasta"
qz.genome <- "/Users/alexwang/0data/0mango/genome/qz.fasta"
fj11.genome <- "/Users/alexwang/0data/0mango/genome/fj11.fasta"
gd10.genome <- "/Users/alexwang/0data/0mango/genome/gd10.fasta"
gz15.genome <- "/Users/alexwang/0data/0mango/genome/gz15.fasta"
hn47.genome <- "/Users/alexwang/0data/0mango/genome/hn47.fasta"
c.acutatum.genome <- "/Users/alexwang/0data/0mango/oth_colletotrichum/c.acutatum/c.acutatum.fa"
c.destructivum.genome <- "/Users/alexwang/0data/0mango/oth_colletotrichum/c.destructivum/c.destructivum.fa"
c.gloeosporioides.genome <- "/Users/alexwang/0data/0mango/oth_colletotrichum/c.gloeosporioides_SMCG1C/c.gloeosporioides.fa"
c.higginsianum.genome <- "/Users/alexwang/0data/0mango/oth_colletotrichum/c.higginsianum_IMI349063/c.higginsianum.fa"
c.lentis.genome <- "/Users/alexwang/0data/0mango/oth_colletotrichum/c.lentis/c.lentis.fa"
c.sp.JS367.genome <- "/Users/alexwang/0data/0mango/oth_colletotrichum/c.sp.JS-367/c.sp.JS-367.fa"
c.truncatum.genome <- "/Users/alexwang/0data/0mango/oth_colletotrichum/c.truncatum/c.truncatum.fa"

yn56.ltr <- "/Users/alexwang/0data/0mango/LTR/yn56.harvest.scn"
yn55.ltr <- "/Users/alexwang/0data/0mango/LTR/yn55.harvest.scn"
qz.ltr <- "/Users/alexwang/0data/0mango/LTR/qz.harvest.scn"
fj11.ltr <- "/Users/alexwang/0data/0mango/LTR/fj11.harvest.scn"
gd10.ltr <- "/Users/alexwang/0data/0mango/LTR/gd10.harvest.scn"
gz15.ltr <- "/Users/alexwang/0data/0mango/LTR/gz15.harvest.scn"
hn47.ltr <- "/Users/alexwang/0data/0mango/LTR/hn47.harvest.scn"
c.acutatum.ltr <- "/Users/alexwang/0data/0mango/LTR/c.acutatum.harvest.scn"
c.destructivum.ltr <- "/Users/alexwang/0data/0mango/LTR/c.destructivum.harvest.scn"
c.gloeosporioides.ltr <- "/Users/alexwang/0data/0mango/LTR/c.gloeosporioides_SMCG1C.harvest.scn"
c.higginsianum.ltr <- "/Users/alexwang/0data/0mango/LTR/C.higginsianum_IMI349063.harvest.scn"
c.lentis.ltr <- "/Users/alexwang/0data/0mango/LTR/c.lentis.harvest.scn"
c.sp.JS367.ltr <- "/Users/alexwang/0data/0mango/LTR/c.sp.JS-367.harvest.scn"
c.truncatum.ltr <- "/Users/alexwang/0data/0mango/LTR/c.truncatum.harvest.scn"

LTR_insert <- function(strain, genome, ltr){    #### result from LTR_harvest
  LTR_harvest.tab <- read.table(ltr)
  colnames(LTR_harvest.tab) <- c("s_ret", "e_ret", "l_ret", "s_lLTR", "e_lLTR", "l_lLTR", 
                                 "s_rLTR", "e_rLTR", "l_rLTR", "sim", "seq_id")
  chrom <- paste0(strain,"_",LTR_harvest.tab$seq_id+1)
  LTR_harvest.tab$seq_id <- chrom
  
  genome.fasta <- Rsamtools::FaFile(genome)
  l.ltr <- GRanges(seqnames = LTR_harvest.tab$seq_id, 
                   ranges = IRanges(start = LTR_harvest.tab$s_lLTR, end = LTR_harvest.tab$e_lLTR))
  r.ltr <- GRanges(seqnames = LTR_harvest.tab$seq_id, 
                   ranges = IRanges(start = LTR_harvest.tab$s_rLTR, end = LTR_harvest.tab$e_rLTR))
  inter.gr <- GRanges(seqnames = LTR_harvest.tab$seq_id, 
                      ranges = IRanges(start = LTR_harvest.tab$e_lLTR, end = LTR_harvest.tab$s_rLTR))
  # showMethods("getSeq")
  # methods(class="FaFile")
  l.ltr.seq <- getSeq(genome.fasta, l.ltr)
  r.ltr.seq <- getSeq(genome.fasta, r.ltr)
  ## predict ORF of inter region
  inter.seq <- getSeq(genome.fasta, inter.gr)
  
  
  pair_align <- pairwiseAlignment(l.ltr.seq, r.ltr.seq)
  seq1 <- as.character(pattern(pair_align))
  seq2 <- as.character(subject(pair_align))
  res <- list()
  # res2 <- list()
  for (i in 1:length(seq1)) {
    x <- list(strsplit(seq1, "")[[i]], strsplit(seq2, "")[[i]])
    xx <- as.DNAbin(x)
    d <- dist.dna(xx)[1]
    # res2[[i]] <- d
    # fungal substitution rate  nucleotides per site per year
    r <- 1.05*10e-9 
    divergent_time <- d/(2*r)
    res[[i]] <- divergent_time/1000000
  }
  rr1 <- unlist(res)
  # rr2 <- unlist(res2)
  # return(list(rr1,rr2))
  return(rr1)
}
LTR_orf <- function(strain, genome, ltr){
  LTR_harvest.tab <- read.table(ltr)
  colnames(LTR_harvest.tab) <- c("s_ret", "e_ret", "l_ret", "s_lLTR", "e_lLTR", "l_lLTR", 
                                 "s_rLTR", "e_rLTR", "l_rLTR", "sim", "seq_id")
  chrom <- paste0(strain,"_",LTR_harvest.tab$seq_id+1)
  LTR_harvest.tab$seq_id <- chrom
  genome.fasta <- Rsamtools::FaFile(genome)
  inter.gr <- GRanges(seqnames = LTR_harvest.tab$seq_id, 
                      ranges = IRanges(start = LTR_harvest.tab$e_lLTR, end = LTR_harvest.tab$s_rLTR))
  ## predict ORF of inter region
  inter.seq <- getSeq(genome.fasta, inter.gr)
  names(inter.seq) <- paste0(names(inter.seq),"_LTR",1:length(inter.seq))
  orf <- predORF(inter.seq, n = 1, type = "gr", mode = "orf", strand = "both")
  return(orf)
}

## Predict longest ORF for sense strand in each query sequence

# orf <- LTR_orf("yn55", genome = yn55.genome, ltr = yn55.ltr)
# sort(getSeq(inter.seq, orf))


fj11.res <- density(LTR_insert("fj11", genome = fj11.genome, ltr = fj11.ltr))
yn55.res <- density(LTR_insert("yn55", genome = yn55.genome, ltr = yn55.ltr))
yn56.res <- density(LTR_insert("yn56", genome = yn56.genome, ltr = yn56.ltr))
qz.res <- density(LTR_insert("qz", genome = qz.genome, ltr = qz.ltr))
gd10.res <- density(LTR_insert("gd10", genome = gd10.genome, ltr = gd10.ltr))
gz15.res <- density(LTR_insert("gz15", genome = gz15.genome, ltr = gz15.ltr))
hn47.res <- density(LTR_insert("hn47", genome = hn47.genome, ltr = hn47.ltr))

c.acutatum.res <- density(LTR_insert("c.acutatum", genome = c.acutatum.genome, ltr = c.acutatum.ltr))
c.destructivum.res <- density(LTR_insert("c.destructivum", genome = c.destructivum.genome, ltr = c.destructivum.ltr))
c.higginsianum.res <- density(LTR_insert("c.higginsianum", genome = c.higginsianum.genome, ltr = c.higginsianum.ltr))
c.gloeosporioides.res <- density(LTR_insert("c.gloeosporioides", genome = c.gloeosporioides.genome, ltr = c.gloeosporioides.ltr))
c.lentis.res <- density(LTR_insert("c.lentis", genome = c.lentis.genome, ltr = c.lentis.ltr))
c.sp.JS367.res <- density(LTR_insert("c.sp.JS-367", genome = c.sp.JS367.genome, ltr = c.sp.JS367.ltr))
c.truncatum.res <- density(LTR_insert("c.truncatum", genome = c.truncatum.genome, ltr = c.truncatum.ltr))

options(scipen = 200)
plot(fj11.res, main = "Time of insertion of LTR",col="red", 
     lwd=1.5, xlab="MYA",ylim=c(0,1.25), xlim=c(0,10)
     # xlim=c(0,9), ylim=c(0, 0.3)
)
# lines(yn55.res, col="orange", lwd=1.5)
# lines(yn56.res, col="purple", lwd=1.5)
# lines(qz.res, col="darkblue", lwd=1.5)
# lines(gd10.res, col="green", lwd=1.5)
# lines(gz15.res, col="lightblue", lwd=1.5)
# lines(hn47.res, col="darkgreen", lwd=1.5)
lines(c.acutatum.res, col="orange", lwd=1.5)
lines(c.destructivum.res, col="purple", lwd=1.5)
lines(c.higginsianum.res, col="black", lwd=1.5)
lines(c.gloeosporioides.res, col="darkblue", lwd=1.5)
lines(c.lentis.res, col="gold", lwd=1.5)
lines(c.sp.JS367.res, col="green", lwd=1.5)
lines(c.truncatum.res, col="darkgreen", lwd=1.5)
legend("topright", legend = c("fj11", "yn55", "yn56", "qz", "gd10", "gz15", "hn47"),
       col = c("red", "orange", "purple", "darkblue", "green", "lightblue", "darkgreen"),
       pch = "l",bty = "n")


###################### plot LTR Phylo tree
insert_time.fj11 <- LTR_insert("fj11", genome = fj11.genome, ltr = fj11.ltr)
insert_time.yn55 <-LTR_insert("yn55", genome = yn55.genome, ltr = yn55.ltr)
insert_time.yn56 <- LTR_insert("yn56", genome = yn56.genome, ltr = yn56.ltr)
insert_time.qz <- LTR_insert("qz", genome = qz.genome, ltr = qz.ltr)
insert_time.gd10 <- LTR_insert("gd10", genome = gd10.genome, ltr = gd10.ltr)
insert_time.gz15 <- LTR_insert("gz15", genome = gz15.genome, ltr = gz15.ltr)
insert_time.hn47 <- LTR_insert("hn47", genome = hn47.genome, ltr = hn47.ltr)
c(paste0("fj11_LTR", 1:length(insert_time.fj11)),
  paste0("yn55_LTR", 1:length(insert_time.yn55)), 
  paste0("yn56_LTR", 1:length(insert_time.yn56)), 
  paste0("qz_LTR", 1:length(insert_time.qz)), 
  paste0("gd10_LTR", 1:length(insert_time.gd10)),
  paste0("gz15_LTR", 1:length(insert_time.gz15)),
  paste0("hn47_LTR", 1:length(insert_time.hn47)))

insertTime.df <- data.frame(label = c(paste0("fj11_LTR", 1:length(insert_time.fj11)),
                                      paste0("yn55_LTR", 1:length(insert_time.yn55)), 
                                      paste0("yn56_LTR", 1:length(insert_time.yn56)), 
                                      paste0("qz_LTR", 1:length(insert_time.qz)), 
                                      paste0("gd10_LTR", 1:length(insert_time.gd10)),
                                      paste0("gz15_LTR", 1:length(insert_time.gz15)),
                                      paste0("hn47_LTR", 1:length(insert_time.hn47))), 
                            insert_time = c(insert_time.fj11, insert_time.yn55, insert_time.yn56, insert_time.qz, insert_time.gd10, insert_time.gz15, insert_time.hn47))

tree <- read.tree("/Users/alexwang/0data/0mango/genome/LTR/msa/tree")
tree <- full_join(tree, insertTime.df, by="label")
ggtree(tree, layout = "circular", lwd=0.2, aes(color=insert_time))  + 
  scale_color_gradientn(colours=c("red","green",'blue')) + 
  geom_tiplab2(size=0.6) + 
  geom_text(aes(label=node), size=0.5)

split(insertTime.df$label, insertTime.df$insert_time)
groupOTU(tree, insertTime.df)
#####################################################################################################
#####################################################################################################
############# LTR distance result from LTR_retriever !
library(ape)
library(Biostrings)
library(GenomicFeatures)
library(Rsamtools)
library(systemPipeR)
library(ggplot2)

ltr_path <- "/Users/alexwang/0data/0mango/LTR/ltr_retriever/"
genome_path <- "/Users/alexwang/0data/0mango/genome/"

LTR_insert <- function(ltr, genome,type="all"){
  ltr.df <- read.table(ltr)
  if(type=="all"){
    ltr.df <- ltr.df
  }else{
    ltr.df <- ltr.df[ltr.df$V10 == type,]
  }
  a1 <- strsplit(ltr.df$V1, split="[.:]", perl=T)
  chrom <- unlist(lapply(a1, function(x) x[1]))
  s_lLTR <- as.numeric(unlist(lapply(a1, function(x) x[2])))
  e_rLTR <- as.numeric(unlist(lapply(a1, function(x) x[4])))
  
  a2 <- strsplit(ltr.df$V7, split="[.:]", perl=T)
  e_lLTR <- as.numeric(unlist(lapply(a2, function(x) x[2]))) - 1
  s_rLTR <- as.numeric(unlist(lapply(a2, function(x) x[4]))) + 1
  
  genome.fasta <- Rsamtools::FaFile(genome)
  l.ltr <- GRanges(seqnames = chrom, 
                   ranges = IRanges(start = s_lLTR, end = e_lLTR))
  r.ltr <- GRanges(seqnames = chrom, 
                   ranges = IRanges(start = s_rLTR, end = e_rLTR))
  
  l.ltr.seq <- getSeq(genome.fasta, l.ltr)
  r.ltr.seq <- getSeq(genome.fasta, r.ltr)
  
  pair_align <- pairwiseAlignment(l.ltr.seq, r.ltr.seq)
  seq1 <- as.character(pattern(pair_align))
  seq2 <- as.character(subject(pair_align))
  res <- list()
  for (i in 1:length(seq1)) {
    x <- list(strsplit(seq1, "")[[i]], strsplit(seq2, "")[[i]])
    xx <- as.DNAbin(x)
    d <- dist.dna(xx)[1]
    # res2[[i]] <- d
    # fungal substitution rate  nucleotides per site per year
    r <- 1.05*10e-9 
    divergent_time <- d/(2*r)
    res[[i]] <- divergent_time/1000000
  }
  rr1 <- unlist(res)
  return(rr1)
}
# read.table(ltr)[c(60,84),]
LTR_insert(ltr = "/Users/alexwang/0data/0mango/LTR/ltr_retriever/Cfr1.fasta.pass.list", 
           genome ="/Users/alexwang/0data/0mango/genome/Cfr1.fasta" )

strains=c("fj11", "yn55", "yn56", "gd10", "hn32","yn32","qz", "hn47","gz14","Cfr1","Cfr2", "gz15", 
          "Cgl", "plg3", "hn23", "Ctr", "yn31", "Csp", "hn42", "Cac", "yn51", "yn33", "Cde", "Chi", "Cle")
all.list <- list()

for (i in 1:length(strains)) {
  aa = LTR_insert(ltr=paste0(ltr_path, strains[i], ".fasta.pass.list"),
                             genome=paste0(genome_path, strains[i], ".fasta"), type="all")
  all.list[[strains[i]]] = aa
}

LTR.n <-  as.numeric(lapply(all.list, function(x) length(x)))
insertAge.df <- data.frame(age=unlist(all.list), 
                           strain=factor(rep(strains, times = LTR.n),levels = strains)
)
mean(insertAge.df[insertAge.df$strain %in% c("Cfr1", "Cfr2", "gz14", "qz", "hn47"), 1])
## volin-boxplot plot
ggplot(insertAge.df, aes(x=strain,y=age)) + geom_violin(scale = "width") + 
  geom_boxplot(fill="orange",width=0.25) + xlab("") + scale_x_discrete(labels=paste0(strains, "\nn=",LTR.n)) + 
  ylab("Insertion age (million years)") + theme_classic()
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
# LTR_insert("/Users/alexwang/0data/0mango/LTR/ltr_retriever/yn56.fasta.pass.list", 
#            genome = "/Users/alexwang/0data/0mango/genome/yn56.fasta", "Gypsy")

## densityplot
ydj.ltr.res <- density(LTR_insert(ltr, genome))
options(scipen = 200)
plot(ydj.ltr.res, main = "Time of insertion of LTR",col="purple", 
     lwd=1.5, xlab="Age (MYA)",ylim=c(0,4.2), xlim=c(0,3)
     # xlim=c(0,9), ylim=c(0, 0.3)
)
##########################################################################################################
############################## Transposon remains, RIP index, RIP exists
library(GenomicFeatures)
library(Rsamtools)
getRIPgc <- function(strain){
  LTR.masked <- read.table(paste0("/Users/alexwang/0data/0mango/LTR/LTR_repeatmasker/",strain,".fasta.out.gff"),fill=T)
  LTR.masked <- LTR.masked[LTR.masked$V5-LTR.masked$V4 > 400, ]
  LTR.masked.gr <- GRanges(seqnames = LTR.masked$V1, 
                           ranges = IRanges(start=as.numeric(LTR.masked$V4), 
                                            end=as.numeric(LTR.masked$V5)),
                           strand = LTR.masked$V7)
  
  LTR.intact <- read.table(paste0("/Users/alexwang/0data/0mango/LTR/ltr_retriever/",strain,".fasta.pass.list.gff3"))
  LTR.intact <- LTR.intact[LTR.intact$V3 == "repeat_region",]
  LTR.intact.gr <- GRanges(seqnames = LTR.intact$V1,
                           ranges = IRanges(start = LTR.intact$V4,
                                            end = LTR.intact$V5),
                           strand = gsub(pattern = "\\?","*",LTR.intact$V7))
  #### keep only the ranges in LTR.masked that do not overlap intact LTR
  LTR.remains <- subsetByOverlaps(LTR.masked.gr, LTR.intact.gr, invert = TRUE)
  #### calculate GC content
  genome.fasta <- FaFile(paste0("/Users/alexwang/0data/0mango/genome/",strain,".fasta"))
  seqs.remains <- getSeq(genome.fasta, LTR.remains)
  sum(letterFrequency(seqs.remains, "GC")) / sum(letterFrequency(seqs.remains, "ATCG"))
  seq.remains.gc <- sort(letterFrequency(seqs.remains, "GC") / letterFrequency(seqs.remains, "ATCG"))
  
  seqs.intact <- getSeq(genome.fasta, LTR.intact.gr)
  sum(letterFrequency(seqs.intact, "GC")) / sum(letterFrequency(seqs.intact, "ATCG"))
  seqs.intact.gc <-  sort(letterFrequency(seqs.intact, "GC") / letterFrequency(seqs.intact, "ATCG"))
  return(list(remains = seq.remains.gc, 
              intact =seqs.intact.gc) )
}

r <- getRIPgc("cfr2")
t.test(r[[1]], r[[2]], alternative = "less")

marksig <- function(x1, x2, y, p, seg_len, offset) {
  lines(c(x1, x2), c(y, y), lwd = 2)
  lines(c(x1, x1), c(y, y - seg_len), lwd = 2)
  lines(c(x2, x2), c(y, y - seg_len), lwd = 2)
  text(x1+(x2-x1)/2, y + offset, p, cex = 1.5)
}
add_axis <- function(x1, x2, y, label, tick_len, offset) {
  lines(c(x1, x2), c(y, y), lwd = 2)
  lines(c(x1, x1), c(y, y+tick_len), lwd = 2)
  lines(c(x2, x2), c(y, y+tick_len), lwd = 2)
  text(x1+(x2-x1)/2, y-offset, labels = label, cex = 1)
}
myboxplot <- function(strain, name, sig, ylab = ""){
  r <- getRIPgc(strain)
  x1 = 1
  x2 = 2
  tick_len = 0.01
  offset = 0.01
  boxplot(r, 
          col = ggsci::pal_lancet('lanonc',alpha = 0.6)(2),
          ylim=c(0.1, 0.7),
          axes = F,
          outline = F,
          frame = F,
          xlab = "", ylab = "")
  marksig(x1, x2, 0.7, sig,tick_len, offset)
  add_axis(x1, x2, 0.05, name, tick_len, offset)
}
pdf("LTR_GC.pdf", width = 13, height = 5)
par(mar = c(5,2,1,1), mfrow=c(1, 10), xpd = T)
myboxplot(strain = "yn55", name = "YN55", "***")
axis(2, lwd = 2)
myboxplot(strain = "fj11", name = "FJ11", "***")
myboxplot(strain = "qz", name = "QZ", "***")
myboxplot(strain = "hn47", name = "HN47", "***")
legend(1, 0.3, 
       legend = c("LTR wreckage", "Intact LTR"), 
       fill = ggsci::pal_lancet('lanonc',alpha = 0.6)(2), 
       bty = "n",
       border = NA,
       cex = 1.5)
myboxplot(strain = "gz14", name = "GZ14", "***")
myboxplot(strain = "cfr1", name = "Cfr1", "***")
myboxplot(strain = "cfr2", name = "Cfr2", "***")
myboxplot(strain = "yn56", name = "YN56", "ns")
myboxplot(strain = "gd10", name = "GD10", "ns")
myboxplot(strain = "gz15", name = "GZ15", "*")

dev.off()






slt <- stringi::stri_split(LTR.masked$V9, regex = "[.:=_]")
subject.gr <- GRanges(seqnames = unlist(lapply(slt, function(x) paste0(x[8],"_",x[9]))),
                      ranges = IRanges(start = as.numeric(unlist(lapply(slt, function(x) x[10]))),
                                       end = as.numeric(unlist(lapply(slt, function(x) x[12])))),
                      strand = "+")
subject.gr2 <- GRanges(seqnames = unlist(lapply(slt, function(x) paste0(x[8],"_",x[9]))),
                      ranges = IRanges(start = as.numeric(unlist(lapply(slt, function(x) x[10]))),
                                       end = as.numeric(unlist(lapply(slt, function(x) x[12])))),
                      strand = "-")

pairwiseAlignment(pattern = getSeq(genome.fasta, LTR.masked.gr),
                  subject = getSeq(genome.fasta, subject.gr))


pa.pos <- pairwiseAlignment(pattern = getSeq(genome.fasta, LTR.masked.gr),
                  subject = getSeq(genome.fasta, subject.gr))
pa.neg <- pairwiseAlignment(pattern = getSeq(genome.fasta, LTR.masked.gr),
                           subject = getSeq(genome.fasta, subject.gr2))
ch <- function(x){
  if(score(pa.pos[x]) >= score(pa.neg[x])){
    mismatch.table <- mismatchTable(pa.pos[x])
  }else{
    mismatch.table <- mismatchTable(pa.neg[x])
  }
  CT_mut.rat <- nrow(mismatch.table[mismatch.table$SubjectSubstring == "C" & mismatch.table$PatternSubstring == "T", ]) / nrow(mismatch.table)
  return(CT_mut.rat)
}
sapply(1:length(pa.pos), ch)




aa = mismatchTable(pa.neg)
aa$PatternSubstring
aa$SubjectSubstring
mismatchSummary(pa.pos[[1]])
pa.pos
pairwiseAlignment(pa.pos[1])

UseMethod("PairwiseAlignments")
methods(pa.pos)
class(pa.pos)
Class
pairwiseAlignment()
?pairwiseAlignment




