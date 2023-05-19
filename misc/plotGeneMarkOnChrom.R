library(karyoploteR)
library(GenomicFeatures)
library(magrittr)

plotGeneOnChrom <- function(name, gs_file, GTF_file, gene_list, repeat_file=NULL){
  gs <- read.table(gs_file, header = T)
  gs <- gs[order(gs$size),]
  gs.gr <- GRanges(seqnames = gs$chr, ranges = IRanges(start=1, end = gs$size), col="red")
  ## bundle regions
  pp = getDefaultPlotParams(plot.type = 1)
  pp$leftmargin = 0.05
  pp$topmargin = 80
  pp$bottommargin = 15
  pp$ideogramheight = 20
  pp$data1inmargin = 10
  pp$data1outmargin = 0
  
  #### read GTF
  txdb <- makeTxDbFromGFF(file = GTF_file)
  all.gene <- genes(txdb)
  kp = plotKaryotype(genome = gs.gr, plot.params = pp, cex=0.6)
  gene.gr <- GRanges(seqnames = gsub(paste0(name,"_"),"",seqnames(all.gene)), 
                     ranges = IRanges(start = start(all.gene), end = end(all.gene)),
                     strand = strand(all.gene),
                     gene.id = all.gene$gene_id)
  ######################### Effectors
  target.id <- read.table(gene_list)[,1]
  target.gr <- all.gene[all.gene$gene_id %in% target.id]
  target.gr <- GRanges(seqnames = gsub(paste0(name,"_"),"",seqnames(target.gr)), 
                         ranges = IRanges(start = start(target.gr), end = end(target.gr)),
                         strand = strand(target.gr),
                         gene.id = target.gr$gene_id)
  kpPlotRegions(kp, data = target.gr,
                r0 = 0, r1 = 0.8, col = "darkblue", avoid.overlapping = TRUE, lwd=1)
  #### plot repeat position
  if(!is.null(repeat_file)){
    cat("loadind RepeatMasker file ...", "\n")
    repe <- read.table(repeat_file)[,c(1,4,5,7)]
    colnames(repe) <- c("ref", "s", "e", "strand")
    repe$ref <- gsub(pattern = paste0(name,"_"), replacement = "",repe$ref)
    kp <- kpPlotRegions(kp, data = GRanges(seqnames = repe$ref,
                                           ranges = IRanges(start = repe$s, end = repe$e)),
                        r0 = 0, r1 = 0.8, col = "red", border = NA, avoid.overlapping = FALSE,lwd=0.5)
  }
}
strain = "Chi"
gs_file <- paste0("/Users/alexwang/0data/0mango/genome/mummer/genome/",strain,".genome.txt")
GTF_file <- paste0("/Users/alexwang/0data/0mango/GTF/",strain,"_EVM.gtf")
gene_list <- paste0("/Users/alexwang/0data/0mango/inf_associated_gene/secretome/",strain,"_effector.geneID")
repeat_file <- paste0("/Users/alexwang/0data/0mango/LTR/LTR_repeatmasker/",strain,".fasta.out.gff")
plotGeneOnChrom(name = strain, gs_file, GTF_file, gene_list)


