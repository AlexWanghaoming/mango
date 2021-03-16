library(edgeRun)
args <- commandArgs(trailingOnly = TRUE)
# strain <- "yn56"
strain <- args[1]
raw.reads <- read.table(paste0("/Users/alexwang/0data/0mango/transcriptome/",strain,"_featureCounts.tsv"), header = T, stringsAsFactors = FALSE, row.names = 1)[,-c(1,2,3,4,5)]

if(strain %in% c("yn56", "gd10")){
  y <- DGEList(counts = raw.reads, group = rep(c("coni", "hypha", "inf3d", "inf5d"), each=3))
}else if(strain == "gz15"){
  y <- DGEList(counts = raw.reads, group = rep(c("hypha","inf5d"), times=c(2,3)))
}else{
  y <- DGEList(counts = raw.reads, group = rep(c("hypha","inf5d"), each=3))
}
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, ]
logcpm <- cpm(y, prior.count=2, log=TRUE)

y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y, method = "TMM")
y <- estimateCommonDisp(y, verbose = TRUE)
y <- estimateTagwiseDisp(y)

# treatment <- "inf3d"
# case <- "hypha"

if(strain %in% c("gd10", "yn56")){
  # pairs = list(c("inf3d", "inf5d"),
  #              c("hypha", "inf5d"),
  #              c("hypha", "inf3d"),
  #              c("coni", "inf5d"),
  #              c("coni", "inf3d"))
  pairs = list(c("coni", "inf5d"),
               c("coni", "inf3d"))
}else{
  pairs = list(c("hypha", "inf5d"))
}

lapply(pairs, function(x){
  z <- UCexactTest(y, pair = x)
  de <- decideTestsDGE(z, adjust.method = "BH", p.value = 0.05, lfc = 1)
  tbl <- tibble::rownames_to_column(z$table, var = "gene")
  
  detags <- rownames(y)[as.logical(de)]
  uptags <- rownames(y)[as.logical(de>0)]
  downtags <- rownames(y)[as.logical(de<0)]
  
  ## output differetial expression genes
  outdir <- "/Users/alexwang/0data/0mango/transcriptome/DEG/"
  write.table(tbl,paste0(outdir, strain, "_", x[2], "_Vs_",x[1],".tbl"), quote = F, row.names = F)
  write.table(detags, paste0(outdir, strain,"_", x[2], "_Vs_",x[1],".geneID"), quote = F, col.names = F, row.names = F)
  write.table(uptags, paste0(outdir, strain,"_", x[2], "_Vs_",x[1],".up.geneID"), quote = F, col.names = F, row.names = F)
  write.table(downtags, paste0(outdir, strain,"_", x[2], "_Vs_",x[1],".down.geneID"), quote = F, col.names = F, row.names = F)
})


