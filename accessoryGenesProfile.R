# codon preference analysis: CAI

## shell script
## cd /Users/alexwang/0data/0mango/accessory/codon_analysis
## awk 'NR==FNR{a[$1]=$1}NR>FNR{RS=">";FS="\n";if ($1 in a) {printf ">%s",$0}}' ../HGT/putative_HGT.geneID all.cds > putative_HGT.cds.fasta
### source activate emboss
### Use all chromosomes of Colletotrichum build codon usage table
# cusp -sequence all.cds -outfile Colletotrichum.cusp
### HGT genes
##cai -seqall putative_HGT.cds.fasta -cfile Colletotrichum.cusp -outfile HGT.cai
# core chromosomes
# cai -seqall all.cds -cfile Colletotrichum.cusp -outfile Colletotrichum.cai

core.cai <- read.table("/Users/alexwang/0data/0mango/accessory/codon_analysis/Colletotrichum.cai")
HGT.cai <- read.table("/Users/alexwang/0data/0mango/accessory/codon_analysis/HGT.cai")
gz23.cai <- read.table("/Users/alexwang/0data/0mango/accessory/codon_analysis/gz23_HGT.cai")

core.density <- density(core.cai$V4)
HGT.density <- density(HGT.cai$V4)

boxplot(core.cai$V4,HGT.cai$V4, gz23.cai$V4)
t.test(core.cai$V4, HGT.cai$V4)
plot(core.density, lwd=2, col = "gray", ylim=c(0, 10), xlab="Codon adaptation index (CAI)", 
     ylab="Percentage of genes (%)", main="")
lines(HGT.density, lwd=1.5, col="orange")

##########################################################
library(magrittr)
# up exp gene percentage of Accessory chromosome genes
strain = "fj11"
sm.geneID <- read.table(paste0("/Users/alexwang/0data/0mango/inf_associated_gene/SM/",strain,"_SM.geneID"), fill = T, header = T)[,1]

effector.geneID <- read.table(paste0("/Users/alexwang/0data/0mango/inf_associated_gene/secretome/",strain,"_effector.geneID"))[,1]
cazy.geneID <- read.table(paste0("/Users/alexwang/0data/0mango/inf_associated_gene/cazy/",strain,"_cazy.geneID"),fill = T)[,1]
up <- read.table(paste0("/Users/alexwang/0data/0mango/transcriptome/DEG/", strain, "_inf5d_Vs_hypha.up.geneID"))[,1]
down <- read.table(paste0("/Users/alexwang/0data/0mango/transcriptome/DEG/", strain, "_inf5d_Vs_hypha.down.geneID"))[,1]
core.genes <- read.table(paste0("/Users/alexwang/0data/0mango/accessory/core/", strain, "_core.geneID"))[,1]
accessory.gene <- read.table(paste0("/Users/alexwang/0data/0mango/accessory/", strain, "_accessory.geneID"))[,1]
## select candidate knock out gene
intersect(down, accessory.gene) %>% intersect(., sm.geneID)
intersect(up, accessory.gene) %>% intersect(., sm.geneID)

intersect(down, accessory.gene) %>% intersect(.,effector.geneID)
intersect(up, accessory.gene) %>% intersect(., effector.geneID)

intersect(down, accessory.gene) %>% intersect(.,cazy.geneID)
intersect(up, accessory.gene) %>% intersect(., cazy.geneID)

### plot percentage of up, down, others on accessory chromsomers
Up.accessory <- length(intersect(up, accessory.gene))
Down.accessory <- length(intersect(down, accessory.gene))
Other.acessory <- length(accessory.gene) - Up.accessory - Down.accessory
Up.core <- length(intersect(up, core.genes))
Down.core <- length(intersect(down, core.genes))
Other.core <- length(core.genes) - Up.core - Down.core
t = data.frame(accessory = c(Up.accessory, Down.accessory, Other.acessory),
           core = c(Up.core, Down.core, Other.core))
t.xtabs <- xtabs(t)
# chisq.test
chisq.test(t.xtabs)

#### plot group stack barplot
prob.accessory <- prop.table(t$accessory)
prob.core <- prop.table(t$core)
ma <- matrix(c(prob.accessory, prob.core),ncol = 2)
colnames(ma) <- c("accessory", "core")
## set color
cc <- c(rgb(112,47,160,maxColorValue = 255),rgb(8,186,255,maxColorValue = 255), rgb(255,128,2, maxColorValue = 255))
par(mar=c(1.5, 5, 5, 2))
barplot(ma, horiz = T, space = 0.5,col = cc,  las=1, xaxt="n")
legend(x = 0, y = 0.5, legend = c(">2-fold", "<-2-fold", "others"), 
       col = cc, xpd = T, pch=15, bty = "n", ncol = 3)
axis(side = 3, cex=0.5)
# text(x = 0.5, y=1.7, "*** X-squared p<2.2e-16")
title("Percentage of genes")

########### GO enricherment analysis
library(clusterProfiler)
godb <- read.delim("/Users/alexwang/0data/db/godb.txt",header = T,col.names = c("term","def", "ontology", "name"))
term2name <- godb[,c(1,4)]
term2gene <- read.table("/Users/alexwang/0data/0mango/EVM_protein/interproscan/gd10_gene2go.txt", header = T, col.names = c("gene", "term"))[,c(2,1)]
er <- enricher(accessory.gene, TERM2GENE = term2gene, TERM2NAME = term2name)
clusterProfiler::dotplot(er)


#######################################################******************** plot heatmap of strains with different pathogenicity

orthologues.df <- read.table("/Users/alexwang/0data/0mango/1_OrthoFinder/Orthologues_Dec28/Orthologues/Orthologues_yn56_protein/yn56_protein__v__gd10_protein.csv", fill = T, header = T, sep = "\t")
tpm1 <- read.table("/Users/alexwang/0data/0mango/transcriptome/tpm/yn56_tpm.txt")
tpm2 <- read.table("/Users/alexwang/0data/0mango/transcriptome/tpm/gd10_tpm.txt")

## select one-one mapping
orthologues.df <- orthologues.df[intersect(grep(",", orthologues.df$yn56_protein, invert = T), 
                                           grep(",", orthologues.df$gd10_protein, invert = T)), ]
getTpm <- function(x) {
  geneID1 = x[2]
  geneID2 = x[3]
  t1 <- as.numeric(tpm1[geneID1, c(7,8,9,10,11,12)])
  t2 <- as.numeric(tpm2[geneID2, c(7,8,9,10,11,12)])
  return(c(t1, t2))
  
}
orthologue.tpm.mat <- t(apply(orthologues.df, 1, getTpm))
colnames(orthologue.tpm.mat) <- c("yn56_inf3_1", "yn56_inf3_2", "yn56_inf3_3",
                                  "yn56_inf5_1", "yn56_inf5_2", "yn56_inf5_3",
                                  "gd10_inf3_1", "gd10_inf3_2", "gd10_inf3_3",
                                  "gd10_inf5_1", "gd10_inf5_2", "gd10_inf5_3")

## read up-regulated geneID
up.yn56 <- read.table("/Users/alexwang/0data/0mango/transcriptome/DEG/yn56_inf5d_Vs_coni.up.geneID")[,1]
up.gd10 <- read.table("/Users/alexwang/0data/0mango/transcriptome/DEG/gd10_inf5d_Vs_coni.up.geneID")[,1]

## select up-regulated gene pairs
idx <- which(orthologues.df$yn56_protein %in% up.yn56 & orthologues.df$gd10_protein %in% up.gd10)
orthologue.tpm.mat <- orthologue.tpm.mat[idx,]

gene_pair <- orthologues.df[idx,c(2,3)]
orthologue.tpm.df <- as.data.frame(orthologue.tpm.mat) %>% mutate(yn56_inf5_mean=(yn56_inf5_1 + yn56_inf5_2 + yn56_inf5_3)/3, 
                                                                  gd10_inf5_mean=(gd10_inf5_1 + gd10_inf5_2 + gd10_inf5_3)/3
) %>% .[,13:14]
res <- cbind(gene_pair, orthologue.tpm.df)
accessory.geneID <- read.table("/Users/alexwang/0data/0mango/accessory/accessory.geneID")[,1]
res <- res[res$gd10_protein %in% accessory.geneID,]
write.table(res, file = "/Users/alexwang/0data/0mango/transcriptome/tpm/inf5dVsCoinUp_tpm.txt", quote = F, sep = "\t", row.names = F)
# boxplot(orthologue.tpm.df$yn56_inf3_mean, orthologue.tpm.df$gd10_inf3_mean, ylim=c(0,400))
# t.test(orthologue.tpm.df$yn56_inf3_mean, orthologue.tpm.df$gd10_inf3_mean)

######################################################### boxplot of genepairs dn/ds, yn00
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
  text(x1+(x2-x1)/2, y-offset, labels = label, cex = 1.5)
}
# scale vector to range(mi, ma)
custom_scale <- function(v, mi, ma){
  (v - min(v)) * (ma-mi) / (max(v)-min(v)) + mi
}

########## CSEP vs random
CSEP.yn00 <- read.table("/Users/alexwang/0data/0mango/inf_associated_gene/secretome/paml/CSEP.yn00.txt", header = T)[,c(2,3,4)]
CSEP_control.yn00 <- read.table("/Users/alexwang/0data/0mango/inf_associated_gene/secretome/paml/CSEP_control.yn00.txt", header = T)[, c(2,3,4)]
CSEP.df <- cbind(CSEP.yn00, CSEP_control.yn00)
colnames(CSEP.df) <- c("omega", "dn", "ds", "omega.control", "dn.control", "ds.control")
CSEP.df <- custom_scale(CSEP.df, mi = 0, ma = 20)

par(xpd = T, 
    mfrow=c(1,3),
    mar=c(5.1, 5, 4.1, 0.5))
x1 = 1
x2 = 2
y = -0.08
tick_len = 0.02
offset = 0.05
boxplot(CSEP.df$omega, CSEP.df$omega.control,
        ylim = c(0, 1),
        col = c("#00BFFF", "#FF7F50"),
        axes = F,
        outline = F,
        frame = F,
        ylab = "Normalized ratio",
        cex.lab = 1.5
        )
add_axis(x1, x2, y, "dN/dS", tick_len, offset)
axis(2, lwd = 2, cex = 1)
t.test(CSEP.df$omega, CSEP.df$omega.control, alternative = "greater")
marksig(1,2, 1.1, "***", 0.02, 0.02)

boxplot(CSEP.df$dn, CSEP.df$dn.control,
        ylim = c(0, 1),
        col = c("#00BFFF", "#FF7F50"),
        axes = F,
        outline = F,
        frame = F
)
add_axis(x1, x2, y, "dN", tick_len, offset)
t.test(CSEP.df$dn, CSEP.df$dn.control, alternative = "greater")
marksig(1,2, 1.1, "****", 0.02, 0.02)

boxplot(CSEP.df$ds, CSEP.df$ds.control,
        ylim = c(0, 1),
        col = c("#00BFFF", "#FF7F50"),
        axes = F,
        outline = F,
        frame = F
)
add_axis(x1, x2, y, "dS", tick_len, offset)
t.test(CSEP.df$ds, CSEP.df$ds.control, alternative = "greater")
marksig(1,2, 1.1, "****", 0.02, 0.02)
legend(0, 1, legend = c("CSEP", "Others"), 
       fill = c("#00BFFF", "#FF7F50"), 
       bty = "n",
       border = NA,
       cex = 1.5)

###################### accessory vs core genes
accessory.yn00 <- read.table("/Users/alexwang/0data/0mango/inf_associated_gene/secretome/paml/accessory/accessory.yn00.txt", header = T)[,c(2,3,4)]
accessory_control.yn00 <- read.table("/Users/alexwang/0data/0mango/inf_associated_gene/secretome/paml/accessory/accessory_control.yn00.txt", header = T)[,c(2,3,4)]
# remove genepairs whose dN/dS = 99
accessory.yn00 <- accessory.yn00[accessory.yn00$omega != 99, ]
accessory_control.yn00 <- accessory_control.yn00[accessory_control.yn00$omega != 99, ]

accessory.df <- cbind(accessory.yn00, accessory_control.yn00[1:nrow(accessory.yn00), ])
colnames(accessory.df) <- c("omega", "dn", "ds", "omega.control", "dn.control", "ds.control")
accessory.df <- custom_scale(accessory.df, mi = 0, ma = 4)

par(xpd = T, 
    mfrow=c(1,3),
    mar=c(5.1, 5, 4.1, 0.5))
x1 = 1
x2 = 2
y = -0.08
tick_len = 0.02
offset = 0.05
boxplot(accessory.df$omega, accessory.df$omega.control,
        ylim = c(0, 1),
        col = c("#00BFFF", "#FF7F50"),
        axes = F,
        outline = F,
        frame = F,
        ylab = "Normalized ratio",
        cex.lab = 1.5
)
add_axis(x1, x2, y, "dN/dS", tick_len, offset)
axis(2, lwd = 2, cex = 1)
t.test(accessory.df$omega, accessory.df$omega.control, alternative = "greater")
marksig(1,2, 1.1, "*****", 0.02, 0.02)

boxplot(accessory.df$dn, accessory.df$dn.control,
        ylim = c(0, 1),
        col = c("#00BFFF", "#FF7F50"),
        axes = F,
        outline = F,
        frame = F
)
add_axis(x1, x2, y, "dN", tick_len, offset)
t.test(accessory.df$dn, accessory.df$dn.control, alternative = "greater")
marksig(1,2, 1.1, "*****", 0.02, 0.02)

boxplot(accessory.df$ds, accessory.df$ds.control,
        ylim = c(0, 1),
        col = c("#00BFFF", "#FF7F50"),
        axes = F,
        outline = F,
        frame = F
)
add_axis(x1, x2, y, "dS", tick_len, offset)
t.test(accessory.df$ds, accessory.df$ds.control)
marksig(1,2, 1.1, "*****", 0.02, 0.02)
legend(0, 1, legend = c("Accessory", "Core"), 
       fill = c("#00BFFF", "#FF7F50"), 
       bty = "n",
       border = NA,
       cex = 1.5)



