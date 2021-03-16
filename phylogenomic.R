## .fasta -> .phy
# python3 convertFmt.py
## trim MSA
# trimal -in ~/0data/0mango/1_OrthoFinder/Orthologues_Dec28/Alignments/SpeciesTreeAlignment.phy 
#    -out ~/0data/0mango/1_OrthoFinder/Orthologues_Dec28/Alignments/SpeciesTreeAlignment.trimAl.phy -automated1
## build ML tree
#phyml -i /disk/alpha/mango/SpeciesTreeAlignment.trimAl.phy -d aa -m LG -f m -v e -a e -o tlr

## root ML tree
#tree <- ggtree::read.tree("/Users/alexwang/0data/0mango/cafe_res/SpeciesTreeAlignment.trimAl.phy_phyml_tree.r8s.txt")
#tree <- phytools::midpoint.root(tree)
#treeio::write.tree(tree, file = "/Users/alexwang/0data/0mango/cafe_res/SpeciesTreeAlignment.trimAl.phy_phyml_tree.r8s.rooted.txt")

### infer gene family evolution
#python2 cafetutorial_prep_r8s.py -i SpeciesTreeAlignment.trimAl.phy_phyml_tree.r8s.rooted.txt -o r8s_ctl_file.txt -s 3360996 -p 'Chi_protei,Cac_protei' -c 54

### build ultrametric tree using r8s in docker container
#docker run -it --rm -v $(pwd):/data shkao/r8s:1.81 r8s -b -f /data/r8s_ctl_file.txt > r8s_tmp.txt

### run CAFE
# awk 'BEGIN{OFS="\t"}NR==1{$NF="";gsub("^\t","",$0);gsub("_protein","",$0);print "Desc","Family ID",$0;next}{$NF="";print "(null)",$0}' Orthogroups.GeneCount.csv > cafe.data
# cafe cafetutorial_run.sh
# visualization
# python3 CAFE_fig.py resultfile.cafe -pb 0.05 -pf 0.05 --dump test/ -g pdf --count_all_expansions
# python2 cafetutorial_report_analysis.py -i resultfile.cafe -o summary_cafe
library(ggtree)
library(ggplot2)
library(ggsci)
library(ggpubr)
tree <- read.tree(file = "/Users/alexwang/0data/0mango/SpeciesTreeAlignment.trimAl.phy_phyml_tree.txt")
tree <- phytools::midpoint.root(tree)
tree$tip.label <- stringr::str_replace_all(string = tree$tip.label, 
                pattern = tree$tip.label,
                replacement = c("C.sp.JS-367","C.gigasporum HN42*","C.truncatum","C.gloeosporioides","C.gloeosporioides GZ15*",
                                "C.asianum FJ11*","C.asianum YN55*", "C.endophytica YN32*", "C.tropicale HN32*",
                                "C.siamense GD10*", "C.siamense YN56*","C.fructicola HN47*", "C.karstii GZ14*",
                                "C.fructicola 1104","C.fructicola Nara_gc","C.fructicola QZ3*","C.horii PLG3*","C.cordylinicola HN23*",
                                "C.cliviae YN31*","C.musae GZ23*", "C.higginsianum","C.destructivum","C.lentis",
                                "C.liaoningense YN33*","C. scovillei YN51*","C.acutatum"))
tree1 <- 
  ggtree(tree) + geom_hilight(node = 32, "steelblue") + 
  geom_tiplab(align = T,linetype = "dotted", size=3) + xlim(0,0.23)+ 
  geom_point(color='black', size=0.8) 
# + geom_treescale() 
# ## check node lables
# tree1 + geom_text(aes(label=node),hjust=-.3,color="red") 
# ## get subset tree
# keep = c("Cfr1_protein", "Cfr2_protein", "gz14_protein", "qz_protein", "hn47_protein")
# subTree1 <- treeio::drop.tip(tree, tree$tip.label[!tree$tip.label %in% keep])
# ggtree(subTree1) + geom_tiplab(align = T, offset = 0.003,linetype = "dotted")+xlim(0, 0.03)+geom_point(color='black', size=1)
# ##############################################################################
library(readxl)
# change color here
cc <- ggsci::pal_lancet('lanonc',alpha = 0.6)(9)
info <- read_xlsx("~/Desktop/1genome_assembly_infomation.xlsx", 
                  col_names = T, skip = 1,col_types = c("text",rep("numeric",35)))[c(3,25,26,15,6,5,4,22,16,2,1,7,10,24,23,14,18,21,19,13,11,9,12,8,17,20),
                                           c(1,4,10,11,12,13,15,16,17,18,19,26,35,36)]
colnames(info) <- c("Pathogen","Genome", "Gene", "Core","Share","Unique","Retroelement", "DNA transposon","Unclassified", "Simple repeat",
                    "CAZymes","SM_cluster","Secretome","Effector")
positions <- info$Pathogen
# genome size
genome.p <- ggplot(data=info) + geom_bar(aes(x=Pathogen,y=Genome),fill=cc[1],stat="identity")+ 
              geom_text(aes(x=Pathogen,y=20,label=Genome),size=3) + scale_x_discrete(limits=rev(positions)) +
              scale_y_continuous(position = "right", expand =c(0,0)) + coord_flip() +
              theme( axis.title.y = element_blank(), axis.title.x = element_text(size=10),
                     axis.text= element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     panel.background = element_blank())  

# gene content
gene.p <- ggplot(data=info) + geom_bar(aes(x=Pathogen,y=Gene),stat="identity",fill=cc[7])+
            geom_text(aes(x=Pathogen,y=8000,label=Gene),size=3) + scale_x_discrete(limits=rev(positions)) + 
            coord_flip() + scale_y_continuous(position = "right", expand =c(0,0)) + 
            theme( axis.title.y = element_blank(), axis.title.x = element_text(size=10),
                   axis.text= element_blank(),axis.line = element_blank(), 
                   axis.ticks = element_blank(), panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(), panel.background = element_blank())  
# proteins 
protein.df <- reshape2::melt(info[,c(1,4,5,6)])
protein.df$variable <- factor(protein.df$variable, levels = c("Unique","Share","Core"))
protein.p <- ggplot(data=protein.df) + geom_bar(aes(x=Pathogen, y=value, fill=variable, group=variable), stat="identity", position = "stack") + 
  scale_x_discrete(limits=rev(positions)) + coord_flip() + scale_y_continuous(position = "right", expand = c(0,0)) + 
  ylab("No. of proteins") + guides(color = guide_legend(override.aes = list(size = 0.03))) + 
  theme(legend.title = element_blank(), legend.position = c(0.4, 0.6),
        legend.key.size = unit(0.3,'cm'),
        axis.title.y = element_blank(), axis.title.x = element_text(size=10),
        axis.text.y = element_blank(),axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),axis.line.x = element_line(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        panel.background = element_blank())
protein.p <- protein.p + scale_fill_lancet(alpha = 0.6)

# Repeat content
repeat.df <- reshape2::melt(info[,c(1,7,8,9,10)])
repeat.df$variable <- factor(repeat.df$variable, levels = c("Simple repeat","DNA transposon","Unclassified","Retroelement"))
repeat.p <- ggplot(data=repeat.df) +
  geom_bar(aes(x=Pathogen, y=value, fill=variable, group=variable), stat="identity", position = "stack") + 
  scale_x_discrete(limits=rev(positions)) + coord_flip() + scale_y_continuous(position = "right", expand = c(0,0)) + 
  ylab("Repeat content (%)") + guides(color = guide_legend(override.aes = list(size = 0.1))) + 
  theme(legend.title = element_blank(), legend.position = c(0.7, 0.6),
        legend.key.size = unit(0.3,'cm'),
        axis.title.y = element_blank(), axis.title.x = element_text(size=10),
        axis.text.y = element_blank(),axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),axis.line.x = element_line(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        panel.background = element_blank())
repeat.p <- repeat.p + scale_fill_lancet(alpha = 0.6)

# CAZymes
CAZY.p <- ggplot(data=info) + geom_bar(aes(x=Pathogen,y=CAZymes),stat="identity",fill=cc[5])+
  geom_text(aes(x=Pathogen,y=300,label=CAZymes),size=3) + scale_x_discrete(limits=rev(positions)) + 
  coord_flip() + scale_y_continuous(position = "right", expand =c(0,0)) + 
  theme( axis.title.y = element_blank(), axis.title.x = element_text(size=10),
         axis.text= element_blank(),axis.line = element_blank(), 
         axis.ticks = element_blank(), panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(), panel.background = element_blank()) 

# secretome
secretome.p <- ggplot(data=info) + geom_bar(aes(x=Pathogen,y=Secretome),stat="identity",fill=cc[4])+
  geom_text(aes(x=Pathogen,y=1000,label=Secretome),size=3) + scale_x_discrete(limits=rev(positions)) + 
  coord_flip() + scale_y_continuous(position = "right", expand =c(0,0)) + 
  theme( axis.title.y = element_blank(), axis.title.x = element_text(size=10),
         axis.text= element_blank(),axis.line = element_blank(), 
         axis.ticks = element_blank(), panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(), panel.background = element_blank())  

# Secondary Metabolic cluster
SM.p <- ggplot(data=info) + geom_bar(aes(x=Pathogen,y=SM_cluster),stat="identity",fill=cc[6])+
  geom_text(aes(x=Pathogen,y=20,label=SM_cluster), size=3) + scale_x_discrete(limits=rev(positions)) + 
  coord_flip() + scale_y_continuous(position = "right", expand =c(0,0)) + 
  theme( axis.title.y = element_blank(), axis.title.x = element_text(size=10),
         axis.text= element_blank(),axis.line = element_blank(), 
         axis.ticks = element_blank(), panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(), panel.background = element_blank())  

# Effector
effector.p <- ggplot(data=info) + geom_bar(aes(x=Pathogen,y=Effector),stat="identity",fill=cc[9])+
  geom_text(aes(x=Pathogen,y=200,label=Effector),size=3) + scale_x_discrete(limits=rev(positions)) + 
  coord_flip() + scale_y_continuous(position = "right", expand =c(0,0)) + 
  theme( axis.title.y = element_blank(), axis.title.x = element_text(size=10),
         axis.text= element_blank(),axis.line = element_blank(), 
         axis.ticks = element_blank(), panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(), panel.background = element_blank())

cowplot::plot_grid(tree1, genome.p, repeat.p, gene.p, protein.p, CAZY.p, secretome.p, SM.p, effector.p, 
                   nrow = 1, align = T, rel_widths = c(5,1,2,1,2,1,1,1,1))


## Cazymes heatmap
require(ggplotify)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
va <- as.matrix(info[,12:17])
as.vector(va)
# heat.p <- pheatmap(mat = info[,12:17], scale = "column", cluster_rows = F, cluster_cols = F,
#           display_numbers = va, show_rownames = F, show_colnames = F, legend = F)

# heat.p <- Heatmap(mat_scaled, cluster_rows = F, cluster_columns = F, column_names_side = "top",
#                   column_title = "CAZY", column_title_side = "top",show_heatmap_legend = F,column_names_rot = T,
#                   col = colorRamp2(c(-2, 0, 2), rev(brewer.pal(n=3, name="RdBu"))),cell_fun = function(j, i, x, y, width, height, fill) {
#                     grid.text(sprintf("%d", as.matrix(info[,12:17])[i, j]), x, y, gp = gpar(fontsize = 10))
#                   }
#                   )
# heat.p <- as.ggplot(heat.p)

mat_scaled = data.frame(apply(info[,12:17], 2, scale))
mat_scaled$Pathogen <- info$Pathogen
mat_scaled <- reshape2::melt(mat_scaled)
heat.p <- ggplot(mat_scaled, aes(x=Pathogen, y=variable, fill=value)) + guides(fill=F) + geom_tile() + 
            geom_text(label=as.vector(va), size=2) + coord_flip() + scale_fill_gradient2(low = "navy",mid = "white", high = "firebrick") + 
            scale_x_discrete(limits=rev(positions)) + ylab("CAZY") +scale_y_discrete(position = "right", expand =c(0,0)) + 
            theme( axis.title.y = element_blank(), axis.title.x = element_text(size=10),
                   axis.text.y = element_blank(), axis.text.x = element_text(),axis.line = element_blank(), axis.ticks = element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                   panel.background = element_blank())  


#######################################################################################################################################################
############################# accessory chromosomes syntenic heatmap

library(GenomicFeatures)
library(ggtree)
library("ggplot2")
GTF_file = "/Users/alexwang/0data/0mango/synteny/accessory/accessory.all.gtf"
LAST.simple = "/Users/alexwang/0data/0mango/synteny/accessory/accessory.accessory.anchors.new"

txdb <- makeTxDbFromGFF(file = GTF_file)
all.gene <- genes(txdb)
last.simple.df <- read.table(LAST.simple)
miniChr.id <- unique(seqnames(all.gene))

# store syntenic scores
miniChr.mat <- matrix(0, length(miniChr.id), length(miniChr.id))
colnames(miniChr.mat) <- miniChr.id
rownames(miniChr.mat) <- miniChr.id

# count syntenic gene numbers
for(k in 1:nrow(last.simple.df)){
  qry = as.character(seqnames(all.gene[all.gene$gene_id ==  last.simple.df[k, "V1"]]))
  target = as.character(seqnames(all.gene[all.gene$gene_id == last.simple.df[k, "V2"]]))
  
  if (qry != target) {
    miniChr.mat[qry, target] = miniChr.mat[qry, target] + 1
    miniChr.mat[target, qry] = miniChr.mat[target, qry] + 1
  }
}

## calculate syntenic scores
for (i in miniChr.id) {
  for (j in miniChr.id) {
    if (i != j) {
      c1 = NROW(all.gene[seqnames(all.gene) == i])
      c2 = NROW(all.gene[seqnames(all.gene) == j])
      
      num = miniChr.mat[i, j]
      
      miniChr.mat[i, j] = ((num/c1) * (num/c2)) / (num/c1 + num/c2)
    }else{
      # The score on the diagonal is 0.5
      miniChr.mat[i, j] = 0.5
    }
  }
}

miniChr.mat[is.nan(miniChr.mat)] = 0
id.order = c(paste0("qz_",11:13),
       paste0("Cfr1_", 11:12), 
       paste0("Cfr2_", 12:13), 
       paste0("gz14_", 12:14),
       paste0("hn47_", 11:12), 
       paste0("yn56_", 11:17),
       paste0("gd10_", 11:14),
       paste0("hn32_", 11:13), 
       paste0("yn32_", c(12:17, 19,20)),
       paste0("yn55_", 11:13),
       paste0("fj11_", 11:16), 
       paste0("gz15_", 11:16),
       paste0("Cgl_", c(2,5,6,7,8,9,13)),
       paste0("hn23_", 12:19), 
       paste0("plg3_", 11:17),
       paste0("gz23_", 11:14), 
       paste0("Cde_", 56),
       paste0("Chi_", 24:25), 
       paste0("Cle_", c(6, 8, 13)), 
       paste0("Cac_", 21),
       paste0("yn51_", 15))

miniChr.mat <- miniChr.mat[id.order, id.order]

library(ComplexHeatmap)
library(circlize)
# set colors
cc <- c("white", rgb(255,128,2, maxColorValue = 255))
# set annotations
my_sample <- factor(unlist(lapply(strsplit(id.order, "_"), function(x) x[1])), 
                    levels = unique(unlist(lapply(strsplit(id.order, "_"), function(x) x[1]))))
my_sample_anno <- rowAnnotation(Strain = my_sample,
                                show_annotation_name = FALSE)

heat.p <- Heatmap(miniChr.mat, cluster_rows = F, cluster_columns = F, column_title_side = "top", show_heatmap_legend = T,
        column_names_rot = 45, col=colorRamp2(c(0,0.5), colors = cc), column_names_centered = T,
        ## set for each cell
        rect_gp = gpar(col = "grey", lwd = 0.5),
        border = F, 
        # font of colnames, rownames
        row_names_gp = gpar(fontsize=6), column_names_gp = gpar(fontsize = 5),
        # split rows and annotations
        left_annotation = my_sample_anno, row_split = my_sample, row_title = NULL, gap = unit(1, "mm"),
        # set legends
        heatmap_legend_param = list(title="Syntenic score", 
                                    at=c(0,0.1,0.2,0.3,0.4,0.5), 
                                    labels=c(0,0.1,0.2,0.3,0.4,0.5),
                                    legend_height = unit(3, "cm")))
# generate ggplot2 object
heat.gg <- ggplotify::as.ggplot(heat.p)
# add phylogenetic tree
tree <- read.tree(file = "/Users/alexwang/0data/0mango/SpeciesTreeAlignment.trimAl.phy_phyml_tree.txt")
cladeA = c("yn33_prote", "Csp_protei", "hn42_prote", "Ctr_protei", "yn31_prote")
tree <- ape::drop.tip(tree, cladeA)
tree <- phytools::midpoint.root(tree)
tree.gg <- ggtree(tree, branch.length = "none", size=0.8)
            # geom_text(aes(label=node),color="red", size=3) + 
            # geom_tiplab(align = T,linetype = "dotted", size=3)

cowplot::plot_grid(plotlist = list(tree.gg, heat.gg), align = "h", rel_widths = c(1,3), nrow = 1)
