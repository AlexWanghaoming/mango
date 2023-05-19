#!/bin/bash
# awk 'NR==FNR{a[$1]=$1}NR>FNR{RS=">";FS="\n";if ($1 in a) {printf ">%s",$0}}' ../HGT/putative_HGT.geneID all.cds > putative_HGT.cds.fasta
# source activate emboss
# Use all chromosomes of Colletotrichum build codon usage table
cusp -sequence all.cds -outfile Colletotrichum.cusp
# HGT genes
cai -seqall putative_HGT.cds.fasta -cfile Colletotrichum.cusp -outfile HGT.cai
# core chromosomes
cai -seqall all.cds -cfile Colletotrichum.cusp -outfile Colletotrichum.cai
# gz23
awk 'NR==FNR{a[$1]=$1}NR>FNR{RS=">";FS="\n";if ($1 in a) {printf ">%s",$0}}' ../HGT/gz23_HGT.geneID all.cds > gz23_HGT.cds.fasta
cai -seqall gz23_HGT.cds.fasta -cfile Colletotrichum.cusp -outfile gz23_HGT.cai
