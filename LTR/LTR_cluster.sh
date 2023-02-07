#!/bin/bash
DIR=/Users/alexwang/0data/0mango/genome

## extract internal LTR sequence from LTR retriever result and cluster them using vmatcvh
for i in yn56 yn55 hn47 fj11 gd10 gz15 qz;do
	awk 'BEGIN{FS=OFS="\t";i=1}{split($1,a,":");split($7,b,"[\\.:]");gsub("\\?","+",$9);print a[1],b[2],b[4],"LTR"i"_"$10,".",$9;i++}' ${i}.fasta.pass.list > ${i}_LTR_internal.bed
	seqkit subseq --bed ${i}_LTR_internal.bed ${DIR}/${i}.fasta > ${i}_LTR_internal.fasta
done

# cluster usingn vmatch
#awk '/>/{gsub(">.* ",">",$0);print $0;next}{print}' LTR_internal.fasta > LTR_internal_rename.fasta
## clustering
#~/software/vmatch-2.3.1-Darwin_i386-64bit/mkvtree -db LTR_internal.fasta -dna -pl -allout -v
#~/software/vmatch-2.3.1-Darwin_i386-64bit/vmatch -dbcluster 95 7 Cluster-pred-chrAll -p -d -seedlength 50 -l 1101 -exdrop 9 -showdesc 100 LTR_internal.fasta

