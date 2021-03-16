#!/bin/bash
GTF=/Users/alexwang/0data/0mango/GTF
GENOME=/Users/alexwang/0data/0mango/genome
for i in qz Cfr1 Cfr2 hn47 gz14 yn56 gd10 hn32 yn32 yn55 fj11 gz15 Cgl hn23 plg3 Ctr gz23 yn31 Csp hn42 Cac yn51 yn33 Cde Chi Cle;do
	## accessory gene id
	awk -v strain="$i" -F "\t" 'NR==FNR && $1 ~ strain{split($2,id,",")}NR>FNR && $3~/gene/{for(i in id){if($1 == strain"_"id[i]){split($9,gene,"\"");print gene[2]}}}' /Users/alexwang/0data/0mango/accessory/accessory_chrom.ID ${GTF}/${i}_EVM.gtf > ${i}_accessory.geneID	

	## accessory genome fasta
	awk -v strain="$i" -F "\t" 'NR==FNR && $1 ~ strain{split($2,id,",");for(i in id){vv[id[i]]=id[i]}}NR>FNR{RS=">";FS="\n";split($0,a,"_");if(a[2] in vv){printf ">%s",$0}}' accessory_chrom.ID ${GENOME}/${i}.fasta > ${i}_accessory.genome.fasta
	
	## accessory GTF
	awk -v strain="$i" -F "\t" 'NR==FNR && $1 ~ strain{split($2,id,",");for(i in id){vv[id[i]]=id[i]}}NR>FNR{split($1,a,"_");if(a[2] in vv){print}}' accessory_chrom.ID ${GTF}/${i}_EVM.gtf > ${i}_accessory.gtf

	## count gene ortholog on chromosomes
#	python3 OrthoCountOnChrom.py ${i}_accessory.geneID ${i}_accessoryGeneCopy.count.txt

done
# select accessory specific orthogroups
awk 'NR==FNR{id[$1]=$1}NR>FNR{total=0;for(i=2;i<=NF;i++){if($i in id){total=total+1}};if(total==NF-1){print $0}}' accessory.geneID ../1_OrthoFinder/Orthogroups.txt > Orthogroups_accessory.txt 
# select specific genes on accessory chromosomes
awk '{for(i=2;i<=NF;i++){print $i}}' Orthogroups_accessory.txt > accessory.specific.geneID
# count for each strain
awk '{s=substr($1,1,5);a[s]++}END{for(i in a){print i,a[i]}}' accessory.specific.geneID
# merge accessory genome
cat qz_accessory.genome.fasta Cfr1_accessory.genome.fasta Cfr2_accessory.genome.fasta gz14_accessory.genome.fasta hn47_accessory.genome.fasta yn56_accessory.genome.fasta gd10_accessory.genome.fasta hn32_accessory.genome.fasta yn32_accessory.genome.fasta yn55_accessory.genome.fasta fj11_accessory.genome.fasta gz15_accessory.genome.fasta Cgl_accessory.genome.fasta hn23_accessory.genome.fasta plg3_accessory.genome.fasta gz23_accessory.genome.fasta Cac_accessory.genome.fasta yn51_accessory.genome.fasta Cde_accessory.genome.fasta Chi_accessory.genome.fasta Cle_accessory.genome.fasta > accessory.genome.all.fasta
cat *_accessory.gtf > accessory.all.gtf
####################################################################################### 
# merge all proteins
cat /Users/alexwang/0data/0mango/EVM_protein/*_protein.fasta > core/all.protein.fasta
# extract accessory proteins
awk 'NR==FNR{a[$1]=$1}NR>FNR{RS=">";FS="\n";if ($1 in a) {printf ">%s",$0}}' accessory.geneID core/all.protein.fasta > accessory.all.protein.fasta
# extract core proteins
awk -F "\t" '$3~/gene/{split($9, id, "\"");print id[2]}' ~/0data/0mango/GTF/*_EVM.gtf | grep -v -f accessory.geneID > core/core.geneID
awk 'NR==FNR{a[$1]=$1}NR>FNR{RS=">";FS="\n";if ($1 in a) {printf ">%s",$0}}' core/core.geneID core/all.protein.fasta > core/core.all.protein.fasta
# blastp accessory on core proteins
makeblastdb -dbtype prot -in core/core.all.protein.fasta -parse_seqids
blastp -db core/core.all.protein.fasta -query accessory.all.protein.fasta -outfmt 6 -evalue 1e-5 -out res.blastp.txt
# get accessory specific proteins
awk '{print $1}' res.blastp.txt | sort | uniq | sort accessory.geneID - - | uniq -u > accessory_sp.blastp.geneID
awk 'NR==FNR{a[$1]=$1}NR>FNR{RS=">";FS="\n";if ($1 in a) {printf ">%s",$0}}' accessory_sp.blastp.geneID accessory.all.protein.fasta > accessory_sp.blastp.protein.fasta
# get potential HGT genes
awk '$12 !~ /N\/A/ && $12 !~ /Colletotrichum/{print}' HGT/accessory_sp.blastp.protein.tbl > HGT/accessory_sp.blastp.protein.outColletotrichum.tbl
# select best hits
awk '{if($1 in gene){next}else{gene[$1]=$1;print $0}}' accessory_sp.blastp.protein.outColletotrichum.tbl > putative_HGT.blastp.tbl
#####################################################################################################################################
# select gz23 potential regions           CMUGZ_13712 - CMUGZ_13718
seqkit subseq --chr gz23_12 -r 155317:193760 gz23_accessory.genome.fasta > > CMUGZ_13712-CMUGZ_13724.fasta
