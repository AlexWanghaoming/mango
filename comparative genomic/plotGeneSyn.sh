#!/bin/bash
GTF=/Users/alexwang/0data/0mango/GTF
GENOME=/Users/alexwang/0data/0mango/genome
## extract cds,pep
#for i in gd10 yn56 yn55 fj11 qz hn47 gz15 Csp Cac Cde Cgl Cle Chi Ctr gz14 gz23 hn23 hn32 hn42 plg3 yn31 yn32 yn33 yn51 Cfr1 Cfr2;do
#	python -m jcvi.formats.gff bed --type=gene --key=gene_id ${GTF}/${i}_EVM.gtf -o ${i}.bed
#	gffread ${GTF}/${i}_EVM.gtf -g ${GENOME}/${i}.fasta -x ${i}.cds.fasta
#	awk '/>/{split($0,a,"=");print ">"a[2];next}{print}' ${i}.cds.fasta > ${i}.cds
#	rm ${i}.cds.fasta
	## pep
#	cp /Users/alexwang/0data/0mango/EVM_protein/${i}_protein.fasta ${i}.pep
#done

## find synteny relationship,this step must have cds fasta file in current directory cds fasta file must have the same headers as bed file above
##### outCg strains
#strains=(Ctr gz23 yn31 Csp hn42 Cac yn51 yn33 Cde Chi Cle)
##### inCg strains
#strains=(qz Cfr1 Cfr2 gz14 hn47 yn56 gd10 hn32 yn32 yn55 fj11 gz15 Cgl hn23 plg3)
##### represent strains
strains=(qz yn56 hn32 yn32 yn55 gz15 hn23 plg3)
ll=${#strains[@]} 
for i in `seq 0 $[${ll}-2]`;do
	## find ortholog
	qry=${strains[${i}]}
	target=${strains[$[${i}+1]]}
	#### nucl
	#python -m jcvi.compara.catalog ortholog ${qry} ${target} --cscore=.99
	#### prot
#	python -m jcvi.compara.catalog ortholog --dbtype=prot ${qry} ${target} --cscore=.99
	# bundle blocks
#	python -m jcvi.compara.synteny screen --minspan=10 --simple ${qry}.${target}.anchors ${qry}.${target}.anchors.new
#	## highlight accessory ribbon
	if [ ! -f "${qry}.${target}.anchors.simple.bak" ];then
		mv ${qry}.${target}.anchors.simple ${qry}.${target}.anchors.simple.bak
		awk 'NR==FNR{gene[$1]=$1}NR>FNR{for(i=1;i<=NF-2;i++){flag=0;if($i in gene){$0="g*"$0;print $0;flag=1;break}}if(flag!=1){print $0}}' /Users/alexwang/0data/0mango/accessory/accessory.geneID ${qry}.${target}.anchors.simple.bak > ${qry}.${target}.anchors.simple
	fi
done

## generate seqids
##rm seqids
#for i in `seq 0 $[${ll}-1]`;do
#	awk -F ">" '/>/{{printf$2","}}END{print ""}' ${GENOME}/${strains[${i}]}.fasta >> seqids
#done
#sed 's/,$//g' seqids > seqids.outCg
#rm seqids

## before plot, you should manual setting .seqid file to many lines
#python -m jcvi.graphics.karyotype --figsize=20x15 seqids2 layout.cg
