#!/bin/bash
set -e
############## secretome, effectors
for i in qz hn47 fj11 yn55 yn56 gz15 gd10 gz14 gz23 yn51 hn42 yn31 yn32 yn33 gz23 hn32 hn23 Chi Cle Cgl Cde Ctr Csp Cac;do
    #seqkit split -f -p 8 ../${i}_protein.fasta
    #find ../${i}_protein.fasta.split -name "*.fasta" | parallel -j 4 deeploc -f {} -o {}_deepLoc
    #cat ../${i}_protein.fasta.split/${i}_protein*_deepLoc.txt > ../${i}_deepLoc.txt
    #awk '$2~/Extracellular/{print $1}' ../${i}_deepLoc.txt > ${i}_secretome.gene


    #signalp -fasta ../${i}_protein.fasta -mature -gff3 -prefix ${i}_signalp
	# predict transmembrane helices
    #~/software/tmhmm-2.0c/bin/tmhmm ${i}_signalp_mature.fasta > transme_helices.txt
    # predict subcellular localization
    #~/software/targetp-2.0/bin/targetp -fasta ${i}_signalp_mature.fasta -mature -gff3 -prefix ${i}_targetP

    ## remove transmembrane helices
    #grep "Number of predicted TMHs:  0" transme_helices.txt | awk '{print $2}' | sort | uniq | grep -f - ${i}_signalp.gff3 > secretome_rmHelices.txt

    ## remove mitochondrial localization
    #awk '$0~/Mitochondrion/{print $1}' ${i}_targetP.gff3 | grep -v -f - secretome_rmHelices.txt | awk -F "\t" '{print $1}' | sort | uniq > ${i}_secretome_signalP.geneID
    #rm secretome_rmHelices.txt
    #rm transme_helices.txt
    #rm ${i}_signalp*
    #rm ${i}_targetP*

    #cat ${i}_secretome_signalP.geneID ${i}_secretome.gene | sort | uniq > ${i}.TotalSecretome.geneID
	awk '/>/{print"";printf$0"|";next}{printf$0}' ~/0data/0mango/EVM_protein/${i}_protein.fasta | grep -f ${i}.TotalSecretome.geneID - | awk 'BEGIN{RS="\|"}{print}' > ${i}.TotalSecretome.fasta
	python ~/software/EffectorP_2.0/Scripts/EffectorP.py -i ${i}.TotalSecretome.fasta -o ${i}.effector.txt -E ${i}_effector.fasta
	grep "Effector probability" ${i}.effector.txt | awk -F "\|" '{print $1}' | sort | uniq > ${i}_effector.geneID

	## find CSEPs
	export LD_LIBRARY_PATH=/stor9000/apps/users/NWSUAF/2018055070/sf/glibc_2.14/lib
	#~/sf/diamond blastp --db /stor9000/apps/users/NWSUAF/2018055070/database/nr --query /stor9000/apps/users/NWSUAF/2018055070/data/mango/gene_anno_blast/secretome/${i}.TotalSecretome.fasta --out ${i}_TotalSecretome_outGenus.tbl --outfmt 6 qseqid sseqid mismatch gapopen length evalue staxids sscinames sskingdoms --salltitles --evalue 1e-5 --block-size 10 --threads 20 --taxon-exclude 5455
	awk '{print $1}' ${i}_TotalSecretome_outGenus.tbl | sort | uniq | grep -v -f - ${i}.TotalSecretome.geneID > ${i}.CSEP.geneID
done
#############################################
### get Second Metabolic backbone genes
for i in yn55 yn56 gz15 hn47 fj11 qz gd10;do
	#unzip ${i}.zip -d ${i}_antiSMASH
	#python3 /Users/alexwang/1.script/asm/CanS41/antismash_gbk2fasta.py ${i}_antiSMASH/*.region*.gbk > ${i}_antiSMASHcoreGenes.fasta
	#awk -F ";" '/>/{a[$2]++}END{for(i in a){print i,a[i]}}' ${i}_antiSMASHcoreGenes.fasta > ${i}_antiSMASH.summary
	## replace gene ID
	#awk 'NR==FNR{if($3~/transcript/){FS="\"";id[$4]=$2}}NR>FNR{if($0 ~ />/){split($0,a,";|>");printf ">"id[a[2]]";"a[3]"\n"}else{print}}' ~/0data/0mango/GTF/${i}_EVM.gtf ${i}_antiSMASHcoreGenes.fasta > ${i}_antiSMASHcoreGenes2.fasta
	#rm ${i}_antiSMASHcoreGenes.fasta
	#awk -F "[>;]" '/>/{print$2}' ${i}_antiSMASHcoreGenes2.fasta > ${i}_SM.geneID
	awk -F "[>;]" '/>/{print$2"\t"$3}' ${i}_antiSMASHcoreGenes2.fasta > ${i}_SM.geneClass
done
#############################################
### get CAZymes
for i in qz yn56 yn55 gd10 hn47 fj11 gz15;do
	awk -v strain="$i" 'BEGIN{total["GH"]=0;total["GT"]=0;total["CE"]=0;total["AA"]=0;total["CBM"]=0;total["PL"]=0}FNR>1 && $6>=2{printf $1"\t";array["GH"]=0;array["GT"]=0;array["CE"]=0;array["AA"]=0;array["CBM"]=0;array["PL"]=0;for(i=2;i<=4;i++){for(k in array){if($i~k){array[k]++}}};for(k in array){if(array[k]>=2){printf k",";total[k]++}};printf "\n"}END{print "GH:"total["GH"]"\nGT:"total["GT"]"\nCE:"total["CE"]"\nAA:"total["AA"]"\nCBM:"total["CBM"]"\nPL:"total["PL"] > sprintf("%s_cazy.summary",strain);close(sprintf("%s_cazy.summary",strain))}' ${i}_cazy.txt > ${i}_cazy.geneID
done
