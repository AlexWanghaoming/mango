#!/bin/bash

nucmer=/Users/alexwang/software/mummer4/bin/nucmer
minMatchLength=5000
GENO_DIR=/Users/alexwang/0data/0mango/genome
telo_DIR=/Users/alexwang/0data/0mango/genome/telo
strains=(qz yn56 yn55 gz15 gd10 fj11 hn47)
for strain in ${strains[@]}
do
    remain=${strains[@]/${strain}}
    for qry in ${remain}
    do
        echo "${strain}_${qry}"
        #echo "step1: run nucmer"
        #${nucmer} --maxmatch --maxgap=500 --mincluster=100 -t 6 --prefix=${strain}_${qry} ${GENO_DIR}/${strain}.fasta ${GENO_DIR}/${qry}.fasta 
        #echo "show-coords -d -l -L ${minMatchLength} -r -c -T ${strain}_${qry}.delta > ${strain}_${qry}.tab2.txt"
        #show-coords -d -l -L ${minMatchLength} -r -c -T ${strain}_${qry}.delta > ${strain}_${qry}.tab2.txt
        
        echo "step2: bundle synteny blocks and assign telomere. "
        echo "python3 mummer_bundle_telo.py ${strain}_${qry}.tab2.txt ${telo_DIR}/${qry}.telo.txt"
        python3 mummer_bundle_telo.py ${strain}_${qry}.tab2.txt ${telo_DIR}/${qry}.telo.txt
    done
done
