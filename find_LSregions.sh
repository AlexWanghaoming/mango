#!/bin/bash
set -eu
## This script is used to define Lineage Specific genomic regions and conserved regions of 7 strains

bundle_dir=/Users/alexwang/0data/0mango/genome/mummer
SHELL_DIR=$(cd "$(dirname "$0")";pwd)
OUT_DIR=/Users/alexwang/0data/0mango/genome/LS_region

for i in `ls ${bundle_dir}/*.bundle_telo3.txt`
do
    fn=$(basename ${i})
    sed -n '2,$p' ${i} | awk 'BEGIN{FS=OFS="\t"}$3-$2>20000{print $1,$2,$3,"-"}' > ${SHELL_DIR}/${fn%.txt}.bed

done

for strain in gd10 gz15 hn47 qz yn55 yn56 fj11
do
    cat ${strain}*.bundle_telo3.bed > ${strain}_merge.bed
    bedtools sort -i ${strain}_merge.bed > ${strain}_merge_sort.bed
    bedtools genomecov -i ${strain}_merge_sort.bed -g ${strain}.genome -bga > ${strain}_disjoint.bed

    #select conserved and LS region
    awk '$4<=1' ${strain}_disjoint.bed > ${OUT_DIR}/${strain}.specific.txt
    awk '$4>=5' ${strain}_disjoint.bed > ${OUT_DIR}/${strain}.conserved.txt
done

