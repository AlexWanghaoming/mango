#!/bin/bash

# core genes count
#awk 'FNR==1{print}{q=1;for(i=2;i<=NF;i++){q=q*$i};if(q!=0){print}}' Orthogroups.GeneCount.csv | awk '{for(i=2;i<=NF;i++){if(FNR==1){id[i]=$i}else{a[i]=a[i]+$i}}}END{for(i=2;i<=NF;i++){print a[i]"\t"id[i-1]}}'

# strain-specific genes count
#awk -F "\t" 'NR==4 || NR==10|| NR==1{for(i=2;i<=NF;i++){if(NR==1){id[i]=$i};a[i]=a[i]+$i}}END{for(i=2;i<=NF;i++){print a[i],id[i]}}' Statistics_PerSpecies.csv

# C.gloeosporioides Complex specific orthogroups
awk 'NR==1{print}NR>1&&$2==0&&$3==0&&$7==0&&$8==0&&$9==0&&$10==0&&$15==0&&$18==0&&$22==0&&$24==0&&$25==0{print}' Orthogroups.GeneCount.csv > CGSC.orthogroups.GeneCount.csv
awk 'NR>1{print $1}' CGSC.orthogroups.GeneCount.csv |xargs -n 1 -I {} grep {} Orthogroups.txt > CGSC.orthogroups.csv

# C.gloeosporioides Complex orthogroups conserved in all strains
awk 'NR==1{print}NR>1&&$2==0&&$3==0&&$7==0&&$8==0&&$9==0&&$10==0&&$15==0&&$18==0&&$22==0&&$24==0&&$25==0 && $4!=0&&$5!=0&&$6!=0&&$11!=0&&$12!=0&&$13!=0&&$14!=0&&$16!=0&&$17!=0&&19!=0&&$20!=0&&$21!=0&&$23!=0{print}' Orthogroups.GeneCount.csv > CGSC_conserved.orthogroups.GeneCount.csv
awk 'NR>1{print $1}' CGSC_conserved.orthogroups.GeneCount.csv |xargs -n 1 -I {} grep {} Orthogroups.txt > CGSC_conserved.orthogroups.csv

## Effectors OGs
grep -f /Users/alexwang/0data/0mango/inf_associated_gene/secretome/all_effector.geneID Orthogroups.txt | sort | uniq > effector_Orthogroups.txt




