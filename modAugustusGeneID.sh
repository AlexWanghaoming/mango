#!/bin/bash

# rename augustus gene id and extract protein data
grep -v "^#" braker/augustus.hints.gtf | awk 'BEGIN{FS=OFS="\t"}$9!~/m1-g/{print $0}'|sed 's/gd10_//g' | awk 'BEGIN{FS=OFS="\t"}$3!~/gene/{$5=$5-1;print;next}{print}'\
	| sort -k1,1n -k4,4n -k5,5nr -k3,3r | awk 'BEGIN{FS=OFS="\t"}{$1="gd10_"$1;if($3!~/gene/){$5=$5+1};print}'\
	| awk 'BEGIN{FS=OFS="\t";i=0}$3~/gene/{i++;zz=sprintf("%05d",i);id[$9]="CSIGD"zz;gsub($9,id[$9], $9);print;next}\
	{split($9,b,"\"|\\.");gsub(b[2],id[b[2]],$9);print}' > gd10.gtf
gffread gd10.gtf -g ../../../c.g/gd10.fasta -y gd10.protein.aa

