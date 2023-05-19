#!/bin/bash
## change RNA-seq bam chr id using two columns id_map file
set -eu

for i in *.bam;do
	strain=${i%%_*}
	PREFIX=${i%_sorted.bam}
	id_map=/disk/alpha/mango/id_map/${strain}.id
	echo ${PREFIX} ...
	# change BAM header
	python3 /disk/alpha/mango/id_map/changeBamHeader.py ${i} > ${PREFIX}.header
	# change BAM chr ID
	samtools view ${i} -h | awk -v name="$strain" 'BEGIN{FS=OFS="\t"}NR==FNR{a[name"_"$1]=name"_"$2}NR>FNR{if($1!~/@/){$3=a[$3];print}else{print}}' ${id_map} - \
		| samtools view -bS - | samtools reheader ${PREFIX}.header -  | samtools sort - ${PREFIX}.sort
done
