minimap2 -ax map-pb -t 30 ../genome/fj11.fasta ../FJ11-1A.filtered_subreads.fa.gz | samtools view -@ 10 -bS - | samtools sort -@ 10 -o fj11_pbSorted.bam -
minimap2 -ax map-pb -t 16 ../genome/qz.fasta /stor9000/apps/users/NWSUAF/2018055070/data/mango/qz_subreads_merged.fasta | samtools view -@ 10 -bS - | samtools sort -@ 10 -o qz_pbSorted.bam -
minimap2 -ax map-pb -t 30 ../genome/gz15.fasta ../GZ15-1A.filtered_subreads.fa.gz | samtools view -@ 10 -bS - | samtools sort -@ 10 -o gz15_pbSorted.bam -
minimap2 -ax map-pb -t 30 ../genome/gd10.fasta ../GD10-1A.filtered_subreads.fa.gz | samtools view -@ 10 -bS - | samtools sort -@ 10 -o gd10_pbSorted.bam -
minimap2 -ax map-pb -t 20 ../genome/yn55.fasta ../YN55-1A.filtered_subreads.fa.gz | samtools view -@ 10 -bS - | samtools sort -@ 10 -o yn55_pbSorted.bam -
minimap2 -ax map-pb -t 30 ../genome/yn56.fasta ../YN56-1-1A.filtered_subreads.fa.gz | samtools view -@ 10 -bS - | samtools sort -@ 10 -o yn56_pbSorted.bam -
minimap2 -ax map-pb -t 30 ../genome/hn47.fasta ../HN47-2A.filtered_subreads.fa.gz | samtools view -@ 10 -bS - | samtools sort -@ 10 -o hn47_pbSorted.bam -
## calculate coverage in bins
for i in yn55 qz hn47 gz15 gd10 fj11 yn56;do
	mosdepth -b 5000 -n -t 4 ${i}_pbSorted.bam
	gzip -d ${i}_pb.regions.bed.gz
done
