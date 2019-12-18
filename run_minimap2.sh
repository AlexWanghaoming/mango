# map raw reads to assembly
minimap2 -ax map-pb -t 30 yn56_polished.fasta /stor9000/apps/users/NWSUAF/2018055070/data/mo/YN56-1-1A.filtered_subreads.fq.gz | samtools view -@ 10 -bS - | samtools sort -@ 10 -o aln_sorted.bam -
