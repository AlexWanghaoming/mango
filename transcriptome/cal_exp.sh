#!/bin/bash
strain=yn56
# run featureCounts
/stor9000/apps/users/NWSUAF/2018055070/sf/subread-2.0.1-Linux-x86_64/bin/featureCounts -t CDS -p -s 2 -g gene_id -T 20 -a /stor9000/apps/users/NWSUAF/2018055070/data/mango/GTF/${strain}_primaryTrans.gtf -o ${strain}_readsCount.raw.tsv *.bam
# run stringTie
for i in `ls *.sorted.bam`;do
	~/sf/stringtie-2.1.3b.Linux_x86_64/stringtie ${i} -G /stor9000/apps/users/NWSUAF/2018055070/data/mango/GTF/${strain}_primaryTrans.gtf -o ${i%.sorted.bam}_stringTieOut.gtf --rf -p 20
done
## manage stringTie output
python3 getStringTieTPM.py *_stringTieOut.gtf
