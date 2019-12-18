# normalize rawreads and contigs name
perl -pe 's/>[^\$]*$/">Seg" . ++$n ."\n"/ge' /stor9000/apps/users/NWSUAF/2018055070/data/mango/FJ11-1A.filtered_subreads.fasta > res/raw_reads.fasta
perl -pe 's/>[^\$]*$/">Seg" . ++$n ."\n"/ge' /stor9000/apps/users/NWSUAF/2018055070/data/mango/fj11/finisherSC/precheck/fj11_keep_contigs.fasta > res/contigs.fasta

python /stor9000/apps/users/NWSUAF/2018055070/sf/finishingTool-2.1/finisherSC.py -par 16 ./res /stor9000/apps/users/NWSUAF/2018055070/sf/mummer/bin
