#!/bin/bash
pbmm2 align -j 20 /stor9000/apps/users/NWSUAF/2018055070/data/qz_falcon/finisherSC/res/improved3.fasta bam.fofn | samtools sort > qz_sorted.bam
pbindex qz_sorted.bam
samtools faidx /stor9000/apps/users/NWSUAF/2018055070/data/qz_falcon/finisherSC/res/improved3.fasta
arrow qz_sorted.bam -r /stor9000/apps/users/NWSUAF/2018055070/data/qz_falcon/finisherSC/res/improved3.fasta -o qz_polished.fasta -j 20

