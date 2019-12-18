#!/bin/bash
#cat `find ~/data/qz-3/ -name "*.fasta"` > qz_merge.fasta
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/stor9000/apps/users/NWSUAF/2018055070/sf/glibc_2.14/lib

canu -correct -p qz -d qz_canu_assembly_corMhap genomeSize=60m useGrid=false corMhapSensitivity=normal minReadLength=2000 minOverlapLength=500 -pacbio-raw /stor9000/apps/users/NWSUAF/2018055070/data/qz_merge.fasta

canu -trim \
        -p qz -d qz_canu_assembly_corMhap maxThreads=24 genomeSize=60m minReadLength=2000 minOverlapLength=500\
        -pacbio-corrected /stor9000/apps/users/NWSUAF/2018055070/data/qz.correctedReads.fasta.gz

# if coverage > 90 set correctedErrorRate=0.013
canu -assemble \
        -p qz -d qz_canu_assembly_corMhap maxThreads=16 genomeSize=60m correctedErrorRate=0.035\
        redMemory=16 oeaMemory=16 -pacbio-corrected /stor9000/apps/users/NWSUAF/2018055070/data/qz_canu_assembly_corMhap/qz.trimmedReads.fasta.gz
