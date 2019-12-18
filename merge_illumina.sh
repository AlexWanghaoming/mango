#!/bin/bash
# used for merge all latest and previous Illumina fq.gz data of strains.
for id in {1..48};do
    echo $id
    PATH1=/home/sad7_disk/data/F19FTSSCWLJ0595_METiczR/Upload_F19FTSSCWLJ0595_METiczR/BGI_result/Separate
    PATH2=/home/wanghm/Downloads
    for i in 1 2;do
        echo "strain $id-$i"
        a=`find $PATH1 -type f -name "${id}.*${i}.fq.gz"` 
        b=`find $PATH2 -type f -name "${id}[^0-9]*${i}.fq.gz"`
        cat ${a} ${b} > ${id}_${i}_merged.fq.gz 
    done
done
