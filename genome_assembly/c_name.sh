#!/bin/bash

for i in `ls *.fasta`;do
    PRE=${i%_polish*}
    awk -v p="$PRE" '/>/{i++;gsub(/.*/,">"p"_"i,$0);print$0;next}{print$0}' ${i} > ${PRE}.fasta
done
