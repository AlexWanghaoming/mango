#!/bin/bash

# remove big file to save storage
cd $mecat_path
find ./ -type f -size +500M -not -name "*.fasta" |xargs rm

# make bam.fofn
find ./ -name "*.bam" |xargs -n 1 -I {} readlink -f {} > bam.fofn
