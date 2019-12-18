awk 'BEGIN{OFS="\t";print "ref","s","e","type"}$1~/[0-9]/{gsub("[a-zA-Z|]","",$5);print $5+1,$6,$7,$11}' gz15_polished.fasta.out > gz15.repeat.bed
