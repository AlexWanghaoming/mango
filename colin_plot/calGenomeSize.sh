awk 'BEGIN{i=1}/>/&&NR>1{print ""}{printf "%s",/^>/ ? i++" " : $0}' yn56.fasta | awk 'BEGIN{print "chr\tsize"}{print$1"\t"length($2)}' > genome/yn56.genome.txt
