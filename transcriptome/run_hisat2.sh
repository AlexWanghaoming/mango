# strand-specific RNAseq
strain=gd10
array=(S-1-1_FRRB192000669-1a S-1-2_FRRB192000670-1a S-1-3_FRRB192000671-1a)
length=${#array[@]}
for ((i=0;i<${length};i++))
do
hisat2 -x /stor9000/apps/users/NWSUAF/2018055070/data/mango/rna_seq/index/${strain} -p 24 --rna-strandness RF --min-intronlen 20 --max-intronlen 4000 --no-unal --no-temp-splicesite --novel-splicesite-outfile ${strain}_coni_$((${i}+1)).splicesite --novel-splicesite-infile ${strain}_coni_$((${i}+1)).splicesite -1 ${array[$i]}_1.clean.fq.gz -2 ${array[$i]}_2.clean.fq.gz | samtools view -bS - | samtools sort -o ${strain}_coni_$((${i}+1)).sorted.bam - 
done

array=(S-3-1_FRRB192000675-1a S-3-2_FRRB192000676-1a S-3-3_FRRB192000677-1a)
length=${#array[@]}
for ((i=0;i<${length};i++))
do
hisat2 -x /stor9000/apps/users/NWSUAF/2018055070/data/mango/rna_seq/index/${strain} -p 24 --rna-strandness RF --min-intronlen 20 --max-intronlen 4000 --no-unal --no-temp-splicesite --novel-splicesite-outfile ${strain}_hypha_$((${i}+1)).splicesite --novel-splicesite-infile ${strain}_hypha_$((${i}+1)).splicesite -1 ${array[$i]}_1.clean.fq.gz -2 ${array[$i]}_2.clean.fq.gz | samtools view -bS - | samtools sort -o ${strain}_hypha_$((${i}+1)).sorted.bam - 
done

array=(MG14-1_FRRB190107869-1a MG14-2_FRRB190107870-1a MG14-3_FRRB190107871-1a)
length=${#array[@]}
for ((i=0;i<${length};i++))
do
hisat2 -x /stor9000/apps/users/NWSUAF/2018055070/data/mango/rna_seq/index/${strain} -p 24 --rna-strandness RF --min-intronlen 20 --max-intronlen 4000 --no-unal --no-temp-splicesite --novel-splicesite-outfile ${strain}_inf3d_$((${i}+1)).splicesite --novel-splicesite-infile ${strain}_inf3d_$((${i}+1)).splicesite -1 ${array[$i]}_1.clean.fq.gz -2 ${array[$i]}_2.clean.fq.gz | samtools view -bS - | samtools sort -o ${strain}_inf3d_$((${i}+1)).sorted.bam - 
done

array=(MG16-1_FRRB190107875-1a MG16-2_FRRB190107876-1a MG16-3_FRRB190107877-1a)
length=${#array[@]}
for ((i=0;i<${length};i++))
do
hisat2 -x /stor9000/apps/users/NWSUAF/2018055070/data/mango/rna_seq/index/${strain} -p 24 --rna-strandness RF --min-intronlen 20 --max-intronlen 4000 --no-unal --no-temp-splicesite --novel-splicesite-outfile ${strain}_inf5d_$((${i}+1)).splicesite --novel-splicesite-infile ${strain}_inf5d_$((${i}+1)).splicesite -1 ${array[$i]}_1.clean.fq.gz -2 ${array[$i]}_2.clean.fq.gz | samtools view -bS - | samtools sort -o ${strain}_inf5d_$((${i}+1)).sorted.bam - 
done
