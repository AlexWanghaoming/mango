#!/bin/bash
set -e 
ref=/disk/alpha/mango/c.g/hn47.fasta
for i in `ls *.bam`;do
	## denovo assemble transcript
	stringtie ${i} -p 10 -o ${i%.sorted.bam}.stingtie.assembly.gtf
done

#MASKED_GENOME=/disk/alpha/mango/masked.genome/hn47_upper.fasta.masked
#ortho_protein=/disk/alpha/mango/odb9_fungi/odb9_fungi.fasta

# 1. Augustus Abinitio prediction using braker2 integrating evidence from protein,RNA-seq
#braker.pl --fungus --cores 14 --etpmode --softmasking --gff3 --genome=${MASKED_GENOME} --prot_seq=${ortho_protein} --bam=RNAseq/67-1_24a.sorted.bam,RNAseq/67-1_48a.sorted.bam,RNAseq/67-1_8a.sorted.bam,RNAseq/C-24-1.sorted.bam,RNAseq/S-24-1.sorted.bam,RNAseq/C-8-1.sorted.bam,RNAseq/S-8-1.sorted.bam,RNAseq/DON_rep2.sorted.bam,RNAseq/MeOH_rep1.sorted.bam,RNAseq/PG_rep1.sorted.bam

# 2. GeneMark_ES Abinitio prediction
python3 ~/bin/gtfUtils.py -i retrain/braker/GeneMark-ES/genemark.gtf -reformat > genemark_ES.gtf
awk 'BEGIN{FS=OFS="\t"}{gsub("\"","",$9);gsub("GeneMark.hmm","GeneMark_ES",$2);if($3~/gene/){gsub("gene_id ","ID=",$9)}else{if($3~/transcript/){gsub("transcript","mRNA", $3);gsub("gene_id ","Parent=",$9);gsub("transcript_id ","ID=",$9)}else{i=i+1;gsub("gene_id ","ID=",$9);gsub("_g","_g"i,$9)gsub("transcript_id ", "Parent=", $9)}};print}' genemark_ES.gtf > genemark_ES.gff3

# 3. extract candicate coding regions from denovo transcripts assembly
# merge transcript assemblies
stringtie --merge -o stringtie_assembly.merge.gtf -p 8 *stingtie.assembly.gtf
trans_asm=/disk/alpha/mango/HN47-2A/rnaseq/stringtie_assembly.merge.gtf
~/EvidenceModeler/TransDecoder-TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl ${trans_asm} ${ref} > transcripts.fasta
/home/wanghm/EvidenceModeler/TransDecoder-TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl ${trans_asm} > transcripts.gff3
/home/wanghm/EvidenceModeler/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t transcripts.fasta
/home/wanghm/EvidenceModeler/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t transcripts.fasta
/home/wanghm/EvidenceModeler/TransDecoder-TransDecoder-v5.5.0/util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3

########### run EVMmodeler using transDecoder output; braker2 output;geneMark_ES output as gene_predictions.gff3
##########                using geneMark_ET output as transcript_alignments.gff3

cat transcripts.fasta.transdecoder.genome.gff3 retrain/braker/braker.gff3 genemark_ES.gff3 > gene_predictions.gff3
#ln gene_predict_STAR_hints/genemark_ET.gff3 transcript_alignments.gff3
~/EvidenceModeler/EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome ${ref} --gene_predictions gene_predictions.gff3 \
		--segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
~/EvidenceModeler/EVidenceModeler-1.1.1/EvmUtils/write_EVM_commands.pl --genome ${ref} --weights /disk/alpha/mango/HN47-2A/rnaseq/weights.txt --gene_predictions gene_predictions.gff3 \
				      --output_file_name evm.out --partitions partitions_list.out >  commands.list
~/EvidenceModeler/EVidenceModeler-1.1.1/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log
~/EvidenceModeler/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
~/EvidenceModeler/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome ${ref}
find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff3

# reformat GTF 
gffread -T EVM.all.gff3 -o EVM.all.gtf
python3 ~/bin/gtfUtils.py -i EVM.all.gtf -reformat > EVM.all.reformat.gtf
python3 ~/bin/gtfUtils.py -i EVM.all.reformat.gtf -r CFRHN_ -o hn47_EVM.gtf
gffread -y hn47_EVM.protein.fasta -T hn47_EVM.gtf -g ${ref}
awk '/>/{split($2,a,"=");print ">"a[2];next}{gsub("\\.","",$0);print}' hn47_EVM.protein.fasta > hn47_protein.fasta
rm hn47_EVM.protein.fasta
