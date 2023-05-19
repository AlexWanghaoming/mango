#!/bin/bash
BuildDatabase -name yn56_db -engine ncbi ./yn56_polished.fasta
RepeatModeler -engine ncbi -pa 4 -database yn56_db
~/software/RepeatMasker/util/queryRepeatDatabase.pl -species fungi > ./fungi_repeats.lib
cat fungi_repeats.lib RM*/consensi.fa.classified > ./myGenome.custom.repeat.lib

RepeatMasker -engine crossmatch -lib myGenome.custom.repeat.lib -dir yn56_rm -pa 4 -no_is yn56_polished.fasta
## for gene prediction
# RepeatMasker --engine ncbi -xsmall -nolow -lib myGenome.custom.repeat.lib -dir yn56_ncbiSoftLTR -pa 5 -no_is yn56.fasta
