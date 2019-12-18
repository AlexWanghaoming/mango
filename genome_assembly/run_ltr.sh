#LTRharvest
~/genometools-1.5.9/bin/gt suffixerator \
	  -db yn56_polished.fasta \
	    -indexname yn56 \
		  -tis -suf -lcp -des -ssp -sds -dna
~/genometools-1.5.9/bin/gt ltrharvest \
	  -index yn56 \
	    -similar 85 -vic 60 -seed 30 -seqids yes \
		  -minlenltr 100 -maxlenltr 1000 -mintsd 4 -maxtsd 20 > yn56.harvest.scn
# LTR_FINDER
#~/LTR_Finder/source/ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C yn56_polished.fasta > yn56.finder.scn

# LTR_retriever
#~/LTR_retriever-master/LTR_retriever -genome yn56_polished.fasta -inharvest yn56.harvest.scn -infinder yn56.finder.scn -threads 10
