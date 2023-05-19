for i in hn47 fj11 qz yn56 gd10 gz15 yn55;do

	GENOME=/disk/alpha/mango/LTR/result/${i}.fasta
	#LTRharvest
	~/genometools-1.5.9/bin/gt suffixerator \
	  -db ${GENOME} \
	    -indexname ${i} \
		  -tis -suf -lcp -des -ssp -sds -dna
	~/genometools-1.5.9/bin/gt ltrharvest \
	  -index ${i} \
	    -similar 85 -vic 10 -seed 20 -seqids yes \
		  -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 > ${i}.harvest.scn
	# LTR harvest Non-canonical LTR-RT candidates
	~/genometools-1.5.9/bin/gt ltrharvest \
	 -index ${i}  \
	    -similar 85 -vic 10 -seed 20 -seqids yes \
		 -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6  -motif TGCA -motifmis 1 > ${i}.harvest.TGCA.scn

	# LTR_FINDER
	~/LTR_Finder/source/ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.85 ${GENOME} > ${i}.finder.scn

	# LTR_retriever integrate results of LTR_harvest and LTR_Finder
	source activate LTR_retriever
	~/LTR_retriever-master/LTR_retriever -genome ${GENOME} -inharvest ${i}.harvest.TGCA.scn -infinder ${i}.finder.scn -nonTGCA ${i}.harvest.scn -threads 2

done
