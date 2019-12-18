# nucmer --maxgap=500 --mincluster=100 --prefix=gz15_qz ../gz15.fasta ../qz.fasta
nucmer --maxgap=500 --mincluster=100 --prefix=gz15_yn55 ../gz15.fasta ../yn55.fasta
nucmer --maxgap=500 --mincluster=100 --prefix=gz15_yn56 ../gz15.fasta ../yn56.fasta
nucmer --maxgap=500 --mincluster=100 --prefix=gz15_hn47 ../gz15.fasta ../hn47.fasta
nucmer --maxgap=500 --mincluster=100 --prefix=gz15_fj11 ../gz15.fasta ../fj11.fasta
nucmer --maxgap=500 --mincluster=100 --prefix=gz15_gd10 ../gz15.fasta ../gd10.fasta

# output table format
# show-coords -d -l -L 10 -r -c -T gz15_qz.delta > gz15_qz.tab.txt
show-coords -d -l -L 10 -r -c -T gz15_yn55.delta > gz15_yn55.tab.txt
show-coords -d -l -L 10 -r -c -T gz15_yn56.delta > gz15_yn56.tab.txt
show-coords -d -l -L 10 -r -c -T gz15_hn47.delta > gz15_hn47.tab.txt
show-coords -d -l -L 10 -r -c -T gz15_fj11.delta > gz15_fj11.tab.txt
show-coords -d -l -L 10 -r -c -T gz15_gd10.delta > gz15_gd10.tab.txt
