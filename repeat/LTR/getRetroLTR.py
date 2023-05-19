
from Bio import SeqIO
import os, time
import sys

def openFasta(faFile:str) -> tuple:
	if not os.path.exists(faFile):
		raise IOError('Fasta file does not exists: %s' % faFile)
	fa = SeqIO.parse(faFile,format="fasta")
	seq_dict = {}
	for seq in fa:
		seq_dict[seq.name] = seq.seq
	return faFile,seq_dict

def getSeq(seq_dict, chr, s, e):
	chrom = header + "_" + str(int(chr)+1)
	seq = seq_dict[chrom]
	seq = seq[int(s):(int(e)-1)]
	return seq
# fa = SeqIO.parse("/Users/alexwang/0data/0mango/oth_colletotrichum/c.truncatum/c.truncatum.fa", format="fasta")

header = sys.argv[1]
#_, seq_dict = openFasta("/Users/alexwang/0data/0mango/oth_colletotrichum/c.truncatum/c.truncatum.fa")
_, seq_dict = openFasta(sys.argv[2])
i = 0
#with open("/Users/alexwang/0data/0mango/genome/LTR/c.truncatum.harvest.scn") as f:
outf = "/Users/alexwang/0data/0mango/genome/LTR/retroele_{}.fa".format(header)
ff = open(outf, "w+")
with open(sys.argv[3]) as f:
	for idx, line in enumerate(f):
		if not line.startswith("#"):
			i = i + 1
			s_ret, e_ret, l_ret, s_lLTR, e_lLTR, l_lLTR, s_rLTR, e_rLTR, l_rLTR, sim, seq_id = line.split()
			inter_seq = getSeq(seq_dict, seq_id, e_lLTR, s_rLTR)
			ff.writelines(">"+header+ "_LTR" + str(i)+"\n"+str(inter_seq)+"\n")

#			outf2 = "/Users/alexwang/0data/0mango/genome/LTR/phylip_test/clustalOutfile%d.fa" % (idx-11)
#			clustalomega_cline = ClustalOmegaCommandline(infile=outf , outfile=outf2, outfmt='phylip', verbose=True, auto=False)
#			clustalomega_cline()
			# time.sleep(100)
