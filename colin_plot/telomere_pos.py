#!/usr/local/bin/python3
import re,sys
from Bio import SeqIO


def seq_info(fasta_path):
	with open(fasta_path,"r") as f:
		l = SeqIO.parse(f,"fasta")
		for seq in l:
			yield seq


def revcom_seq(seq):
	seq = seq.upper()
	intab = "ATCGN"
	outtab = "TAGCN"
	trantab = str.maketrans(intab,outtab)
	return seq[::-1].translate(trantab)


def find_telo(path, *telo_string):   # 用星号传入元祖
	string = [revcom_seq(i) for i in telo_string] + list(telo_string)
	res = {}
	for seq in seq_info(path):
		na = int(seq.name.split("_")[1])
		sequence = str(seq.seq).upper()
		seq_len = len(seq.seq)
		l, r = float('inf'), float('inf')
		for telo in string:
			l_search = re.search(telo, sequence[:max_telo_distance])
			r_search = re.search(telo[::-1], sequence[:-max_telo_distance:-1])
			if l_search:
				l = min(l, l_search.start()+1)
			if r_search:
				r = min(r, r_search.start())
		res[na] = (seq_len, l if l!=float('inf') else 'nd', seq_len-r if r!=float('inf') else "nd")
	return res


if __name__ == '__main__':
	fa = sys.argv[1]
	telo_str1 = "CCCTAA"*2
	# telo_str2 =
	max_telo_distance = 50000
	res = find_telo(fa, telo_str1)
	with open("%s.telo.txt" % sys.argv[1].split(".")[0],"w") as final:
		final.writelines("\t".join(["ref","chr_len", "left", "right"])+"\n")
		final.writelines(str(key)+"\t"+str(value[0])+"\t"+str(value[1])+"\t"+str(value[2])+"\n" for key,value in res.items())

