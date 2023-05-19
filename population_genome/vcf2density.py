import sys
"""
para1: genome.txt 2 columns
para2: vcf file
para3: output file
"""
# gs = open("/Users/alexwang/0data/0mango/genome/yn55.genome.txt").readlines()
gs = open(sys.argv[1]).readlines()
bs = 25000

def makeWindows(kk, bs=bs):
	id,size = kk.split()
	size = int(size)
	a = []
	start = 1
	while start<size:
		end = start+bs-1
		if end + bs >= size:
			a.append([start, size, 0])
			break
		else:
			a.append([start,end,0])
		start = start+bs
	return a

def main():
	tt = list(map(makeWindows, gs[1:]))
	with open(sys.argv[2]) as f:
		for rec in f:
			for i in gs[1:]:
				chr_id = int(i.split()[0])
				if not rec.startswith("#"):
					chrom, pos, *remain = rec.split()
					if chrom == "yn55_"+str(chr_id):
						idx = int(pos)//bs
						try:
							tt[chr_id-1][idx][2]+=1
						except IndexError:
							tt[chr_id-1][-1][2]+=1
	# write records
	res = open(sys.argv[3],"w+")
	for i in range(len(tt)):
		chr_group = tt[i]
		for rec in chr_group:
			rec.insert(0, "yn55_" + str(i+1))
			res.writelines("\t".join(map(str,rec)) + "\n")
	res.close()

if __name__ == '__main__':
	# tt = list(map(makeWindows, gs[1:]))
	# print(tt)
	main()




