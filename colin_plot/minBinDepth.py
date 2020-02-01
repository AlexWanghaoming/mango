from itertools import islice
import sys

bs = 5000

with open(sys.argv[1]) as f:
	repo = {}
	for li in f:
		chrom,pos,cov = li.split()
		pos,cov = int(pos),int(cov)
		if not chrom in repo.keys():
			repo[chrom] = []
			if pos != 1:
				for i in range(1, pos):
					repo[chrom].append([i,0])
		repo[chrom].append([pos, cov])

	res = open("%s.txt" % sys.argv[2],"w+")
	for key in repo.keys():
		ll = iter(repo[key])    # islice对迭代器切片
		while True:
			bin = list(islice(ll, bs))     # 每次对迭代器消耗5000个值
			if bin:
				s,e,me,ma,mi = bin[0][0],\
				               bin[-1][0],\
				               sum([i[1] for i in bin])/bs,\
				               max(i[1] for i in bin),\
				               min(i[1] for i in bin)
				res.writelines("\t".join(map(str,[key,s,e,me,ma,mi]))+"\n")
			else:
				break
	res.close()




