import re, sys

def parsingGFFattr(attr_string):
	attr_dict = {}
	attr_split = re.findall(r'.*? ".*?";', attr_string)
	for i in attr_split:
		if i:
			k, v = i.strip().split(" ",1)
			if v.endswith(";"):
				v = v[:-1]
			v = v.strip('"')
			attr_dict[k] = v
	return attr_dict

def getTpm(gtf, gene_dict):
	for fileIdx,file in enumerate(gtf):
		with open(file) as f:
			for line in f:
				if "TPM" in line and "reference_id" in line:
					*other, attr = line.strip().split("\t")
					attr_dict = parsingGFFattr(attr)
					geneid, tpm = attr_dict['ref_gene_id'], attr_dict['TPM']
					gene_dict[geneid][fileIdx] = float(tpm)
	return gene_dict

def main():
	# gtf = ["gd10_coni_2_stringTieOut.gtf", "gd10_inf5d_1_stringTieOut.gtf", "gd10_hypha_1_stringTieOut.gtf"]
	gtf = sys.argv[1:]
	stage = [i.split("_")[1] + i.split("_")[2] for i in gtf]
	# gene num of strain
	n = 17162
	# gene id prefix
	prefix = "CSIGD"
	gene_dict = {prefix + str(i).zfill(5): [0] * len(gtf) for i in range(1, n + 1)}
	res = getTpm(gtf, gene_dict)
	# outfile name
	with open("gd10_exp.txt", "w+") as outf:
		outf.writelines("geneID"+"\t" + "\t".join(stage) + "\n")
		for key in res:
			m = map(lambda x: str(x), res[key])
			outf.writelines(key+"\t" + "\t".join(m) + "\n")

if __name__ == '__main__':
	main()
