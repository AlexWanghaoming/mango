import time
import sys,re
import argparse
from collections import defaultdict
def main(args):
	## Add GO annotation
	#with open("/Users/alexwang/0data/0mango/protein/interproscan/go_annot.txt") as f:
	with open(args.go) as g:
		go_list = g.readlines()
	N=14
	sub_OG = []
	if not args.go:
		# with open("/Users/alexwang/0data/0mango/transcriptome/ortholog_heat/cazy.geneID") as target:
		with open(args.gene) as target:
			target_geneID = list(target.read().splitlines())
		outf1 = open(args.pre + "_count.tsv", "w+")
		outf2 = open(args.pre + "_targetGene.tsv", "w+")
	else:
		outf = open("orthogroups_GO.tsv", "w+")
	og_go = defaultdict(list)
	# with open("/Users/alexwang/0data/0mango/transcriptome/ortholog_heat/Orthogroups.tsv") as f:
	with open(args.ortho) as f:
		for idx, line in enumerate(f):
			if not args.go:
				og = line.strip().split("\t")
				cc = [og[0]]
				gg = [og[0]]
				for strain_og in og[1:]:
					sub_gg = []
					target_count = 0
					gene_list = strain_og.strip().split(", ")
					for ge in gene_list:
						if ge in target_geneID:
							sub_gg.append(ge)
							target_count += 1
					cc.append(target_count)
					if sub_gg:
						gg.append(",".join(sub_gg))
				if sum(cc[1:]) != 0:
					sub_OG.append(line)
					if len(cc)!=N+1:
						cc = cc + [0]*(N+1-len(cc))
					cc = list(map(str, cc))
					outf1.writelines("\t".join(cc)+"\n")
					outf2.writelines("\t".join(gg)+"\n")
			else:
				og = re.split("[\t ,]", line.strip())
				for k in go_list:
					gene, goid, definition, ontology, term = k.strip().split("\t")
					if gene in og[1:]:
						if not term in og_go[og[0]]:
							og_go[og[0]].append(term)

	if not args.go:
		with open(args.pre + "_subOG.tsv", "w+") as outf3:
			outf3.writelines(sub_OG)
		outf2.close()
		outf1.close()
	else:
		for key in og_go:
			outf.writelines(key+"\t"+"|".join(og_go[key]) + "\n")


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-og", "--orthogroup", help="input orthogroup tsv from orthoFinder2", dest="ortho")
	parser.add_argument("-target", "--target", help="input target gene id", dest="gene")
	parser.add_argument("-pre", "--prefix", help="input prefix of target genes", dest="pre")
	parser.add_argument("-go", "--GOannotation", help="add GO term to orthogroups", dest="go")
	# parser.add_argument("-s", "--sub", help="output sub orthogroup of target genes", dest="sub", action="store_true")
	args = parser.parse_args()
	main(args)




