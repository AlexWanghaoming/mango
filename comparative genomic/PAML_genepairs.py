from Bio.Phylo.PAML import yn00
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
import os, subprocess
import random
import pandas as pd

def seqAlign(fastaFile, outfile):
	muscle_cline = ClustalOmegaCommandline(infile=fastaFile, outfile=outfile, verbose=True, auto=True)
	print(muscle_cline)
	subprocess.run(str(muscle_cline), shell=True)

def runPAL2NAL(PEP_aln, CDS_aln, outfile):
	pal2nal_command = "~/software/pal2nal.v14/pal2nal.pl -nogap  -output paml %s %s > %s" \
	                  % (PEP_aln, CDS_aln, outfile)
	print(pal2nal_command)
	subprocess.run(pal2nal_command, shell=True)

def write_result(dic, wd, name):
	with open(wd + name, "a+") as res:
		for key in dic:
			res.write(key + "\t" + str(dic[key]) + "\n")

def generate_genepair(gene_path):
	a = pd.read_csv(gene_path, header=None)
	with open("/Users/alexwang/0data/0mango/synteny/all.anchors") as f2:
		genepair = [i.strip().split("\t") for i in f2.readlines()]
		random.shuffle(genepair)
	for gene1, gene2 in genepair:
		if not gene1 in list(a[0]) and (not gene2 in list(a[0])):
			yield gene1,gene2

class Results():
	omega = {}
	dn = {}
	ds = {}

def runYn00(genepairs, allCDSpath, allPEPpath, wd):
	print(Results.omega)

	gene1, gene2 = genepairs
	prefix = gene1 + "-" + gene2
	PEPfile = wd + prefix + ".pep"
	records1 = (r for r in SeqIO.parse(allPEPpath, "fasta") if r.id in genepairs)
	SeqIO.write(records1, PEPfile, "fasta")
	## align proteins sequence
	PEP_aln = wd + prefix + ".pep.aln"
	seqAlign(fastaFile=PEPfile, outfile=PEP_aln)

	# generate CDS sequence
	CDSfile = wd + prefix + ".cds"
	records2 = (r for r in SeqIO.parse(allCDSpath, "fasta") if r.id in genepairs)
	SeqIO.write(records2, CDSfile, "fasta")
	## align CDS using clustalo output .aln file
	CDS_aln = CDSfile + ".aln"
	seqAlign(fastaFile=CDSfile, outfile=CDS_aln)

	## generate .nuc file for yn00,
	nuc_file = wd + prefix + ".nuc"
	runPAL2NAL(PEP_aln=PEP_aln, CDS_aln=CDS_aln, outfile=nuc_file)

	## run yn00 in PAML
	yn00_res = wd + prefix + "_yn00.txt"
	yn = yn00.Yn00(alignment=nuc_file,
	               out_file=yn00_res,
	               working_dir=wd)
	yn.set_options(verbose=True)
	Yn00_results = yn.run(verbose=True)

	Results.omega["{0}-{1}".format(gene1, gene2)] = Yn00_results[gene1][gene2]["YN00"]["omega"]
	Results.dn["{0}-{1}".format(gene1, gene2)] = Yn00_results[gene1][gene2]["YN00"]["dN"]
	Results.ds["{0}-{1}".format(gene1, gene2)] = Yn00_results[gene1][gene2]["YN00"]["dS"]

def main():
	working_dir = "/Users/alexwang/0data/0mango/inf_associated_gene/secretome/paml/accessory/"
	CDS_PATH = "/Users/alexwang/0data/0mango/cds/all.cds"
	PEP_PATH = "/Users/alexwang/0data/0mango/EVM_protein/all.pep"

	size = 2232
	i = 0
	gene = "/Users/alexwang/0data/0mango/accessory/accessory.geneID"
	for genepairs in generate_genepair(gene_path=gene):
		i = i + 1
		if i<= size:
			runYn00(genepairs=genepairs, allCDSpath=CDS_PATH, allPEPpath=PEP_PATH, wd=working_dir)
		else:
			break
	write_result(dic=Results.omega, wd=working_dir, name="accessory_control_omega.txt")
	write_result(dic=Results.dn, wd=working_dir, name="accessory_control_dn.txt")
	write_result(dic=Results.ds, wd=working_dir, name="accessory_control_ds.txt")

if __name__ == '__main__':
	main()

