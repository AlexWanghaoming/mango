from Bio.Phylo.PAML import yn00
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
import os, subprocess
import multiprocessing
import random

def split_list_n_list(origin_list, n):
	"""
	this function is used to split a list to n chunks for multi-process
	:param origin_list:
	:param n:
	:return: split list
	"""
	if len(origin_list) % n == 0:
		cnt = len(origin_list) // n
	else:
		cnt = len(origin_list) // n + 1

	for i in range(0, n):
		yield origin_list[i * cnt:(i + 1) * cnt]

def cdsAlign(CDSfile, outfile):
	muscle_cline = ClustalOmegaCommandline(infile=CDSfile, outfile=outfile, verbose=True, auto=True)
	print(muscle_cline)
	subprocess.run(str(muscle_cline), shell=True)

def runPAL2NAL(PEP_aln, CDS_aln, outfile):
	pal2nal_command = "~/software/pal2nal.v14/pal2nal.pl -nogap  -output paml %s %s > %s" \
	                  % (PEP_aln, CDS_aln, outfile)
	print(pal2nal_command)
	subprocess.run(pal2nal_command, shell=True)

def write_result(dic, wd, name):
	with open(wd + name, "a+") as res:
		# res.write("orthologues\tname\n")
		for key in dic:
			res.write(key + "\t" + str(dic[key]) + "\n")


def runYn00(ogFastaPath2, allCDSpath, wd):
	Omega = {}
	Dn = {}
	Ds = {}
	for ogFastaPath in ogFastaPath2:
		OG_fasta = os.path.basename(ogFastaPath)
		prefix = str(OG_fasta.split(".")[0])
		PEP_aln = wd + prefix + ".pep.aln"
		protein_id = []
		with open(PEP_aln, "w+") as f:
			for seq1 in SeqIO.parse(ogFastaPath, "fasta"):
				id = seq1.id.split("_")[2] + "_" + seq1.id.split("_")[3]
				protein_id.append(id)
				f.write(">" + id + "\n" + str(seq1.seq) + "\n")
		# generate CDS .fasta file
		CDSfile = wd + prefix + ".cds"
		records = (r for r in SeqIO.parse(allCDSpath, "fasta") if r.id in protein_id)
		SeqIO.write(records, CDSfile, "fasta")
		## align coding regions using clustalo output .aln file
		CDS_aln = CDSfile + ".aln"
		cdsAlign(CDSfile=CDSfile, outfile=CDS_aln)
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
		memo = []
		for gene1 in Yn00_results:
			memo.append(gene1)
			for gene2 in Yn00_results[gene1]:
				if not gene2 in memo:
					Omega["{0}-{1}".format(gene1, gene2)] = Yn00_results[gene1][gene2]["YN00"]["omega"]
					Dn["{0}-{1}".format(gene1, gene2)] = Yn00_results[gene1][gene2]["YN00"]["dN"]
					Ds["{0}-{1}".format(gene1, gene2)] = Yn00_results[gene1][gene2]["YN00"]["dS"]
	return Omega, Dn, Ds

def main():
	working_dir = "/Users/alexwang/0data/0mango/inf_associated_gene/secretome/paml/"
	CDS_PATH = "/Users/alexwang/0data/0mango/cds/all.cds"
	alignment_PATH = "/Users/alexwang/0data/0mango/1_OrthoFinder/Orthologues_Dec28/Alignments/"

	with open("/Users/alexwang/0data/0mango/1_OrthoFinder/SingleCopyOrthogroups.txt") as f1:
		singleCopyOG = [i.rstrip() for i in f1.readlines()]
	with open("/Users/alexwang/0data/0mango/1_OrthoFinder/CSEP_Orthogroups.txt") as f2:
	# with open("/Users/alexwang/0data/0mango/accessory/Orthogroups_accessory.txt") as f2:
		effectOG = [i.split(":")[0] for i in f2.readlines()]
	## orthogroups names of singleCopyEffectors
	singleCopyEffectorOG = list(set(singleCopyOG).intersection(set(effectOG)))
	## select random OGs for comparation
	size = len(singleCopyEffectorOG)
	print("random size: {0}".format(size))
	randomSingleCopyOG = random.sample(singleCopyOG, size)
	ll = []
	for PATH, _, fileList in os.walk(alignment_PATH):
		for fastaFileName in fileList:
			# if fastaFileName.split(".")[0] in singleCopyEffectorOG:
			if fastaFileName.split(".")[0] in randomSingleCopyOG:

				ll.append(alignment_PATH + fastaFileName)

	# run yn00 multi-process
	processor = 6
	pool = multiprocessing.Pool(processes=processor)
	a = []
	# split list into chunks
	OG_list_chunk = [k for k in split_list_n_list(ll, processor)]
	for j in range(processor):
		########## asynchronization， not blocking, not ordered results
		a.append(pool.apply_async(func=runYn00, args=(OG_list_chunk[j], CDS_PATH, working_dir)))
	res = [z.get() for z in a]
	########### blocking , ordered results
	# from functools import partial
	# res = pool.map(partial(runYn00, allCDSpath = CDS_PATH, wd = working_dir), ll)
	for tu in res:
		# output result
		write_result(dic=tu[0], wd=working_dir, name="CSEP_control_omega.txt")
		write_result(dic=tu[1], wd=working_dir, name="CSEP_control_dn.txt")
		write_result(dic=tu[2], wd=working_dir, name="CSEP_control_ds.txt")

if __name__ == '__main__':
	main()

