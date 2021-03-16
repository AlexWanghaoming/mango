import sys
import re
import time
# convert RepBase data to fasta
# embl = "/Users/alexwang/0data/db/Libraries/RMRBSeqs.embl"
embl = sys.argv[1]
prefix = re.split(r"[.]", embl)[-2]
out = open(prefix+".fasta", "w")

cc = ""
with open("/Users/alexwang/0data/db/Libraries/RMRBSeqs.embl") as ff:
	for line in ff:
		if line.startswith("ID"):
			info = line.split()
			out.write(">" + info[1] + "\n")
		elif line.startswith(" "):
			strr = "".join(re.findall("[atcgnATCGN]", line))
			cc += strr
		elif line.startswith("//"):
			out.write(cc + "\n")
			cc = ""
		else:
			pass
out.close()
