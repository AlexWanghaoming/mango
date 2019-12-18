import sys
total_assembly_length = 0

threshold = 50
threshold = (threshold*0.01)
fa = sys.argv[1]
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

with open(fa) as fp:
    for name, seq in read_fasta(fp):
        total_assembly_length += len(seq)

with open(fa) as fp:
    x = [(len(seq), name, seq) for name, seq in read_fasta(fp)]
    x = sorted(x, reverse=True)
    count = 0
    for length, name, seq in x:
        count += (length)
        if count >= total_assembly_length*(threshold):
            print("N50: %s bp (%s)" % (length, name))
            break
