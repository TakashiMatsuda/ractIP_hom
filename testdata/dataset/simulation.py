from Bio import SeqIO

short_sequences = []

sequence_name = args[0]

for record in SeqIO.parse(open(sequence_name, "rU"), "clustal") :
    if len(record.seq) < 300 :
        short_sequences.append(record)

print "Found %i sequences" % len(short_sequences)

output_name = args[1]
output_handle = open((output_name), "wU"), "fasta")
####　ないファイル名に対して書き込みを行えるか、Pythonのopen関数にそれが実装されているか

### random anealing
for record in 
