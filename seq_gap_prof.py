#!/usr/bin/env python

import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import sys, close
from collections import defaultdict

f = open(sys.argv[1])
seqs = list(SeqIO.parse(f, "fasta"))

pos_gaps = defaultdict(int)
pos_gaps_count = defaultdict(int)
aln_len = len(seqs[0].seq)

for p in range(aln_len):
	gaps = 0
	for seq in seqs:
		if seq[p] == "-":
			gaps += 1
	pos_gaps[p] = gaps
	pos_gaps_count[gaps] += 1

		
f.close()

new_f = open(sys.argv[1].replace(".fas", "_gap_prof.txt"),'w')


for gc in sorted(pos_gaps_count.keys()):
	print>>new_f, gc, pos_gaps_count[gc]

print>>new_f, "\n"
for p in range(aln_len):
	print>>new_f, int(p+1), pos_gaps[p]

new_f.close()

