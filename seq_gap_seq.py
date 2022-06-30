#!/usr/bin/env python

import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import sys, close
from collections import defaultdict

f = open(sys.argv[1])
seqs = list(SeqIO.parse(f, "fasta"))

seq_gaps = defaultdict(int)
seq_len = len(seqs)
aln_len = len(seqs[0])

for seq in seqs:
	for p in range(aln_len):
		if seq[p] == "-":
			seq_gaps[seq.description] += 1
f.close()

new_f = open(sys.argv[1].replace(".fas", "_gap_seq.txt"),'w')

for seq in seqs:
	print>>new_f, seq.description+","+str(seq_gaps[seq.description])

new_f.close()

