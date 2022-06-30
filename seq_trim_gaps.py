#!/usr/bin/env python

import warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from os import sys, remove, close

parser = argparse.ArgumentParser(description="Trim alignment positions with gaps.")
parser.add_argument('f1', help="Input FASTA file")
args = parser.parse_args()

f = open(args.f1)
gap_pos = set()
for seq in list(SeqIO.parse(f, "fasta")):
	for x in range (0, len(seq)):
		if seq[x].find("-") != -1:
			gap_pos.add(x)

seq_pos = set()
for x in range (0, len(seq)):
	seq_pos.add(x)

trim_pos = sorted(seq_pos.difference(gap_pos))
degapped = []

f.seek(0)
for seq in list(SeqIO.parse(f, "fasta")):
	seq1 = ""
	for x in trim_pos:
		seq1 += seq[x]
	degapped.append(seq1)

f.seek(0)
new_f = open(args.f1.replace(".fas", "_trim.fas"),'w')
count = 0
for seq in list(SeqIO.parse(f, "fasta")):
	new_f.write(">"+seq.description+"\n")
	new_f.write(degapped[count]+"\n")
	count += 1
new_f.close()
f.close()
