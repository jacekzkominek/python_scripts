#!/usr/bin/env python

import warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from os import sys, remove, close

parser = argparse.ArgumentParser(description="Trim alignment positions with identical codons.")
parser.add_argument('f1', help="Input FASTA file")
args = parser.parse_args()

f = open(args.f1)
seqs = list(SeqIO.parse(f, "fasta"))
seq_count = len(seqs)
seq_len = len(seqs[0])
f.seek(0)

trim_singletons = 0
nonid_pos = []
					
for x in range (0, seq_len/3):
	first_codon = ""
	for i, seq in enumerate(list(SeqIO.parse(f, "fasta"))):
		if i == 0:
			first_codon = seq[x*3]+seq[x*3+1]+seq[x*3+2]
			continue
		else:
			current_codon = seq[x*3]+seq[x*3+1]+seq[x*3+2]
			if first_codon != current_codon:
				nonid_pos.append(x)
				f.seek(0)
				break
			if i == seq_count-1:
				f.seek(0)

trimmed = []
f.seek(0)
for seq in list(SeqIO.parse(f, "fasta")):
	seq1 = ""
	for x in nonid_pos:
		seq1 += seq[x*3]+seq[x*3+1]+seq[x*3+2]
	trimmed.append(seq1)

f.seek(0)

new_f = open(args.f1.replace(".fas", "_trim_IDcod.fas"),'w')
	
count = 0
for seq in list(SeqIO.parse(f, "fasta")):
	new_f.write(">"+seq.id+"\n")
	new_f.write(trimmed[count]+"\n")
	count += 1
new_f.close()
f.close()

