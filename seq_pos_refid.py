#!/usr/bin/env python

import warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from os import sys, close

parser = argparse.ArgumentParser(description="I have no clue")
parser.add_argument('f1', help="Input FASTA file")
args = parser.parse_args()

f = open(args.f1)
seqs = list(SeqIO.parse(f, "fasta"))
seq_count = len(seqs)
seq_len = len(seqs[0])

l = open(args.f1+"_scores.txt",'w')
for x in range (0, seq_len):
	ref_res = seqs[0][x]
	score = 1	
	for i, seq in enumerate(seqs):
			if i > 0 and seqs[i][x] == ref_res:
				score += 1
	print>>l,"Position",x+1,ref_res, score, str(score*100/seq_count)

f.close()
l.close()
print


