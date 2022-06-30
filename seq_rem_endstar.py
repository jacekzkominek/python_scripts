#!/usr/bin/env python

import warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import sys, close

parser = argparse.ArgumentParser(description="Remove the trailing asterisk from amino acid sequences.")
parser.add_argument('f1', help="Input FASTA file")
args = parser.parse_args()

f1 = open(args.f1)
seqs = list(SeqIO.parse(f1, "fasta"))
out_seqs = []
print args.f1
for seq in seqs:
	if seq.seq[-1] == "*":
		seq.seq = seq.seq[:-1]
	out_seqs.append(seq)
	
SeqIO.write(out_seqs, args.f1, "fasta")

f1.close()
