#!/usr/bin/env python

import warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import sys, close

parser = argparse.ArgumentParser(description="Sort sequences based on length (ascending by default)")
parser.add_argument('f1', help="Input FASTA file")
parser.add_argument('--reverse', action='store_true', default=False, help="Sort descending")
parser.add_argument('--extract', type=int, default=0, help="Sort descending")
args = parser.parse_args()

f = open(args.f1)
seqs = list(SeqIO.parse(f, "fasta"))

out_seqs = [] 
if args.reverse == False:
	for seq in sorted (seqs, key=lambda s: len(s)):
		if len(seq) >= args.extract:
			out_seqs.append(seq)
elif args.reverse == True:
	for seq in sorted (seqs, key=lambda s: len(s), reverse=True):
		if len(seq) >= args.extract:
			out_seqs.append(seq)

if args.extract == 0:
	SeqIO.write(out_seqs, args.f1.replace(".fas", "_lensort.fas"), "fasta")
else:
	SeqIO.write(out_seqs, args.f1.replace(".fas", "_lensort"+str(args.extract)+".fas"), "fasta")

f.close()

