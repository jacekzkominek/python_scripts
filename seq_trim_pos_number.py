#!/usr/bin/env python

from Bio import SeqIO
import sys, argparse

parser = argparse.ArgumentParser(description="Trim a number of columns from the start or end of alignment")
parser.add_argument('f1', help="Input FASTA file")
parser.add_argument('lencut', type=int, help="Number of columns to trim (negative=count from the end)")
args = parser.parse_args()

filename = args.f1
cut = int(args.lencut)
seqs = []

for seq_record in SeqIO.parse(filename, "fasta"):
	if (cut > 0):
		seq_trimmed = seq_record[cut:]
	else:
		seq_trimmed = seq_record[:cut]
	seqs.append(seq_trimmed)

SeqIO.write(seqs, filename+"_cut", "fasta")

