#!/usr/bin/env python

import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, argparse, os

parser = argparse.ArgumentParser(description="Filter stop codons (\"*\") from protein sequences.")
parser.add_argument('file1', help="Input file")
parser.add_argument('--from_end', action='store_true', default=False, help="Strip from the end of sequences.")
#parser.add_argument('--all_seq', action='store_true', default=False, help="Strip from the entire sequences.")
parser.add_argument('--sep', action='store_true', default=False, help="Separate sequences with mid-sequence stop codons.")
args = parser.parse_args()

#if args.all_seq == True:
#	args.end_seq = False

f = open(args.file1)
seqs = list(SeqIO.parse(f, "fasta"))

nonX_seqs = []
X_seqs = []



for seq in seqs:
	if args.from_end == True:
		if seq.seq.find("*") == len(seq.seq)-1:
			X_seqs.append(SeqRecord(Seq(str(seq.seq[:-1]), IUPAC.protein),id=seq.id,name=seq.name,description=seq.description))
		else:
			X_seqs.append(seq)
	elif args.sep == True:
		if seq.seq.find("*") != -1:
			X_seqs.append(seq)
		else:
			nonX_seqs.append(seq)

f.close()

#for (i,a) in enumerate(X_seqs):
	#print i,a

if args.from_end == True:
	SeqIO.write(X_seqs, args.file1.replace(".fas", "_endX.fas"), "fasta")

elif args.sep == True:
	SeqIO.write(X_seqs, args.file1.replace(".fas", "_X.fas"), "fasta")
	SeqIO.write(nonX_seqs, args.file1.replace(".fas", "_nonX.fas"), "fasta")


