#!/usr/bin/env python

import warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import sys, close

parser = argparse.ArgumentParser(description="Trim sequences with identical names.")
parser.add_argument('f1', help="Input FASTA file")
args = parser.parse_args()

f = open(args.f1)
seqs = list(SeqIO.parse(f, "fasta"))
out_seqs = []
out_seqs_names = set()

for seq in seqs:
	if seq.description in out_seqs_names:
		continue
	else:
		out_seqs_names.add(seq.description)
		out_seqs.append(seq)

f.close()

SeqIO.write(out_seqs, args.f1+"a", "fasta")
