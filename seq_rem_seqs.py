#!/usr/bin/env python

import os, warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import sys, close

parser = argparse.ArgumentParser(description="Remove the trailing asterisk from amino acid sequences.")
parser.add_argument('f1', help="Input FASTA file")
parser.add_argument('f2', help="Input Textfile")
args = parser.parse_args()

rem_seqs = []
with open(args.f2) as file2:
	for l in file2:
		rem_seqs.append(l.strip())

out_seqs = []
for seq in list(SeqIO.parse(args.f1, "fasta")):
	if seq.description not in rem_seqs:
		out_seqs.append(seq)
	
SeqIO.write(out_seqs, args.f1+"_rem.fas", "fasta")


