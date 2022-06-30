#!/usr/bin/env python

import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, argparse, os

parser = argparse.ArgumentParser(description="Filter stop codons (\"*\") from protein sequences.")
parser.add_argument('file1', help="Input file")
parser.add_argument('str', nargs="*", help="String to find")
args = parser.parse_args()

seqs = list(SeqIO.parse(open(args.file1), "fasta"))
out_seqs = []
for seq in seqs:
	for s in args.str:
		if seq.description.find(s) != -1: 
			out_seqs.append(seq)

for s in out_seqs:
	print ">"+s.description
	print s.seq

