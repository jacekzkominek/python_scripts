#!/usr/bin/env python

import os, warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import sys, close

parser = argparse.ArgumentParser(description="Remove the trailing asterisk from amino acid sequences.")
parser.add_argument('f1', help="Input FASTA file")
args = parser.parse_args()

f = open(args.f1)
# ~ if args.f1 == "all":
	# ~ for f in os.listdir("."):
out_seqs = []
passed = set()
# ~ if f.find("full.pep") != -1:
for seq in SeqIO.parse(f, "fasta"):
	if str(seq.seq) not in passed:
		out_seqs.append(seq)
		passed.add(str(seq.seq))
SeqIO.write(out_seqs, args.f1+"_dedup.fas", "fasta")


