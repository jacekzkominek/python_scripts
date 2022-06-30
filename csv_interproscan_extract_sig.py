#!/usr/bin/env python

import os
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description="Extract genes with specific annotation signature")
parser.add_argument("input_files", nargs="+", help="CSV files with input data")
parser.add_argument("--sigs", nargs="+", help="Signature(s)")
#parser.add_argument("--ren", action="store_true", help="Rename seqs based on source file")
args = parser.parse_args()

files = args.input_files
genes = set()
ext = []
for f in files:
	#FIND GENES WITH SPECIFIC SIGNATURE
	with open(f, "rb") as f1:
		hits = {}
		for l in f1:
			if len(l.split("\t")) >= 5:
				////////////////////for s in sigs:
					if l.split("\t")[4].strip() in args.sigs:
						genes.add(l.split()[0])
	
	#EXTRACT FOUND GENES
	with open(f[:-4], "rb") as f1:
		for s in SeqIO.parse(f1, "fasta"):
			if s.id in genes:
				s2 = s
				s.id = os.path.splitext(f)[0]
				ext.append(s2)
	
SeqIO.write(ext, "sig_"+args.sig.replace(":","_")+"_ext.fas", "fasta")

