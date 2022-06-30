#!/usr/bin/env python

import warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from os import sys, close

parser = argparse.ArgumentParser(description='Merge sequences from multiple FASTA files into one.')
parser.add_argument('files', nargs='*', help="Subject FASTA files")
parser.add_argument('--out', nargs='?', default="merged_seq.fas", help="Output file")
args = parser.parse_args()

files = []
for f in args.files:
	files.append(SeqIO.parse(f, "fasta"))
	
appendix = []
for f in files:
	for s in f:
		appendix.append(s)
		
SeqIO.write(appendix, args.out, "fasta")


