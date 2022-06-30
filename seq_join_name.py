#!/usr/bin/env python

import os
import sys
import argparse
import fnmatch
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Join all seqs in current directory')
parser.add_argument('fasta_file', default="", nargs='?', help="Fasta file with sequence data.")
parser.add_argument('ext', default="fas", nargs='?', help="")
parser.add_argument('select', default="fas", nargs='?', help="")
args = parser.parse_args()

seqs = []

if args.fasta_file == "all":
	for fname in os.listdir("."):
		if os.path.isfile(fname):
			if fnmatch.fnmatch (fname, "*."+args.ext):
				for seq in SeqIO.parse(fname, "fasta"):
					if seq.id == args.select:
						seq.id = id=fname[:-len("."+args.ext)]
						seq.description=""
						seqs.append(seq)
						break

SeqIO.write(seqs, "joined.fas", "fasta")
