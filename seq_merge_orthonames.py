#!/usr/bin/env python

import os, sys, warnings, argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='Merge sequences from multiple FASTA by name.')
parser.add_argument('dir', help="Dir with FASTA files")
parser.add_argument('--list', nargs="*", help="List of seq types")
args = parser.parse_args()

seq_dict = defaultdict(list)
for f in os.listdir(args.dir):
	seqs = list(SeqIO.parse(args.dir+"/"+f, "fasta"))
	for s in seqs:
		for lseq in args.list:
			if s.description.find(lseq) != -1:
				seq_dict[lseq].append(s)
				break

for lseq in args.list:
	SeqIO.write(seq_dict[lseq], lseq+".fas", "fasta")

