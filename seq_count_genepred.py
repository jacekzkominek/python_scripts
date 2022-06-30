#!/usr/bin/env python

import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, argparse, os

parser = argparse.ArgumentParser(description="Filter stop codons (\"*\") from protein sequences.")
parser.add_argument('file1', help="Input file")
args = parser.parse_args()

if args.file1 == "all":
	for f in sorted(os.listdir(os.getcwd())):
		snap = 0
		aug = 0
		gm = 0
		m = 0
		if f.find(".fas") != -1:
			fh = open(f)
			seqs = list(SeqIO.parse(fh,"fasta"))
			for s in seqs:
				if s.id.find("snap") != -1:
					snap += 1
				elif s.id.find("augustus") != -1:
					aug += 1
				elif s.id.find("genemark") != -1:
					gm += 1
				elif s.id.find("maker") != -1:
					m += 1
				else:
					print s.id
			print f,len(seqs),str(snap),str(aug),str(gm),str(m)
