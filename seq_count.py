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
		if f.find(".fas") != -1:
			fh = open(f)
			seqs = list(SeqIO.parse(fh,"fasta"))
			print f,len(seqs)
