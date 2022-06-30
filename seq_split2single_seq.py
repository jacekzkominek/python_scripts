#!/usr/bin/env python

import warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import sys, close

parser = argparse.ArgumentParser(description="Split alignment into individual files.")
parser.add_argument('f1', help="Input FASTA file")
args = parser.parse_args()

f = open(args.f1)
seqs = list(SeqIO.parse(f, "fasta"))

out_seqs = []

for seq in seqs:
	SeqIO.write(seq, seq.description+".fas", "fasta")
	#SeqIO.write(seq, seq.description.replace("|","_").replace("/","_")+".fas", "fasta")

f.close()

