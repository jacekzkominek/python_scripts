#!/usr/bin/env python

import warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import sys, close

parser = argparse.ArgumentParser(description="Reverse order of sequences.")
parser.add_argument('f1', help="Input FASTA file")
args = parser.parse_args()

f = open(args.f1)
seqs = list(SeqIO.parse(f, "fasta"))

out_seqs = [] 
for seq in seqs:
	out_seqs.insert(0,seq)

SeqIO.write(out_seqs, args.f1.replace(".fas", "_rev.fas"), "fasta")

f.close()

