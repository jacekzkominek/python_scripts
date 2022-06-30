#!/usr/bin/env python

import os, warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import sys, close

parser = argparse.ArgumentParser(description="Remove the trailing asterisk from amino acid sequences.")
parser.add_argument('f1', help="Input FASTA file")
args = parser.parse_args()

f1 = open(args.f1)
seqs = list(SeqIO.parse(f1, "fasta"))
out_seqs = []
for seq in seqs:
	print ">"+seq.description+"\n"+str(seq.seq).replace("-","")
f1.close()

