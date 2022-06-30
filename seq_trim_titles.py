#!/usr/bin/env python

import warnings, argparse, os, sys, copy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="Trim sequences with identical names.")
parser.add_argument('f1', help="Input FASTA file")
parser.add_argument('--trim', default="50", help="Length to trim")
args = parser.parse_args()

out_seqs = []

with open(args.f1) as f:
	for seq in list(SeqIO.parse(f, "fasta")):
		seq2 = copy.deepcopy(seq)
		seq2.title = ""
		seq2.id = ""
		seq2.description = seq2.description[:int(args.trim)]
		out_seqs.append(seq2)

SeqIO.write(out_seqs, args.f1.replace(".fas","_trim"+args.trim+".fas"), "fasta")
