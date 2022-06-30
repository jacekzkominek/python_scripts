#!/usr/bin/env python

import warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import sys, close

parser = argparse.ArgumentParser(description="Sort sequences based on the ortholog type.")
parser.add_argument('f1', help="Input FASTA file")
args = parser.parse_args()

f = open(args.f1)
seqs = list(SeqIO.parse(f, "fasta"))
out_seqs = dict()

for seq in seqs:
	orto = seq.description[4:seq.description.find("|")]
	if (orto in out_seqs.keys()):
		out_seqs[orto].append(seq)
	else:
		new_orto = []
		new_orto.append(seq)
		out_seqs[orto] = new_orto
		
for orto in out_seqs.keys():
	SeqIO.write(out_seqs[orto], orto+".fas", "fasta")

f.close()


