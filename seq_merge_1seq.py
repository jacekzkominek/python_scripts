#!/usr/bin/env python

import os, sys, warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='Merge sequences from multiple FASTA files into one.')
parser.add_argument('file', help="Subject FASTA files")
args = parser.parse_args()

seqs = list(SeqIO.parse(os.getcwd()+"/"+args.file, "fasta"))
seq_out = ""
for s in seqs:
	seq_out += s.seq

with open(args.file+"_out","w") as f:
	f.writelines(">"+args.file+"_merged\n")
	f.writelines(seq_out+"\n")

