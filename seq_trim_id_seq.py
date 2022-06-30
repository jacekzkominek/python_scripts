#!/usr/bin/env python

import warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import sys, close
from collections import defaultdict

parser = argparse.ArgumentParser(description="Trim sequences from identical organisms (OS species annotation).")
parser.add_argument('f1', help="Input FASTA file")
args = parser.parse_args()

f = open(args.f1)
seqs = list(SeqIO.parse(f, "fasta"))

out_seqs = []
out_dict = defaultdict(list)

for seq in seqs:
	title = seq.description
	title= title[title.find("OS=")+3:]
	t1 = title.split(" ")[0]+" "+title.split(" ")[1]
	if t1 not in out_dict.keys():
		out_seqs.append(seq)
		out_dict[t1].append(str(seq.seq))
	else:
		if str(seq.seq) not in out_dict[t1]:
			out_seqs.append(seq)
			out_dict[t1].append(str(seq.seq))

#print len(out_seqs)
SeqIO.write(out_seqs, args.f1.replace(".fas", "_trim_id_seq.fas"), "fasta")

f.close()

