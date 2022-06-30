#!/usr/bin/env python

import os, sys, argparse
from Bio import SeqIO
from collections import defaultdict

parser = argparse.ArgumentParser(description="Inspect results of InterProScan.")
parser.add_argument("blast_file", help="Input file")
parser.add_argument("seq_file", help="Input file with refseqs")
parser.add_argument("dir", help="Input dir with seqs")
parser.add_argument("--eval", default="1E-50", help="Threshold evalue")
parser.add_argument("--extract", nargs="?", help="Refhits to extract")
parser.add_argument("--extract_all", default=False, action="store_true", help="Extract all refhits")
args = parser.parse_args()

seq_ids = []
with open(args.seq_file) as f:
	for s in SeqIO.parse(f,"fasta"):
		seq_ids.append(s.id)
	
#blast = []
genes = defaultdict(lambda: defaultdict(list))
with open(args.blast_file) as f:
	source = ""
	for l in f:
		if len(l.split()) == 1:
			source = l.strip()
		else:
			query = l.split()[0]
			target = l.split()[1]
			ev = float(l.split()[10])
			bit = float(l.split()[11])
			if ev >= float(args.eval):
				genes[query][source].append(target)

out_seqs = defaultdict(list)
for f in sorted(os.listdir(args.dir)):
	for g in genes.keys():
		if f in genes[g].keys():
			with open(args.dir+"/"+f) as f1:
				for s in SeqIO.parse(f1,"fasta"):
					if s.id in genes[g][f]:
						s2 = s
						s2.description = ""
						s2.id = ""
						s2.name = ""
						s2.id = f+"_"+s.id.replace(" ","_")
						out_seqs[g].append(s2)

for g in out_seqs:
	SeqIO.write(out_seqs[g], g+"_out.fas", "fasta")
