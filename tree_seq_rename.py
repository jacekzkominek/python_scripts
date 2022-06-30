#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
import os, sys, argparse
from collections import defaultdict 

parser = argparse.ArgumentParser(description="Crossref Uniprot sequences and output only those present in all species")
parser.add_argument('seq', help="Input FASTA file")
parser.add_argument('tree', help="Input NWK file")
parser.add_argument('out_tree', help="Output NWK file")
args = parser.parse_args()

seqs = []

for seq in SeqIO.parse(args.seq, "fasta"):
	seqs.append(seq.description)

check = defaultdict(bool)
for s in seqs:
	check[s] = False
	
out = []
with open(args.tree) as f:
	for l in f:
		l = l.strip()
		found = False
		if l.find("'") != -1:
			sp = l[l.find("'")+1:l.rfind("'")]
			for s in seqs:
				if s.find(sp.replace(" ","_")) != -1:
					out.append(l.replace(sp,s))
					found = True
					check[s] = True
					break
			if found == False:
				out.append(l)
				print l,"- on tree, not found in seq!"
		else:
			out.append(l)

for c in check:
	if check[c] == False:
		print c,"- in seq, not on tree"

outf = open(args.out_tree,"w")
for l in out:
	outf.write(l+"\n")
	
