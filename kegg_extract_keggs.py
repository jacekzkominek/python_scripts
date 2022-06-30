#!/usr/bin/env python

import os, sys, argparse
from collections import defaultdict
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Obtain distribution of KEGG modules')
parser.add_argument("f", help="Files with input data.")
parser.add_argument("--kegglist", default="", help="KEGG list file")
parser.add_argument("--kegg", nargs="+", help="KEGG to extract")

args = parser.parse_args()

keggs = []
if args.kegglist != "":
	with open(args.kegglist) as f:
		for l in f:
			l = l.strip()
			for k1 in l.split():
				keggs.append(k1)
else:
	keggs = args.kegg

keggs_out = []
species = []
with open(args.f) as kf:
	for l in kf:
		l = l.strip()
		k = l.split(",")[0]
		if k == "KEGG":
			species = l.split(",")[1:]
		else:
			if k in keggs:
				keggs_out.append(l.split(","))

print("\t".join(["Species"]+[k[0] for k in keggs_out]))

for i,s in enumerate(species):
	print("\t".join([s]+[k[i+1] for k in keggs_out]))
