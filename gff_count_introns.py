#!/usr/bin/env python

import os, sys, argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='')
parser.add_argument('dir', help="Dir.")
args = parser.parse_args()

for file1 in sorted(os.listdir(args.dir)):
	with open(args.dir+"/"+file1, "rb") as f:
		genes = defaultdict(list)
		for l in f:
			l = l.strip()
			contig = l.split()[0]
			cds = l.split()[2]
			start = l.split()[3]
			end = l.split()[4]
			gene = l.split()[-1]
			if cds == "CDS":
				genes[gene].append(start+"_"+end)
		
		counts = defaultdict(int)		
		for g in genes:
			counts[len(genes[g])] += 1
		print file1,
		for c in sorted(counts):
			print counts[c],
		print


