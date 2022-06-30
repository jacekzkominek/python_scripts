#!/usr/bin/env python

import warnings
import sys, argparse, os, subprocess
from collections import defaultdict

parser = argparse.ArgumentParser(description="")
#~ parser.add_argument("query", help="Query FASTA")
parser.add_argument("dir", help="Target dir")
#~ parser.add_argument("--tblastn", default=False, action="store_true", help="Use TBLASTN instead of BLASTP")
#~ parser.add_argument("--blastn", default=False, action="store_true", help="Use BLASTN instead of BLASTP")
#~ parser.add_argument("--eval", default="1E-5", help="Evalue")
args = parser.parse_args()

trnas = defaultdict(lambda : defaultdict(int))
trnas_list = set()
for f1 in sorted(os.listdir(args.dir)):
	if f1.find(".gff") != -1:

		with open(args.dir+"/"+f1, "r") as f:
			for l in f:
				gene = l.split()[2]
				if gene == "gene":
					trna = l.split()[-1][l.split()[-1].find("noncoding-")+10:l.split()[-1].find("-gene-")]
					trnas[f1][trna] += 1
					trnas_list.add(trna)

print "",
for trna in sorted(trnas_list):
	print trna,
print 
for f in sorted(trnas.keys()):
	print f,
	for trna in sorted(trnas_list):
		print trnas[f][trna],
	print
