#!/usr/bin/env python

import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, argparse, os
from collections import defaultdict

parser = argparse.ArgumentParser(description="Filter stop codons (\"*\") from protein sequences.")
parser.add_argument('files', nargs="*", help="Input file")
parser.add_argument('--genes', nargs="*", help="Input file")
args = parser.parse_args()

genes = defaultdict(list)
sp = set()
for f in args.files:
	for s in list(SeqIO.parse(open(f),"fasta")):
		target = ""
		if s.description.find("yHMPu") == -1:
			target = str(s.description.split("_")[1]+"_"+s.description.split("_")[2])
		else:
			target = str(s.description.split("_")[1]+"_"+s.description.split("_")[2]+"_"+s.description.split("_")[3])
		genes[s.description.split("_")[0]].append(target)
		sp.add(target)

print "\t"
for s in sorted(sp):
	print "\t"+s,
print

for g in genes:
	print g,
	for s in sorted(sp):
		if s in genes[g]:
			print "\t"+"1",
		elif s not in genes[g]:
			print "\t"+"0",
	print
		
