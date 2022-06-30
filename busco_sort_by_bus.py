#!/usr/bin/env python

import warnings, argparse, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import sys, close
from collections import defaultdict

parser = argparse.ArgumentParser(description=".")
parser.add_argument('flist', nargs='+', help="Input FASTA file")
args = parser.parse_args()

files = []
if args.flist == ["all"]:
	for f in os.listdir("."):
		#if f.find("full.pep") != -1:
		if f.find("full.pep_dedup") != -1:
			files.append(f)
else:
	files = args.flist
	
bus = defaultdict(list)
for f in files:
	seqs = SeqIO.parse(f, "fasta")
	print f
	for s in seqs:
		s.description=f
		s.name=s.name.strip()
		bus[s.name].append(s)	

os.makedirs(os.getcwd()+"/bus_seqs")
for b in bus:
	SeqIO.write(bus[b], "./bus_seqs/"+str(len(bus[b]))+"_"+b+".fas", "fasta")
