#!/usr/bin/env python

import os, sys, argparse
from collections import defaultdict
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Obtain distribution of KEGG modules')
parser.add_argument("f", nargs="+", help="Files with input data.")
parser.add_argument("--inputlist", action="store_true", default=False, help="Input file contains file list")
parser.add_argument("--kegglist", default="", help="KEGG list file")
parser.add_argument("--kegg", nargs="+", help="KEGG to extract")

args = parser.parse_args()

files = []
if args.inputlist == True:
	with open(args.f[0]) as inputfile:
		for l in inputfile:
			files.append(l.strip().split()[0])
else:		
	files = args.f

seq_ids = defaultdict(list)
seq_out = defaultdict(list)

keggs = []
if args.kegglist != "":
	with open(args.kegglist) as f:
		for l in f:
			keggs.append(l.strip())
else:
	keggs = args.kegg

print "Extracting keggs"," ".join(keggs)
sys.stdout.flush()

files_shortlist = set()
for i,f in enumerate(files):
	with open(f) as f1:
		print(str(i+1)+"/"+str(len(files))+" "+f)
		sys.stdout.flush()
		for l in f1:
			l = l.strip()
			for k in keggs:
				if len(l.split()) >= 2 and l.split()[1] == k:
					seq_ids[k].append([f,l.split()[0]])
					files_shortlist.add(f)

for i,f in enumerate(files_shortlist):
	seqs = SeqIO.parse(f.replace(".txt",".fas"),"fasta")
	print(str(i+1)+"/"+str(len(files_shortlist))+" "+f.replace(".txt",".fas"))
	sys.stdout.flush()
	for s in seqs:
		for k in keggs:
			for ks in seq_ids[k]:
				if ks[0] == f and s.id == ks[1]:
					s.id = f.split("/")[-1]+"_"+s.id
					seq_out[k].append(s)
					break

for k in seq_out.keys():
	SeqIO.write(seq_out[k],k+"_seqs.fas","fasta")
