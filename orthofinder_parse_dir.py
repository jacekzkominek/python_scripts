#!/usr/bin/env python

import warnings
import sys, argparse, os, subprocess
from collections import defaultdict 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="")
parser.add_argument("query", help="Query FASTA")
parser.add_argument("dir", help="Target dir")
parser.add_argument("--fasta_dir", default="", help="Dir with FASTA files")
args = parser.parse_args()

queries = []
with open(args.query) as q:
	for l in q:
		if l.find(">") != -1:
			queries.append(l.split()[0][1:])

print "Input\t"+"\t".join(sorted(queries))#+"\tTotal"			

q_map = defaultdict(lambda : defaultdict(list))
for d in sorted(os.listdir(args.dir)):
	for f in sorted(os.listdir(args.dir+"/"+d)):
		if f.find("Results") != -1 and os.path.exists(args.dir+"/"+d+"/"+f+"/Orthogroups.csv"):
			with open(args.dir+"/"+d+"/"+f+"/Orthogroups.csv","r") as f1:
				for l in f1:
					if l.split()[0].find("OG") != -1:
						refs = l.split("\t")[1].strip()
						if len(refs.split(",")) == 1:
							for q in queries:
								#~ print refs,q
								if refs == q:
									#~ print "AA"
									q_map[d][q] = len(l.split("\t")[2].split(","))
						#~ else:
							#~ for r in refs.split(","):
								#~ r = r.strip()
								#~ first = True
								#~ for q in queries:
									#~ if r == q:
										#~ if first == True:
											#~ q_map[d][q] == 1
										#~ else:
											#~ q_map[d][q] == 0
		
for sp in sorted(q_map):
	print sp,
	for q in sorted(queries):
		if q not in q_map[sp]:
			print 0,
		else:
			print q_map[sp][q],
	print
