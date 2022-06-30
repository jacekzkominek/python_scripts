#!/usr/bin/env python

from Bio import SeqIO
import sys, argparse
from collections import defaultdict 

parser = argparse.ArgumentParser(description="Crossref Uniprot sequences and output only those present in all species")
parser.add_argument('f1', nargs="+", help="Input FASTA file")
parser.add_argument('--os', action="store_true", default=False, help="")
args = parser.parse_args()

species = set()
multi = list()
s = "Tax="
e = "RepID="
count = defaultdict(int)
for seq_record in SeqIO.parse(args.f1[0], "fasta"):
	d = seq_record.description
	if args.os == True:
		s = "OS="
		e = "GN="
	if d.find(s) != -1:
		sp = d[d.find(s)+len(s):d.find(e)]
		if len(sp.split()) >= 2 and sp.find("sp.") == -1:
			count[sp] += 1
			if sp in species:
				species.discard(sp)
				#~ multi.add(sp)
			if sp.find("(strain") == -1:
				species.add(sp)
			if sp.find("(strain") != -1:
				species.add(sp.split()[0]+" "+sp.split()[1])
			#~ sp = sp.split()[0]+"_"+sp.split()[1]
			#~ if sp not in species:
				#~ species.add(sp)
		
for s in sorted(species):
	print s
	
#~ for s in count:
	#~ print count[s]
