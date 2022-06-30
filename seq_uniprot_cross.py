#!/usr/bin/env python

from Bio import SeqIO
import sys, argparse
from collections import defaultdict 

parser = argparse.ArgumentParser(description="Crossref Uniprot sequences and output only those present in all species")
parser.add_argument('f1', nargs="+", help="Input FASTA file")
parser.add_argument('--os', action="store_true", default=False, help="")
args = parser.parse_args()

species = set()
print args.f1[0]
s = "Tax="
e = "RepID="
for seq_record in SeqIO.parse(args.f1[0], "fasta"):
	d = seq_record.description
	if args.os == True:
		s = "OS="
		e = "GN="
	if d.find(s) != -1:
		sp = d[d.find(s)+len(s):d.find(e)]
		if len(sp.split()) >= 2:
			sp = sp.split()[0]+"_"+sp.split()[1]
			if sp not in species:
				species.add(sp)
				#~ print sp
print len(species)
for f in args.f1[1:]:
	print
	print f
	species_list = set()
	for seq_record in SeqIO.parse(f, "fasta"):
		d = seq_record.description
		if args.os == True:
			s = "OS="
			e = "GN="
		if d.find(s) != -1:
			sp = d[d.find(s)+len(s):d.find(e)]
			if len(sp.split()) >= 2:
				sp = sp.split()[0]+"_"+sp.split()[1]
				if sp not in species_list:
					species_list.add(sp)
					#~ print sp
	species = species.intersection(species_list)
	print len(species)
	#~ break
	
print "\n===INTERSECTION==="
for s in species:
	print s
	

#~ out = defaultdict(list)

#~ SeqIO.write(seqs, filename, "fasta")
